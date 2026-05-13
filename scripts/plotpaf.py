#!/usr/bin/env python3
"""Plot backbone-vs-reference PAF alignment.

Python port of plotpaf_png.R with identical chaining logic and polished visuals.
Only dependency: matplotlib.

Usage: plotpaf.py <input.paf> <output_prefix> [options]
  --reference <genome.fa>   Reference FASTA — auto-generates position map via physlr
  --positions <pos.tsv>     Pre-computed position map (from physlr index-positions)
  --caption                 Show input file path at bottom of each plot

Produces: <output_prefix>.backbone.png and <output_prefix>.reference.png
When --reference or --positions is provided, both plots use genomic bp coordinates.
"""
import os
import sys
from collections import Counter, defaultdict
from bisect import bisect_right
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Polychrome alphabet.colors(26) with iron → indigo fix
# ---------------------------------------------------------------------------
ALPHABET_COLOURS = [
    '#AA0DFE', '#3283FE', '#85660D', '#782AB6', '#565656',
    '#1C8356', '#16FF32', '#F7E1A0', '#3F00FF', '#1CBE4F',
    '#C4451C', '#DEA0FD', '#FE00FA', '#325A9B', '#FEAF16',
    '#F8A19F', '#90AD1C', '#F6222E', '#1CFFCE', '#2ED9FF',
    '#B10DA1', '#C075A6', '#FC1CBF', '#B00068', '#FBE426',
    '#FA0087',
]

ALLOWED_CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
ALLOWED_SET = set(ALLOWED_CHROMS)

# GRCh38 chromosome sizes (bp)
GRCH38_SIZE = {
    'chr1': 248_956_422, 'chr2': 242_193_529, 'chr3': 198_295_559,
    'chr4': 190_214_555, 'chr5': 181_538_259, 'chr6': 170_805_979,
    'chr7': 159_345_973, 'chr8': 145_138_636, 'chr9': 138_394_717,
    'chr10': 133_797_422, 'chr11': 135_086_622, 'chr12': 133_275_309,
    'chr13': 114_364_328, 'chr14': 107_043_718, 'chr15': 101_991_189,
    'chr16': 90_338_345,  'chr17': 83_257_441,  'chr18': 80_373_285,
    'chr19': 58_617_616,  'chr20': 64_444_167,  'chr21': 46_709_983,
    'chr22': 50_818_468,  'chrX': 156_040_895,  'chrY': 57_227_415,
}

# GRCh38 centromere ranges (start, end) from UCSC centromeres.txt
CENTROMERE_RANGE = {
    'chr1':  (121_700_000, 125_100_000), 'chr2':  (91_800_000,  96_000_000),
    'chr3':  (87_800_000,  94_000_000),  'chr4':  (48_200_000,  51_800_000),
    'chr5':  (46_100_000,  51_400_000),  'chr6':  (58_500_000,  62_600_000),
    'chr7':  (58_100_000,  62_100_000),  'chr8':  (43_200_000,  47_200_000),
    'chr9':  (42_200_000,  45_500_000),  'chr10': (38_000_000,  41_600_000),
    'chr11': (51_000_000,  55_800_000),  'chr12': (33_200_000,  37_800_000),
    'chr13': (16_500_000,  18_900_000),  'chr14': (16_100_000,  18_200_000),
    'chr15': (17_500_000,  20_500_000),  'chr16': (35_300_000,  38_400_000),
    'chr17': (22_700_000,  27_400_000),  'chr18': (15_400_000,  21_500_000),
    'chr19': (24_200_000,  28_100_000),  'chr20': (26_200_000,  30_400_000),
    'chr21': (10_900_000,  13_000_000),  'chr22': (13_700_000,  17_400_000),
    'chrX':  (58_100_000,  63_800_000),  'chrY':  (10_300_000,  10_600_000),
}

# ---------------------------------------------------------------------------
# Position map: minimizer index → genomic bp
# ---------------------------------------------------------------------------

def load_position_map(path):
    """Load a position map TSV from `physlr index-positions`.

    Returns dict: chrom -> (mx_indices, bp_positions, total_mx, seq_len)
    where mx_indices and bp_positions are sorted lists for interpolation.
    """
    data = defaultdict(lambda: ([], [], [0], [0]))
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5 or parts[0] == 'chrom':
                continue
            chrom = parts[0]
            mx_idx = int(parts[1])
            bp_pos = int(parts[2])
            total_mx = int(parts[3])
            seq_len = int(parts[4])
            entry = data[chrom]
            entry[0].append(mx_idx)
            entry[1].append(bp_pos)
            entry[2][0] = total_mx
            entry[3][0] = seq_len
    result = {}
    for chrom, (mx_list, bp_list, total_list, seqlen_list) in data.items():
        result[chrom] = (mx_list, bp_list, total_list[0], seqlen_list[0])
    return result


def mx_to_bp(mx_val, mx_indices, bp_positions, paf_qlength=None, map_total=None):
    """Interpolate a minimizer index to genomic bp using the position map.

    If paf_qlength differs from map_total (due to code version differences),
    scales mx_val to the map's coordinate space before interpolating.
    """
    if paf_qlength and map_total and paf_qlength != map_total:
        mx_val = mx_val * (map_total / paf_qlength)

    if not mx_indices:
        return mx_val
    if mx_val <= mx_indices[0]:
        return bp_positions[0]
    if mx_val >= mx_indices[-1]:
        return bp_positions[-1]

    idx = bisect_right(mx_indices, mx_val) - 1
    if idx >= len(mx_indices) - 1:
        return bp_positions[-1]

    m0, m1 = mx_indices[idx], mx_indices[idx + 1]
    b0, b1 = bp_positions[idx], bp_positions[idx + 1]
    if m1 == m0:
        return b0
    frac = (mx_val - m0) / (m1 - m0)
    return b0 + frac * (b1 - b0)


# ---------------------------------------------------------------------------
# PAF I/O
# ---------------------------------------------------------------------------

def read_paf(path):
    records = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            records.append({
                'Qname':  parts[0],
                'Qlength': int(parts[1]),
                'Qstart':  int(parts[2]),
                'Qend':    int(parts[3]),
                'Orientation': parts[4],
                'Tname':  parts[5],
                'Tlength': int(parts[6]),
                'Tstart':  int(parts[7]),
                'Tend':    int(parts[8]),
                'Matches': int(parts[9]),
                'Length':  int(parts[10]),
                'Mapq':   int(parts[11]),
            })
    return records

# ---------------------------------------------------------------------------
# Chaining (faithful to R script)
# ---------------------------------------------------------------------------

def _overlap(a_s, a_e, b_s, b_e):
    return (a_s <= b_s <= a_e or a_s <= b_e <= a_e or
            b_s <= a_s <= b_e or b_s <= a_e <= b_e)


def _interval_dist(a_s, a_e, b_s, b_e):
    if _overlap(a_s, a_e, b_s, b_e):
        return 0
    return max(b_s - a_e, a_s - b_e)


def chain_alignments(records):
    MAPQ_THR = 1
    NBARCODES_THR = 50
    QDIST_THR = 500_000
    TDIST_THR = 10

    filtered = [r for r in records if r['Mapq'] >= MAPQ_THR]
    filtered.sort(key=lambda r: (-r['Tlength'], r['Tname'], r['Tstart']))

    groups = defaultdict(list)
    for r in filtered:
        groups[(r['Qname'], r['Tname'])].append(r)

    all_segments = []

    for (qname, tname), grp in groups.items():
        grp.sort(key=lambda r: r['Tstart'])
        n = len(grp)
        if n < 2:
            continue

        keep = [False] * n
        for i in range(n):
            if i > 0 and _overlap(grp[i-1]['Qstart'], grp[i-1]['Qend'],
                                   grp[i]['Qstart'], grp[i]['Qend']):
                keep[i] = True
            if i < n - 1 and _overlap(grp[i]['Qstart'], grp[i]['Qend'],
                                       grp[i+1]['Qstart'], grp[i+1]['Qend']):
                keep[i] = True
        grp = [r for r, k in zip(grp, keep) if k]
        if not grp:
            continue

        seg_id = 0
        segments_map = defaultdict(list)
        segments_map[0].append(grp[0])
        for i in range(1, len(grp)):
            prev, cur = grp[i-1], grp[i]
            qd = _interval_dist(prev['Qstart'], prev['Qend'],
                                cur['Qstart'], cur['Qend'])
            td = _interval_dist(prev['Tstart'], prev['Tend'],
                                cur['Tstart'], cur['Tend'])
            if qd >= QDIST_THR or td >= TDIST_THR:
                seg_id += 1
            segments_map[seg_id].append(cur)

        for seg in segments_map.values():
            orient_counts = Counter(r['Orientation'] for r in seg)
            all_segments.append({
                'Qname': qname,
                'Qlength': seg[0]['Qlength'],
                'Qstart': min(r['Qstart'] for r in seg),
                'Qend':   max(r['Qend']   for r in seg),
                'Orientation': orient_counts.most_common(1)[0][0],
                'Tname': tname,
                'Tlength': seg[0]['Tlength'],
                'Tstart': min(r['Tstart'] for r in seg),
                'Tend':   max(r['Tend']   for r in seg),
                'Matches': sum(r['Matches'] for r in seg),
                'Length':  sum(r['Length']  for r in seg),
                'Mapq':   sorted(r['Mapq'] for r in seg)[len(seg)//2],
                'Barcodes': len(seg),
            })

    print(f"Chained segments before filter: {len(all_segments)}")
    chained = [s for s in all_segments if s['Barcodes'] >= NBARCODES_THR]
    print(f"Chained segments (min barcodes = {NBARCODES_THR}): {len(chained)}")

    if not chained:
        NBARCODES_THR = 10
        chained = [s for s in all_segments if s['Barcodes'] >= NBARCODES_THR]
        print(f"Fallback (min barcodes = {NBARCODES_THR}): {len(chained)}")

    return chained

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def chrom_sort_key(name):
    s = name.replace('chr', '')
    if s == 'X': return 23
    if s == 'Y': return 24
    try: return int(s)
    except ValueError: return 99


def _apply_theme(ax, fig):
    """Clean theme matching ggplot2 theme_bw()."""
    fig.set_facecolor('white')
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#404040')
        spine.set_linewidth(0.5)
    ax.grid(True, which='major', axis='both', color='#E8E8E8',
            linewidth=0.4, linestyle='-')
    ax.set_axisbelow(True)
    ax.tick_params(axis='both', which='major', direction='out',
                   length=4, width=0.5, color='#404040', labelsize=10, pad=4)
    ax.tick_params(axis='both', which='minor', length=0)


def _make_legend(ax, handles, labels, title, fontsize=9, expand=False):
    """ggplot2-style legend: no frame, fills vertical space, large swatches.
    If expand=True, use R-style large legend filling ~90% of plot height."""
    if expand:
        n = len(labels)
        # R-style: large square swatches, big font, filling the plot height
        fs = 11
        title_fs = 13
        spacing = 0.6
        hheight = 1.6
        hlength = 2.2
        htextpad = 0.8
    else:
        fs = fontsize
        title_fs = 11
        spacing = 0.25
        hheight = 1.2
        hlength = 1.8
        htextpad = 0.6
    leg = ax.legend(
        handles, labels,
        title=title,
        loc='center left',
        bbox_to_anchor=(1.01, 0.5),
        fontsize=fs,
        title_fontsize=title_fs,
        ncol=1,
        frameon=False,
        borderpad=0.3,
        labelspacing=spacing,
        handlelength=hlength,
        handleheight=hheight,
        handletextpad=htextpad,
        borderaxespad=0.3,
    )
    leg.get_title().set_fontweight('bold')
    return leg

# ---------------------------------------------------------------------------
# Plot 1: Backbone coverage
# ---------------------------------------------------------------------------

def plot_backbone(chained, output_prefix, show_caption, input_paf, pos_map=None):
    # Collect all backbone names and their molecule counts
    tinfo = {}
    for s in chained:
        t = s['Tname']
        if t not in tinfo or s['Tlength'] > tinfo[t]:
            tinfo[t] = s['Tlength']
    all_tnames = list(tinfo.keys())

    # Qlength from PAF (for scaling in mx_to_bp)
    qlens_mx = {}
    for s in chained:
        q = s['Qname']
        qlens_mx[q] = max(qlens_mx.get(q, 0), s['Qlength'])

    use_genomic = pos_map is not None

    def to_bp(qname, mx_val):
        """Convert minimizer index to bp using position map, or return as-is."""
        if not use_genomic or qname not in pos_map:
            return mx_val
        mx_indices, bp_positions, map_total, seq_len = pos_map[qname]
        return mx_to_bp(mx_val, mx_indices, bp_positions,
                        paf_qlength=qlens_mx.get(qname), map_total=map_total)

    # Build bp coordinate mapping for each backbone.
    # Each chained segment maps a reference region (Qstart-Qend) to a backbone
    # region (Tstart-Tend). We convert the reference span to bp and lay out
    # segments sequentially along the backbone in bp space.
    backbone_bp = {}  # tname -> {(Tstart,Tend,Qname) -> (bp_start, bp_end)}
    backbone_len_bp = {}  # tname -> total bp length
    for tname in all_tnames:
        segs = sorted(
            [s for s in chained if s['Tname'] == tname],
            key=lambda s: s['Tstart'])
        cum_bp = 0.0
        mapping = {}
        for s in segs:
            q_bp_start = to_bp(s['Qname'], s['Qstart'])
            q_bp_end = to_bp(s['Qname'], s['Qend'])
            seg_bp = abs(q_bp_end - q_bp_start)
            if seg_bp < 1:
                seg_bp = 1
            mapping[(s['Tstart'], s['Tend'], s['Qname'])] = (cum_bp, cum_bp + seg_bp)
            cum_bp += seg_bp
        backbone_bp[tname] = mapping
        backbone_len_bp[tname] = cum_bp if cum_bp > 0 else tinfo[tname]

    # Sort backbones: by estimated bp length (descending) when available,
    # otherwise by molecule count
    if use_genomic:
        tnames = sorted(all_tnames, key=lambda t: (-backbone_len_bp[t], t))
    else:
        tnames = sorted(all_tnames, key=lambda t: (-tinfo[t], t))
    tname_idx = {t: i for i, t in enumerate(tnames)}
    n_paths = len(tnames)

    # Chromosome colors
    qnames = sorted({s['Qname'] for s in chained}, key=chrom_sort_key)
    qcolor = {q: ALPHABET_COLOURS[i % len(ALPHABET_COLOURS)]
              for i, q in enumerate(qnames)}

    fig, ax = plt.subplots(figsize=(13, 10))
    _apply_theme(ax, fig)

    for s in chained:
        ti = tname_idx[s['Tname']]
        key = (s['Tstart'], s['Tend'], s['Qname'])
        if use_genomic and key in backbone_bp.get(s['Tname'], {}):
            x0, x1 = backbone_bp[s['Tname']][key]
        else:
            x0, x1 = s['Tstart'], s['Tend']
        ax.add_patch(Rectangle(
            (x0, ti), x1 - x0, 1,
            facecolor=qcolor.get(s['Qname'], '#999'), edgecolor='none'))

    if use_genomic:
        max_x = max(backbone_len_bp.values()) * 1.01
    else:
        max_x = max(s['Tend'] for s in chained) * 1.01
    ax.set_xlim(0, max_x)
    ax.set_ylim(n_paths, 0)

    if use_genomic:
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f'{x / 1e6:.0f}Mb'))
        ax.set_xlabel('Estimated Size (bp)', fontsize=11, labelpad=7)
    else:
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f'{int(x):,}' if x >= 1 else '0'))
        ax.set_xlabel('Position (minimizer index)', fontsize=11, labelpad=7)

    # Y-axis
    ax.set_ylabel('Backbone (Target)', fontsize=11, labelpad=7)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=12))

    # Legend
    handles = [Rectangle((0, 0), 1, 1, facecolor=qcolor[q]) for q in qnames]
    labels = list(qnames)
    _make_legend(ax, handles, labels, 'Chromosome', expand=True)

    if show_caption:
        fig.text(0.45, 0.004, input_paf, ha='center', fontsize=5,
                 color='#AAAAAA', style='italic')

    plt.tight_layout(rect=[0, 0.01, 0.85, 1])
    out = f'{output_prefix}.backbone.png'
    plt.savefig(out, dpi=150, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close()
    print(f"Saved: {out}")

# ---------------------------------------------------------------------------
# Plot 2: Reference coverage
# ---------------------------------------------------------------------------

def plot_reference(chained, output_prefix, show_caption, input_paf, pos_map=None):
    qnames = sorted({s['Qname'] for s in chained}, key=chrom_sort_key)
    qidx = {q: i + 1 for i, q in enumerate(qnames)}
    n_chroms = len(qnames)

    # Qlength from PAF = minimizer index length
    qlens_mx = {}
    for s in chained:
        q = s['Qname']
        qlens_mx[q] = max(qlens_mx.get(q, 0), s['Qlength'])

    use_genomic = pos_map is not None

    def to_bp(qname, mx_val):
        """Convert minimizer index to bp using position map, or return as-is."""
        if not use_genomic or qname not in pos_map:
            return mx_val
        mx_indices, bp_positions, map_total, seq_len = pos_map[qname]
        return mx_to_bp(mx_val, mx_indices, bp_positions,
                        paf_qlength=qlens_mx.get(qname), map_total=map_total)

    # Chromosome lengths: real bp if position map available, else minimizer index
    if use_genomic:
        qlens = {}
        for q in qnames:
            if q in pos_map:
                qlens[q] = pos_map[q][3]  # seq_len from position map
            else:
                qlens[q] = GRCH38_SIZE.get(q, qlens_mx.get(q, 0))
    else:
        qlens = dict(qlens_mx)

    # Top 25 backbone paths by total barcodes (labels 1-25)
    tbar = defaultdict(int)
    for s in chained:
        tbar[s['Tname']] += s['Barcodes']
    top25 = sorted(tbar, key=lambda t: -tbar[t])[:25]
    top_set = set(top25)

    tlen_map = {}
    for s in chained:
        if s['Tname'] in top_set:
            tlen_map[s['Tname']] = max(tlen_map.get(s['Tname'], 0), s['Tlength'])
    ordered_t = sorted(top25, key=lambda t: (-tlen_map.get(t, 0), t))
    tcolor = {t: ALPHABET_COLOURS[i % len(ALPHABET_COLOURS)]
              for i, t in enumerate(ordered_t)}

    chained_sorted = sorted(chained,
        key=lambda s: (-s['Qlength'], s['Qname'], s['Qstart'], s['Barcodes']))

    # Row geometry: small gap between rows
    gap = 0.06
    bar_h = 1.0 - gap  # bar height per row

    fig, ax = plt.subplots(figsize=(13, 10))
    _apply_theme(ax, fig)

    # Layer 1: Ideogram — alternating gray
    for qname in qnames:
        qi = qidx[qname]
        y_top = qi + gap / 2
        rows = [s for s in chained_sorted if s['Qname'] == qname]
        for j, s in enumerate(rows):
            gray = '#808080' if (qi + j) % 2 == 0 else '#C0C0C0'
            x0 = to_bp(qname, s['Qstart'])
            x1 = to_bp(qname, s['Qend'])
            ax.add_patch(Rectangle(
                (x0, y_top), x1 - x0, bar_h,
                facecolor=gray, edgecolor='none'))

    # Layer 2: Colored backbone segments
    for s in chained_sorted:
        if s['Tname'] in tcolor:
            qi = qidx[s['Qname']]
            y_top = qi + gap / 2
            x0 = to_bp(s['Qname'], s['Qstart'])
            x1 = to_bp(s['Qname'], s['Qend'])
            ax.add_patch(Rectangle(
                (x0, y_top), x1 - x0, bar_h,
                facecolor=tcolor[s['Tname']], edgecolor='none'))

    # Layer 3: Chromosome length dots
    for qname in qnames:
        qi = qidx[qname]
        ax.plot(qlens[qname], qi + 0.5, 'o',
                color='black', markersize=3, zorder=5)

    # Layer 4: Centromere — hatched rectangle (aligned exactly with bars)
    if use_genomic:
        for qname in qnames:
            cen = CENTROMERE_RANGE.get(qname)
            if cen is None:
                continue
            cen_start, cen_end = cen
            ql = qlens.get(qname, 0)
            if cen_start <= ql:
                qi = qidx[qname]
                y_top = qi + gap / 2
                cen_end_clipped = min(cen_end, ql)
                ax.add_patch(Rectangle(
                    (cen_start, y_top), cen_end_clipped - cen_start, bar_h,
                    facecolor='none', edgecolor='#333333',
                    linewidth=0.3, hatch='///',
                    clip_on=True, zorder=8))

    max_x = max(qlens.values()) * 1.04 if qlens else 100
    ax.set_xlim(0, max_x)

    if use_genomic:
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f'{x / 1e6:.0f}Mb'))
        ax.set_xlabel('Genomic Position (bp)', fontsize=11, labelpad=7)
    else:
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f'{int(x):,}' if x >= 1 else '0'))
        ax.set_xlabel('Position (minimizer index)', fontsize=11, labelpad=7)

    ax.set_ylim(n_chroms + 1, 0.5)
    ax.set_yticks([qidx[q] + 0.5 for q in qnames])
    ax.set_yticklabels(qnames, fontsize=10)
    ax.set_ylabel('Chromosome', fontsize=11, labelpad=7)

    # Legend: backbone 1-25 + Other1/Other2 + Centromere + Chrom end
    handles = [Rectangle((0, 0), 1, 1, facecolor=tcolor[t]) for t in ordered_t]
    labels = [str(i + 1) for i in range(len(ordered_t))]
    handles.append(Rectangle((0, 0), 1, 1, facecolor='#808080'))
    labels.append('Other1')
    handles.append(Rectangle((0, 0), 1, 1, facecolor='#C0C0C0'))
    labels.append('Other2')
    if use_genomic:
        handles.append(Rectangle((0, 0), 1, 1, facecolor='none',
                                 edgecolor='#333333', linewidth=0.5, hatch='///'))
        labels.append('Centromere')
    handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                          markersize=5, markeredgecolor='black'))
    labels.append('Chrom. end')

    _make_legend(ax, handles, labels, 'Backbone', fontsize=8, expand=True)

    if show_caption:
        fig.text(0.45, 0.004, input_paf, ha='center', fontsize=5,
                 color='#AAAAAA', style='italic')

    plt.tight_layout(rect=[0, 0.01, 0.85, 1])
    out = f'{output_prefix}.reference.png'
    plt.savefig(out, dpi=150, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close()
    print(f"Saved: {out}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _find_physlr():
    """Locate the physlr binary. Checks common locations."""
    import shutil
    # Check PATH first
    p = shutil.which('physlr')
    if p:
        return p
    # Check relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for candidate in [
        os.path.join(script_dir, '..', 'target', 'release', 'physlr'),
        os.path.join(script_dir, '..', 'target', 'debug', 'physlr'),
        os.path.join(script_dir, 'physlr'),
    ]:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def _generate_position_map(reference, k=32, w=32, step=100):
    """Run physlr index-positions to generate a position map TSV.
    Returns the path to the generated file, or None on failure."""
    import subprocess
    import tempfile
    physlr = _find_physlr()
    if not physlr:
        print("Warning: physlr binary not found. Cannot generate position map.",
              file=sys.stderr)
        print("  Install physlr or provide --positions <pos.tsv> directly.",
              file=sys.stderr)
        return None
    # Write to a temp file next to the reference
    ref_dir = os.path.dirname(os.path.abspath(reference))
    ref_base = os.path.splitext(os.path.basename(reference))[0]
    pos_path = os.path.join(ref_dir, f'{ref_base}.positions.tsv')
    # Reuse if it already exists
    if os.path.isfile(pos_path) and os.path.getsize(pos_path) > 0:
        print(f"Reusing existing position map: {pos_path}")
        return pos_path
    print(f"Generating position map from {reference}...")
    try:
        subprocess.run(
            [physlr, 'index-positions', reference,
             '-o', pos_path,
             '-k', str(k), '-w', str(w), '--step', str(step)],
            check=True)
        print(f"  Wrote position map to {pos_path}")
        return pos_path
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Warning: failed to generate position map: {e}", file=sys.stderr)
        return None


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.paf> <output_prefix> [options]",
              file=sys.stderr)
        print(f"  --reference <genome.fa>  Reference FASTA (auto-generates position map)",
              file=sys.stderr)
        print(f"  --positions <pos.tsv>    Pre-computed position map (from physlr index-positions)",
              file=sys.stderr)
        print(f"  --caption                Show input filename as caption",
              file=sys.stderr)
        sys.exit(1)

    input_paf = sys.argv[1]
    output_prefix = sys.argv[2]
    show_caption = '--caption' in sys.argv

    pos_map = None

    # --positions takes precedence over --reference
    if '--positions' in sys.argv:
        pos_idx = sys.argv.index('--positions')
        if pos_idx + 1 < len(sys.argv):
            pos_file = sys.argv[pos_idx + 1]
            print(f"Loading position map: {pos_file}")
            pos_map = load_position_map(pos_file)
            print(f"  Loaded {len(pos_map)} chromosomes")
    elif '--reference' in sys.argv:
        ref_idx = sys.argv.index('--reference')
        if ref_idx + 1 < len(sys.argv):
            ref_file = sys.argv[ref_idx + 1]
            pos_file = _generate_position_map(ref_file)
            if pos_file:
                print(f"Loading position map: {pos_file}")
                pos_map = load_position_map(pos_file)
                print(f"  Loaded {len(pos_map)} chromosomes")

    records = read_paf(input_paf)
    print(f"Read {len(records)} records")

    records = [r for r in records if r['Qname'] in ALLOWED_SET]
    print(f"After chromosome filter: {len(records)} records")

    records = [r for r in records if r['Tlength'] >= 50]

    records.sort(key=lambda r: (-r['Tlength'], r['Tname'], r['Tstart'], -r['Matches']))
    seen = set()
    deduped = []
    for r in records:
        key = (r['Tname'], r['Tstart'])
        if key not in seen:
            seen.add(key)
            deduped.append(r)
    records = deduped

    chained = chain_alignments(records)
    if not chained:
        print("No segments to plot.")
        sys.exit(1)

    plot_backbone(chained, output_prefix, show_caption, input_paf, pos_map=pos_map)
    plot_reference(chained, output_prefix, show_caption, input_paf, pos_map=pos_map)
    print("Done!")


if __name__ == '__main__':
    main()
