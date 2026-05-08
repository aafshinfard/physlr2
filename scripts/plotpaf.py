#!/usr/bin/env python3
"""Plot backbone-vs-reference PAF alignment.

Generates two plots:
  1. Backbone coverage: backbone paths colored by reference chromosome
  2. Reference coverage: chromosomes colored by backbone path

Faithful to the original Physlr plotpaf.rmd visualization.

Usage: plotpaf.py <input.paf> <output_prefix>
  Produces: <output_prefix>.backbone.png and <output_prefix>.reference.png
"""
import sys
import csv
from collections import defaultdict

def read_paf(path):
    """Read PAF file. Returns list of dicts."""
    records = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                continue
            records.append({
                'qname': parts[0],
                'qlength': int(parts[1]),
                'qstart': int(parts[2]),
                'qend': int(parts[3]),
                'orientation': parts[4],
                'tname': parts[5],
                'tlength': int(parts[6]),
                'tstart': int(parts[7]),
                'tend': int(parts[8]),
                'score': int(parts[9]),
                'length': int(parts[10]),
                'mapq': int(parts[11]),
            })
    return records

def chain_alignments(records, min_barcodes=50):
    """Chain alignments using the same algorithm as the original plotpaf.rmd.

    Steps (matching the R implementation):
    1. Filter by mapq >= 1
    2. Sort by (desc(Tlength), Tname, Tstart)
    3. Group by (Qname, Tname), keep only records that overlap with
       their neighbor on the Q axis
    4. Detect segment boundaries where Qdistance >= 500000 or
       Tdistance >= 10
    5. Summarise each segment: merge Q/T ranges, count records
    6. Filter segments with count >= min_barcodes (default 50)
    """
    mapq_threshold = 1
    qdistance_threshold = 500000
    tdistance_threshold = 10

    def is_overlap(a_start, a_end, b_start, b_end):
        return (a_start <= b_start <= a_end or a_start <= b_end <= a_end or
                b_start <= a_start <= b_end or b_start <= a_end <= b_end)

    def interval_distance(a_start, a_end, b_start, b_end):
        if is_overlap(a_start, a_end, b_start, b_end):
            return 0
        return max(b_start - a_end, a_start - b_end)

    # Step 1: filter by mapq
    filtered = [r for r in records if r['mapq'] >= mapq_threshold]

    # Step 2: sort by desc(tlength), tname, tstart
    filtered.sort(key=lambda r: (-r['tlength'], r['tname'], r['tstart']))

    # Step 3: group by (qname, tname)
    groups = defaultdict(list)
    for r in filtered:
        groups[(r['qname'], r['tname'])].append(r)

    # For each group, keep only records that overlap with neighbor on Q axis,
    # then detect segment boundaries and summarise
    chained = []
    for (qname, tname), recs in groups.items():
        # Already sorted by tstart within group (from global sort)
        # Filter: keep records that overlap with lag or lead on Q axis
        keep = []
        for i, r in enumerate(recs):
            has_overlap = False
            if i > 0:
                prev = recs[i - 1]
                if is_overlap(prev['qstart'], prev['qend'], r['qstart'], r['qend']):
                    has_overlap = True
            if i < len(recs) - 1:
                nxt = recs[i + 1]
                if is_overlap(r['qstart'], r['qend'], nxt['qstart'], nxt['qend']):
                    has_overlap = True
            if has_overlap:
                keep.append(r)

        if not keep:
            continue

        # Detect segment boundaries
        segments = []
        current_seg = [keep[0]]
        for i in range(1, len(keep)):
            prev = keep[i - 1]
            curr = keep[i]
            qdist = interval_distance(prev['qstart'], prev['qend'],
                                      curr['qstart'], curr['qend'])
            tdist = interval_distance(prev['tstart'], prev['tend'],
                                      curr['tstart'], curr['tend'])
            if qdist >= qdistance_threshold or tdist >= tdistance_threshold:
                segments.append(current_seg)
                current_seg = [curr]
            else:
                current_seg.append(curr)
        segments.append(current_seg)

        # Summarise each segment
        for seg in segments:
            if len(seg) < min_barcodes:
                continue
            chained.append({
                'qname': qname,
                'tname': tname,
                'qlength': seg[0]['qlength'],
                'tlength': seg[0]['tlength'],
                'qstart': min(r['qstart'] for r in seg),
                'qend': max(r['qend'] for r in seg),
                'tstart': min(r['tstart'] for r in seg),
                'tend': max(r['tend'] for r in seg),
                'score': sum(r['score'] for r in seg),
                'orientation': max(set(r['orientation'] for r in seg),
                                   key=lambda o: sum(1 for r in seg if r['orientation'] == o)),
                'barcodes': len(seg),
            })

    return chained

def chr_sort_key(name):
    """Sort chromosomes: numeric first, then alpha."""
    name = name.replace('chr', '')
    try:
        return (0, int(name), '')
    except ValueError:
        return (1, 0, name)

def filter_main_chromosomes(chained):
    """Keep only main chromosomes (chr1-22, chrX, chrY, or numeric 1-22, X, Y)."""
    main = set()
    for i in range(1, 23):
        main.add(f'chr{i}')
        main.add(str(i))
    main.update(['chrX', 'chrY', 'X', 'Y'])
    return [r for r in chained if r['qname'] in main]

def get_top_backbones(chained, n=30):
    """Return the top N backbone paths by total score."""
    scores = defaultdict(int)
    for r in chained:
        scores[r['tname']] += r['score']
    top = sorted(scores, key=lambda t: -scores[t])[:n]
    top_set = set(top)
    return [r for r in chained if r['tname'] in top_set], top

def plot_with_matplotlib(records, chained, output_prefix):
    """Generate plots using matplotlib."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import matplotlib.colors as mcolors

    # Color palette (26 distinct colors)
    palette = list(mcolors.TABLEAU_COLORS.values()) + [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
        '#469990', '#dcbeff', '#9A6324', '#800000', '#aaffc3',
        '#808000',
    ]

    # Filter to main chromosomes and top backbones
    main_chained = filter_main_chromosomes(chained)
    if not main_chained:
        main_chained = chained  # fallback if no main chroms found

    top_chained, top_tnames = get_top_backbones(main_chained, n=30)
    if not top_chained:
        print("No data to plot after filtering")
        return

    qnames = sorted(set(r['qname'] for r in main_chained), key=chr_sort_key)
    qname_color = {q: palette[i % len(palette)] for i, q in enumerate(qnames)}

    # --- Plot 1: Backbone coverage (colored by reference chromosome) ---
    tnames = top_tnames
    tname_idx = {t: i for i, t in enumerate(tnames)}

    fig, ax = plt.subplots(figsize=(14, max(6, len(tnames) * 0.35 + 2)))
    for r in top_chained:
        if r['tname'] in tname_idx:
            tidx = tname_idx[r['tname']]
            ax.add_patch(Rectangle(
                (r['tstart'], tidx), max(r['tend'] - r['tstart'], 1), 0.8,
                facecolor=qname_color.get(r['qname'], '#999999'),
                edgecolor='none', alpha=0.8
            ))

    max_tlen = max((r['tlength'] for r in top_chained), default=100)
    ax.set_xlim(0, max_tlen)
    ax.set_ylim(-0.5, len(tnames) + 0.5)
    ax.set_xlabel('Position along backbone (nodes)')
    ax.set_ylabel('Backbone path')
    ax.set_yticks([i + 0.4 for i in range(len(tnames))])
    ax.set_yticklabels(tnames, fontsize=8)
    ax.invert_yaxis()
    ax.set_title('Backbone coverage by reference chromosome (top 30 paths)')

    handles = [Rectangle((0, 0), 1, 1, facecolor=qname_color[q]) for q in qnames]
    ax.legend(handles, qnames, loc='upper right', fontsize=6,
              ncol=max(1, len(qnames) // 8), title='Chromosome')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.backbone.png', dpi=150)
    plt.close()
    print(f'Saved {output_prefix}.backbone.png')

    # --- Plot 2: Reference coverage (colored by backbone path) ---
    tname_color = {t: palette[i % len(palette)] for i, t in enumerate(tnames)}

    fig, ax = plt.subplots(figsize=(14, max(6, len(qnames) * 0.5 + 2)))
    for r in top_chained:
        if r['qname'] in qnames:
            qidx = qnames.index(r['qname'])
            ax.add_patch(Rectangle(
                (r['qstart'], qidx), max(r['qend'] - r['qstart'], 1), 0.8,
                facecolor=tname_color.get(r['tname'], '#999999'),
                edgecolor='none', alpha=0.7
            ))

    qlengths = {}
    for r in main_chained:
        qlengths[r['qname']] = max(qlengths.get(r['qname'], 0), r['qlength'])

    ax.set_xlim(0, max(qlengths.values(), default=100) * 1.05)
    ax.set_ylim(-0.5, len(qnames) + 0.5)
    ax.set_xlabel('Position along reference (minimizer index)')
    ax.set_ylabel('Chromosome')
    ax.set_yticks([i + 0.4 for i in range(len(qnames))])
    ax.set_yticklabels(qnames, fontsize=8)
    ax.invert_yaxis()
    ax.set_title('Reference coverage by backbone path (top 30 paths)')

    handles = [Rectangle((0, 0), 1, 1, facecolor=tname_color[t]) for t in tnames[:26]]
    ax.legend(handles, [f'Path {t}' for t in tnames[:26]], loc='upper right',
              fontsize=6, title='Backbone', ncol=max(1, len(tnames) // 8))

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.reference.png', dpi=150)
    plt.close()
    print(f'Saved {output_prefix}.reference.png')

def plot_ascii(records, chained, output_prefix):
    """Fallback ASCII summary when matplotlib is unavailable."""
    print("matplotlib not available, generating text summary only")
    qnames = sorted(set(r['qname'] for r in chained), key=chr_sort_key)
    tnames = sorted(set(r['tname'] for r in chained), key=lambda x: int(x) if x.isdigit() else 0)

    with open(f'{output_prefix}.summary.txt', 'w') as f:
        f.write("=== Backbone coverage by chromosome ===\n")
        for tname in tnames:
            recs = [r for r in chained if r['tname'] == tname]
            chroms = defaultdict(int)
            for r in recs:
                chroms[r['qname']] += r['score']
            top = sorted(chroms.items(), key=lambda x: -x[1])[:5]
            f.write(f"  Backbone {tname}: {', '.join(f'{c}({s})' for c, s in top)}\n")

        f.write("\n=== Reference coverage by backbone ===\n")
        for qname in qnames:
            recs = [r for r in chained if r['qname'] == qname]
            paths = defaultdict(int)
            for r in recs:
                paths[r['tname']] += r['score']
            top = sorted(paths.items(), key=lambda x: -x[1])[:5]
            f.write(f"  {qname}: {', '.join(f'path{t}({s})' for t, s in top)}\n")

    print(f'Saved {output_prefix}.summary.txt')

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.paf> <output_prefix>", file=sys.stderr)
        sys.exit(1)

    paf_file = sys.argv[1]
    output_prefix = sys.argv[2]
    min_barcodes = int(sys.argv[3]) if len(sys.argv) > 3 else 50

    records = read_paf(paf_file)
    print(f"Read {len(records)} PAF records")

    chained = chain_alignments(records, min_barcodes=min_barcodes)
    print(f"Chained into {len(chained)} segments (min_barcodes={min_barcodes})")

    if not chained:
        print("No chained segments to plot. Try lowering min_barcodes.")
        sys.exit(0)

    try:
        plot_with_matplotlib(records, chained, output_prefix)
    except ImportError:
        plot_ascii(records, chained, output_prefix)

if __name__ == '__main__':
    main()
