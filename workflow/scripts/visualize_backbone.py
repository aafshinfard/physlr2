#!/usr/bin/env python3
"""
Visualize the Physlr physical map (backbone paths).

Generates an SVG visualization of the backbone paths showing:
- Each path as a horizontal chain of molecules
- Path lengths annotated
- Color-coded by path size

Usage: python3 visualize_backbone.py <backbone.path> -o <output.svg>
"""

import argparse
import sys


def read_paths(path):
    paths = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                paths.append(line.split())
    return paths


def generate_svg(paths, output, title="Physical Map"):
    """Generate an SVG visualization of backbone paths."""
    if not paths:
        print("No paths to visualize", file=sys.stderr)
        return

    # Sort by length descending
    paths.sort(key=len, reverse=True)

    # Layout parameters
    margin = 40
    node_r = 3
    node_spacing = 8
    path_spacing = 20
    label_width = 120

    max_path_len = max(len(p) for p in paths)
    width = margin * 2 + label_width + max_path_len * node_spacing + 100
    height = margin * 2 + len(paths) * path_spacing + 40

    # Color scale based on path length
    def path_color(length):
        if length >= 30:
            return "#2166ac"
        elif length >= 20:
            return "#4393c3"
        elif length >= 15:
            return "#92c5de"
        elif length >= 10:
            return "#d1e5f0"
        else:
            return "#f4a582"

    with open(output, "w") as f:
        f.write(f'<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" '
                f'width="{width}" height="{height}" '
                f'viewBox="0 0 {width} {height}">\n')
        f.write(f'<style>\n')
        f.write(f'  text {{ font-family: monospace; font-size: 11px; }}\n')
        f.write(f'  .title {{ font-size: 16px; font-weight: bold; }}\n')
        f.write(f'  .label {{ font-size: 10px; fill: #333; }}\n')
        f.write(f'</style>\n')
        f.write(f'<rect width="{width}" height="{height}" fill="white"/>\n')

        # Title
        f.write(f'<text x="{width//2}" y="{margin - 10}" '
                f'text-anchor="middle" class="title">{title}</text>\n')

        # Summary
        total_mols = sum(len(p) for p in paths)
        f.write(f'<text x="{width//2}" y="{margin + 5}" '
                f'text-anchor="middle" class="label">'
                f'{len(paths)} paths, {total_mols} molecules</text>\n')

        # Draw paths
        y_start = margin + 25
        for i, path in enumerate(paths):
            y = y_start + i * path_spacing
            color = path_color(len(path))

            # Label
            f.write(f'<text x="{margin}" y="{y + 4}" class="label">'
                    f'Path {i+1} ({len(path)})</text>\n')

            # Draw chain
            x_start = margin + label_width
            for j in range(len(path)):
                x = x_start + j * node_spacing
                f.write(f'<circle cx="{x}" cy="{y}" r="{node_r}" '
                        f'fill="{color}" stroke="#333" stroke-width="0.5"/>\n')
                if j > 0:
                    x_prev = x_start + (j - 1) * node_spacing
                    f.write(f'<line x1="{x_prev + node_r}" y1="{y}" '
                            f'x2="{x - node_r}" y2="{y}" '
                            f'stroke="{color}" stroke-width="1.5"/>\n')

        f.write('</svg>\n')

    print(f"SVG written to {output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Backbone path file")
    parser.add_argument("-o", "--output", required=True, help="Output SVG file")
    parser.add_argument("-t", "--title", default="Physical Map",
                        help="Plot title")
    args = parser.parse_args()

    paths = read_paths(args.input)
    generate_svg(paths, args.output, args.title)


if __name__ == "__main__":
    main()
