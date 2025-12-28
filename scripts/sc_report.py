#!/usr/bin/env python3
"""Generate simple single-cell QC report (SVG).

Inputs:
  --qc        obs_qc.tsv from sc-matrix/sc-qc
  --features  features.tsv (optional; counts features)
  --out       output SVG path

Outputs a two-panel SVG: coverage histogram and CG/CHH ratio distribution.
"""
import argparse
import math
import os
from typing import List, Tuple


def read_tsv(path: str) -> Tuple[List[str], List[List[str]]]:
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"Empty TSV: {path}")
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    return header, rows


def parse_qc(path: str):
    header, rows = read_tsv(path)
    col_idx = {name: i for i, name in enumerate(header)}
    def get(col: str, row):
        return float(row[col_idx[col]]) if col in col_idx else 0.0
    coverages = [get("total_coverage", r) for r in rows]
    cg_counts = [get("n_sites_cg", r) for r in rows]
    chh_counts = [get("n_sites_chh", r) for r in rows]
    ratios = []
    for cg, chh in zip(cg_counts, chh_counts):
        if cg == 0 and chh == 0:
            ratios.append(0.0)
        elif chh == 0:
            ratios.append(cg)
        else:
            ratios.append(cg / chh)
    n_sites = [get("n_sites", r) for r in rows]
    return coverages, ratios, n_sites


def count_features(path: str) -> int:
    if not path:
        return 0
    if not os.path.exists(path):
        return 0
    with open(path, "r", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def histogram(values: List[float], bins: int = 20) -> Tuple[List[float], List[float]]:
    if not values:
        return [], []
    vmax = max(values)
    if vmax <= 0:
        vmax = 1.0
    bin_edges = [vmax * i / bins for i in range(bins + 1)]
    counts = [0 for _ in range(bins)]
    for v in values:
        idx = min(int((v / vmax) * bins), bins - 1)
        counts[idx] += 1
    return bin_edges, counts


def render_svg(out_path: str, coverages: List[float], ratios: List[float], n_sites: List[float], feature_count: int) -> None:
    width, height = 900, 450
    padding = 50
    svg = []
    svg.append(f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>")
    svg.append("<style> text { font-family: Helvetica, Arial, sans-serif; font-size: 12px; } </style>")

    # Coverage histogram (left)
    hist_x, hist_y = histogram(coverages, bins=20)
    panel_w = (width - 3 * padding) // 2
    panel_h = height - 2 * padding
    x0, y0 = padding, padding
    svg.append(f"<text x='{x0}' y='{y0 - 15}' font-size='14' font-weight='bold'>Coverage distribution</text>")
    if hist_y:
        max_count = max(hist_y)
        for i, count in enumerate(hist_y):
            bar_h = 0 if max_count == 0 else (count / max_count) * (panel_h - 30)
            x = x0 + i * (panel_w / len(hist_y))
            y = y0 + (panel_h - bar_h)
            w = max(2, panel_w / len(hist_y) - 2)
            svg.append(f"<rect x='{x:.2f}' y='{y:.2f}' width='{w:.2f}' height='{bar_h:.2f}' fill='#4a90e2' />")
        svg.append(f"<text x='{x0}' y='{y0 + panel_h + 20}'>0</text>")
        svg.append(f"<text x='{x0 + panel_w - 20}' y='{y0 + panel_h + 20}'>max={max(hist_x):.0f}</text>")
    else:
        svg.append(f"<text x='{x0 + 20}' y='{y0 + 40}'>No coverage data</text>")

    # Ratio box (right)
    x0 = 2 * padding + panel_w
    svg.append(f"<text x='{x0}' y='{padding - 15}' font-size='14' font-weight='bold'>CG/CHH ratio & features</text>")
    rect_w = panel_w
    rect_h = panel_h
    svg.append(f"<rect x='{x0}' y='{padding}' width='{rect_w}' height='{rect_h}' fill='none' stroke='#ccc' />")
    if ratios:
        avg_ratio = sum(ratios) / len(ratios)
        avg_sites = sum(n_sites) / len(n_sites) if n_sites else 0
        cx = x0 + rect_w * 0.1
        cy = padding + rect_h * 0.3
        svg.append(f"<text x='{cx}' y='{cy}' font-size='24' font-weight='bold'>{avg_ratio:.2f}</text>")
        svg.append(f"<text x='{cx}' y='{cy + 20}' fill='#666'>mean CG/CHH</text>")
        svg.append(f"<text x='{cx}' y='{cy + 60}' font-size='16'>median coverage ~ {sorted(coverages)[len(coverages)//2]:.0f}</text>")
        svg.append(f"<text x='{cx}' y='{cy + 80}' font-size='16'>avg detected sites ~ {avg_sites:.0f}</text>")
    else:
        svg.append(f"<text x='{x0 + 20}' y='{padding + 40}'>No CG/CHH data</text>")
    if feature_count:
        svg.append(f"<text x='{x0 + 10}' y='{padding + rect_h - 15}' font-size='16'>features: {feature_count}</text>")

    svg.append("</svg>")
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\n".join(svg))


def main() -> int:
    parser = argparse.ArgumentParser(description="Single-cell QC report (SVG)")
    parser.add_argument("--qc", required=True, help="obs_qc.tsv path")
    parser.add_argument("--features", default=None, help="features.tsv path (optional)")
    parser.add_argument("--out", required=True, help="output SVG path")
    args = parser.parse_args()

    coverages, ratios, n_sites = parse_qc(args.qc)
    feature_count = count_features(args.features) if args.features else 0
    render_svg(args.out, coverages, ratios, n_sites, feature_count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
