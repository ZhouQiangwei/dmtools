#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

def load_csv(path):
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

def summarize_performance(rows):
    summary = {}
    for row in rows:
        key = (row["benchmark"], row["dataset"])
        if row["status"] != "ok":
            continue
        summary.setdefault(key, []).append(float(row["seconds"]))
    out = []
    for (bench, dataset), vals in summary.items():
        avg = sum(vals) / len(vals)
        out.append({"benchmark": bench, "dataset": dataset, "avg_seconds": f"{avg:.6f}", "n": str(len(vals))})
    return out

def write_table(rows, out_path):
    if not rows:
        return
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    figures_dir = results_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    perf_path = results_dir / "performance.csv"
    storage_path = results_dir / "storage.csv"

    if perf_path.exists():
        perf_rows = load_csv(perf_path)
        perf_summary = summarize_performance(perf_rows)
        write_table(perf_summary, figures_dir / "performance_summary.csv")

    if storage_path.exists():
        storage_rows = load_csv(storage_path)
        write_table(storage_rows, figures_dir / "storage_table.csv")

    meta = {"results_dir": str(results_dir)}
    with (figures_dir / "metadata.json").open("w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)

    print(f"Figures/tables written to {figures_dir}")

if __name__ == "__main__":
    main()
