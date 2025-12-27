#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

def load_yaml(path: Path):
    try:
        import yaml  # type: ignore
    except ImportError:
        raise SystemExit("PyYAML is required: pip install pyyaml")
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)

def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)

def run_cmd(cmd, env=None):
    start = time.perf_counter()
    proc = subprocess.run(cmd, shell=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    elapsed = time.perf_counter() - start
    return proc.returncode, elapsed, proc.stdout, proc.stderr

def capture_hardware(out_dir: Path):
    info = {}
    cmds = {
        "uname": "uname -a",
        "lscpu": "lscpu",
        "mem": "free -h",
        "disk": "lsblk -o NAME,SIZE,TYPE,MOUNTPOINT -p",
        "date": "date -Is",
    }
    for key, cmd in cmds.items():
        rc, _, out, err = run_cmd(cmd)
        info[key] = out.strip() if rc == 0 else err.strip()
    with (out_dir / "hardware.txt").open("w", encoding="utf-8") as fh:
        for key, val in info.items():
            fh.write(f"## {key}\n{val}\n\n")

def drop_caches():
    cmd = "sudo sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches'"
    rc, _, _, err = run_cmd(cmd)
    return rc == 0, err

def file_size_bytes(path: Path):
    try:
        return path.stat().st_size
    except FileNotFoundError:
        return None

def human_size(num):
    if num is None:
        return "N/A"
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if num < 1024:
            return f"{num:.2f}{unit}"
        num /= 1024
    return f"{num:.2f}PB"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config_path = Path(args.config)
    cfg = load_yaml(config_path)

    run_id = cfg.get("run_id") or dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = Path("bench/results") / run_id
    ensure_dir(out_dir)
    ensure_dir(out_dir / "outputs")

    shutil.copy2(config_path, out_dir / "config_used.yaml")
    capture_hardware(out_dir)

    cache_cfg = cfg.get("cache_strategy", {})
    cache_mode = cache_cfg.get("mode", "none")
    drop_cache_enabled = bool(cache_cfg.get("drop_caches", False))

    paths = cfg.get("paths", {})

    # Storage comparison
    storage_cfg = cfg.get("benchmarks", {}).get("storage", {})
    if storage_cfg.get("enabled", False):
        with (out_dir / "storage.csv").open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["dataset", "path", "bytes", "size_human"])
            for key in storage_cfg.get("datasets", []):
                p = Path(paths.get(key, ""))
                size = file_size_bytes(p)
                writer.writerow([key, str(p), size if size is not None else "", human_size(size)])

    perf_path = out_dir / "performance.csv"
    with perf_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["benchmark", "dataset", "iteration", "command", "seconds", "status", "stdout", "stderr"])

        for bench_name, bench_cfg in cfg.get("benchmarks", {}).items():
            if bench_name == "storage" or not bench_cfg.get("enabled", False):
                continue
            datasets = bench_cfg.get("datasets", [])
            repeats = int(bench_cfg.get("repeats", 1))
            command_template = bench_cfg.get("command")
            regions = bench_cfg.get("regions", [])

            for dataset_key in datasets:
                dm_path = paths.get(dataset_key)
                if not dm_path:
                    continue

                def run_with_command(cmd, iteration):
                    if cache_mode == "cold" and drop_cache_enabled:
                        ok, err = drop_caches()
                        if not ok:
                            writer.writerow([bench_name, dataset_key, iteration, cmd, "", "drop_cache_failed", "", err])
                    rc, elapsed, out, err = run_cmd(cmd)
                    status = "ok" if rc == 0 else f"error({rc})"
                    writer.writerow([bench_name, dataset_key, iteration, cmd, f"{elapsed:.6f}", status, out.strip(), err.strip()])

                if bench_name == "region_query" and regions:
                    for region in regions:
                        for i in range(repeats):
                            cmd = command_template.format(
                                dm=dm_path,
                                region=region,
                                regions_bed=paths.get("regions_bed", ""),
                                matrix_regions_bed=paths.get("matrix_regions_bed", ""),
                                output=paths.get("output_dir", out_dir / "outputs"),
                                output_prefix=str(Path(paths.get("output_dir", out_dir / "outputs")) / f"{dataset_key}_{bench_name}"),
                            )
                            run_with_command(cmd, i + 1)
                else:
                    for i in range(repeats):
                        cmd = command_template.format(
                            dm=dm_path,
                            region=regions[0] if regions else "",
                            regions_bed=paths.get("regions_bed", ""),
                            matrix_regions_bed=paths.get("matrix_regions_bed", ""),
                            output=str(Path(paths.get("output_dir", out_dir / "outputs")) / f"{dataset_key}_{bench_name}.tsv"),
                            output_prefix=str(Path(paths.get("output_dir", out_dir / "outputs")) / f"{dataset_key}_{bench_name}"),
                        )
                        run_with_command(cmd, i + 1)

    summary = {
        "run_id": run_id,
        "results_dir": str(out_dir),
        "cache_strategy": cache_cfg,
    }
    with (out_dir / "summary.json").open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(f"Results written to {out_dir}")

if __name__ == "__main__":
    main()
