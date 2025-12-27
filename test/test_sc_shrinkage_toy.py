#!/usr/bin/env python3
import math
import os
import subprocess
import tempfile
from pathlib import Path


def write_mtx(path, nrows, ncols, entries):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("%%MatrixMarket matrix coordinate real general\n")
        handle.write(f"{nrows} {ncols} {len(entries)}\n")
        for row, col, value in entries:
            handle.write(f"{row} {col} {value}\n")


def read_mtx(path):
    entries = {}
    nrows = ncols = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if not parts:
                continue
            if nrows == 0:
                nrows, ncols = int(parts[0]), int(parts[1])
                continue
            row, col, value = int(parts[0]), int(parts[1]), float(parts[2])
            entries[(row, col)] = value
    return nrows, ncols, entries


def main():
    repo_root = Path(__file__).resolve().parents[1]
    dmtools = repo_root / "dmtools"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        mtx_path = tmp / "toy.mC.mtx"
        cov_path = tmp / "toy.cov.mtx"
        barcodes_path = tmp / "toy.barcodes.tsv"
        features_path = tmp / "toy.features.tsv"
        out_prefix = tmp / "toy_out"

        write_mtx(
            mtx_path,
            2,
            2,
            [
                (1, 1, 0),
                (2, 1, 1),
                (1, 2, 0),
                (2, 2, 0),
            ],
        )
        write_mtx(
            cov_path,
            2,
            2,
            [
                (1, 1, 1),
                (2, 1, 1),
                (1, 2, 10),
                (2, 2, 10),
            ],
        )

        barcodes_path.write_text("cell1\ncell2\n", encoding="utf-8")
        features_path.write_text(
            "region1\tchr1\t0\t10\tregion1\nregion2\tchr1\t10\t20\tregion2\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                str(dmtools),
                "sc-shrinkage",
                "--mtx",
                str(mtx_path),
                "--cov",
                str(cov_path),
                "--barcodes",
                str(barcodes_path),
                "--features",
                str(features_path),
                "--out",
                str(out_prefix),
                "--prior-strength",
                "20",
                "--min-cov",
                "1",
            ],
            check=True,
        )

        _, _, entries = read_mtx(tmp / "toy_out.mean.mtx")
        expected_cell1_region1 = 10.0 / 21.0
        expected_cell2_region1 = 11.0 / 21.0
        epsilon = 1e-6

        assert math.isclose(entries[(1, 1)], expected_cell1_region1, rel_tol=0, abs_tol=epsilon)
        assert math.isclose(entries[(2, 1)], expected_cell2_region1, rel_tol=0, abs_tol=epsilon)

        for key in [(1, 2), (2, 2)]:
            value = entries[key]
            assert 0.0 <= value <= 1.0
            assert value <= 1e-4


if __name__ == "__main__":
    main()
