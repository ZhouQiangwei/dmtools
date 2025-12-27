#!/usr/bin/env python3
import math
import statistics
import subprocess
import tempfile
from pathlib import Path


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
        genome_len = tmp / "toy.len"
        bed_path = tmp / "toy.bed"
        bedsimple_path = tmp / "toy.bedsimple"
        dm_path = tmp / "toy.dm"
        idmap_path = tmp / "toy.dm.idmap.tsv"
        matrix_prefix = tmp / "toy_matrix"
        shrink_prefix = tmp / "toy_shrink"

        genome_len.write_text("chr1\t20\n", encoding="utf-8")
        bed_path.write_text(
            "chr1\t0\t10\tregion1\n"
            "chr1\t10\t20\tregion2\n",
            encoding="utf-8",
        )
        bedsimple_path.write_text(
            "chr1\t1\t2\tcell1\t+\tCG\t0\t1\n"
            "chr1\t1\t2\tcell2\t+\tCG\t1\t1\n"
            "chr1\t11\t12\tcell1\t+\tCG\t0\t10\n"
            "chr1\t11\t12\tcell2\t+\tCG\t0\t10\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                str(dmtools),
                "mr2dm",
                "-C",
                "-S",
                "--Cx",
                "--Id",
                "-f",
                "bedsimple",
                "-g",
                str(genome_len),
                "-m",
                str(bedsimple_path),
                "-o",
                str(dm_path),
                "--idmap-out",
                str(idmap_path),
            ],
            check=True,
        )

        subprocess.run(
            [
                str(dmtools),
                "sc-matrix",
                "-i",
                str(dm_path),
                "-o",
                str(matrix_prefix),
                "--bed",
                str(bed_path),
                "--emit-counts",
            ],
            check=True,
        )

        subprocess.run(
            [
                str(dmtools),
                "sc-shrinkage",
                "--mtx",
                f"{matrix_prefix}.mC.mtx",
                "--cov",
                f"{matrix_prefix}.cov.mtx",
                "--barcodes",
                f"{matrix_prefix}.barcodes.tsv",
                "--features",
                f"{matrix_prefix}.features.tsv",
                "--out",
                str(shrink_prefix),
                "--prior-strength",
                "20",
                "--min-cov",
                "1",
            ],
            check=True,
        )

        _, _, mc_entries = read_mtx(f"{matrix_prefix}.mC.mtx")
        _, _, cov_entries = read_mtx(f"{matrix_prefix}.cov.mtx")
        _, _, mean_entries = read_mtx(f"{shrink_prefix}.mean.mtx")

        low_cov_raw = []
        low_cov_mean = []
        for key, cov in cov_entries.items():
            if cov == 0:
                continue
            raw = mc_entries[key] / cov
            mean = mean_entries[key]
            assert 0.0 <= mean <= 1.0
            if math.isclose(cov, 1.0, rel_tol=0, abs_tol=1e-6):
                low_cov_raw.append(raw)
                low_cov_mean.append(mean)

        assert len(low_cov_raw) >= 2
        for raw, mean in zip(low_cov_raw, low_cov_mean):
            assert abs(mean - 0.5) < abs(raw - 0.5)

        raw_var = statistics.pvariance(low_cov_raw)
        mean_var = statistics.pvariance(low_cov_mean)
        assert mean_var < raw_var


if __name__ == "__main__":
    main()
