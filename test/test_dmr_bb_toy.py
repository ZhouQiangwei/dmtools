#!/usr/bin/env python3
import subprocess
import tempfile
from pathlib import Path


def main():
    repo_root = Path(__file__).resolve().parents[1]
    dmtools = repo_root / "dmtools"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        input_path = tmp / "dmr_counts.tsv"
        output_path = tmp / "dmr_out.tsv"

        input_path.write_text(
            "chr1\t0\t10\t80\t100\t20\t100\n"
            "chr1\t10\t20\t50\t100\t50\t100\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                str(dmtools),
                "dmr-bb",
                "--input",
                str(input_path),
                "--out",
                str(output_path),
                "--rho",
                "0.0",
            ],
            check=True,
        )

        lines = output_path.read_text(encoding="utf-8").strip().splitlines()
        assert len(lines) == 3
        header = lines[0].split("\t")
        assert header[:6] == ["chrom", "start", "end", "delta", "p_value", "q_value"]

        first = lines[1].split("\t")
        assert first[0] == "chr1"
        pval = float(first[4])
        assert pval < 1e-6


if __name__ == "__main__":
    main()
