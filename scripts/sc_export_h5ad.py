#!/usr/bin/env python3
import argparse
import os
import sys


def fail(msg: str, code: int = 1) -> None:
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(code)


def load_deps():
    try:
        import anndata as ad  # type: ignore
        import pandas as pd  # type: ignore
        import scipy.io  # type: ignore
    except ImportError as exc:
        fail(
            "Error: missing Python dependencies for h5ad export. "
            "Install with: python3 -m pip install anndata scipy pandas\n"
            f"Details: {exc}",
            code=2,
        )
    return ad, pd, scipy.io


def resolve_matrix_path(prefix: str) -> str:
    candidates = [f"{prefix}.matrix.mtx", f"{prefix}.mtx"]
    for path in candidates:
        if os.path.exists(path):
            return path
    fail(
        f"Error: missing matrix file. Expected {candidates[0]} (or fallback {candidates[1]}).",
        code=1,
    )
    return ""


def read_lines(path: str) -> list[str]:
    if not os.path.exists(path):
        fail(f"Error: missing required file: {path}")
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert sc-matrix bundle to h5ad.")
    parser.add_argument("--prefix", required=True, help="Output prefix for sc-matrix bundle.")
    parser.add_argument("--output", required=True, help="Output .h5ad path.")
    args = parser.parse_args()

    ad, pd, scipy_io = load_deps()

    matrix_path = resolve_matrix_path(args.prefix)
    barcodes_path = f"{args.prefix}.barcodes.tsv"
    features_path = f"{args.prefix}.features.tsv"
    obs_qc_path = f"{args.prefix}.obs_qc.tsv"

    barcodes = read_lines(barcodes_path)
    if not barcodes:
        fail(f"Error: no barcodes found in {barcodes_path}")

    if not os.path.exists(features_path):
        fail(f"Error: missing required file: {features_path}")
    features = pd.read_csv(
        features_path,
        sep="\t",
        header=None,
        names=["feature_id", "chrom", "start", "end", "name"],
        dtype={"feature_id": str, "chrom": str, "name": str},
    )
    if features.empty:
        fail(f"Error: no features found in {features_path}")

    if not os.path.exists(obs_qc_path):
        fail(f"Error: missing required file: {obs_qc_path}")
    obs_qc = pd.read_csv(obs_qc_path, sep="\t")
    if "cell_id" not in obs_qc.columns:
        fail(f"Error: obs_qc.tsv missing required column 'cell_id': {obs_qc_path}")

    matrix = scipy_io.mmread(matrix_path).tocsr()

    if matrix.shape != (len(barcodes), len(features)):
        fail(
            "Error: matrix dimensions do not match barcodes/features. "
            f"Matrix shape {matrix.shape}, barcodes {len(barcodes)}, features {len(features)}."
        )

    obs = pd.DataFrame(index=barcodes)
    obs_qc = obs_qc.set_index("cell_id")
    obs = obs.join(obs_qc, how="left")

    var = features.copy()
    var = var.set_index("feature_id")

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.write_h5ad(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
