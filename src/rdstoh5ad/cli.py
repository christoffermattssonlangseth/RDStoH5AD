from __future__ import annotations

import argparse
import importlib.resources
import os
import subprocess
import sys
from pathlib import Path

from . import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rds2h5ad",
        description="Inspect and convert .rds inputs to .h5ad using an R backend.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--rscript",
        default=os.environ.get("RDS2H5AD_RSCRIPT", "Rscript"),
        help="Path to Rscript. Defaults to RDS2H5AD_RSCRIPT or Rscript on PATH.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect_parser = subparsers.add_parser(
        "inspect",
        help="Inspect the input and print machine-readable JSON.",
    )
    inspect_parser.add_argument("input", help="Path to the input .rds or .RData file.")
    inspect_parser.add_argument(
        "--assay",
        help="Requested assay to inspect. If omitted, the backend chooses a default.",
    )
    inspect_parser.set_defaults(func=run_inspect)

    convert_parser = subparsers.add_parser(
        "convert",
        help="Convert the input object to .h5ad.",
    )
    convert_parser.add_argument("input", help="Path to the input .rds or .RData file.")
    convert_parser.add_argument("output", help="Path to the output .h5ad file.")
    convert_parser.add_argument(
        "--assay",
        help="Requested assay. If omitted, the backend chooses a default assay.",
    )
    convert_parser.add_argument(
        "--x-layer",
        help="Layer or assay to map to AnnData X. Defaults to data/logcounts when present, otherwise counts.",
    )
    convert_parser.add_argument(
        "--layers",
        help="Comma-separated layers/assays to retain in the H5AD file. Defaults to all compatible layers for the selected assay.",
    )
    convert_parser.add_argument(
        "--reduced-dims",
        help="Comma-separated reducedDims/reductions to export. Defaults to all available embeddings.",
    )
    convert_parser.add_argument(
        "--no-spatial",
        action="store_true",
        help="Skip adding inferred spatial coordinates to obsm.",
    )
    convert_parser.set_defaults(func=run_convert)

    validate_parser = subparsers.add_parser(
        "validate",
        help="Read a .h5ad file and print a compact structural summary.",
    )
    validate_parser.add_argument("input", help="Path to the input .h5ad file.")
    validate_parser.set_defaults(func=run_validate)

    return parser


def backend_script_path() -> str:
    resource = importlib.resources.files("rdstoh5ad").joinpath("r").joinpath("rds_to_h5ad.R")
    with importlib.resources.as_file(resource) as path:
        return str(path)


def run_backend(rscript: str, backend_args: list[str]) -> subprocess.CompletedProcess[str]:
    command = [rscript, backend_script_path(), *backend_args]
    return subprocess.run(command, capture_output=True, text=True)


def die(message: str, stderr: str | None = None, exit_code: int = 1) -> int:
    print(message, file=sys.stderr)
    if stderr:
        stderr = stderr.strip()
        if stderr:
            print(stderr, file=sys.stderr)
    return exit_code


def run_inspect(args: argparse.Namespace) -> int:
    backend_args = ["inspect", "--input", str(Path(args.input))]
    if args.assay:
        backend_args.extend(["--assay", args.assay])

    result = run_backend(args.rscript, backend_args)
    if result.returncode != 0:
        return die("Inspection failed.", result.stderr, result.returncode)

    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)
    sys.stdout.write(result.stdout)
    return 0


def run_convert(args: argparse.Namespace) -> int:
    backend_args = [
        "convert",
        "--input",
        str(Path(args.input)),
        "--output",
        str(Path(args.output)),
    ]
    if args.assay:
        backend_args.extend(["--assay", args.assay])
    if args.x_layer:
        backend_args.extend(["--x-layer", args.x_layer])
    if args.layers:
        backend_args.extend(["--layers", args.layers])
    if args.reduced_dims:
        backend_args.extend(["--reduced-dims", args.reduced_dims])
    if args.no_spatial:
        backend_args.append("--no-spatial")

    result = run_backend(args.rscript, backend_args)
    if result.returncode != 0:
        return die("Conversion failed.", result.stderr, result.returncode)

    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)
    return 0


def run_validate(args: argparse.Namespace) -> int:
    backend_args = ["validate", "--input", str(Path(args.input))]
    result = run_backend(args.rscript, backend_args)
    if result.returncode != 0:
        return die("Validation failed.", result.stderr, result.returncode)

    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)
    sys.stdout.write(result.stdout)
    return 0


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
