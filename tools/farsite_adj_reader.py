#!/usr/bin/env python3
"""
farsite_adj_reader.py – Read and apply FARSITE fuel adjustment (.adj) files.

A FARSITE fuel adjustment file specifies a per-fuel-model rate-of-spread (ROS)
multiplier.  This script can:

  1. **Inspect** an existing ``.adj`` file (print the table to stdout).
  2. **Generate** a template ``.adj`` file for a given set of fuel model numbers.
  3. **Apply** the adjustments to a Rothermel parameter CSV and write the
     modified version to a new file (useful for calibrating landscape runs
     without recompiling the solver).

The ``.adj`` file format (whitespace-delimited; ``#`` or ``!`` lines are
comments) ::

    [optional: <num_fuel_models>]
    <fuel_model_number>  <adj_factor>
    ...

The ``adj_factor`` is a multiplicative scale applied to the no-wind, no-slope
Rothermel rate of spread ``R₀``:

    R_adjusted = R_original × adj_factor

**References**
    Finney, M.A. (2004). *FARSITE: Fire Area Simulator—Model Development and
    Evaluation*. USDA Forest Service Research Paper RMRS-RP-4.

Usage examples
--------------
Inspect an existing .adj file::

    python3 farsite_adj_reader.py --adj fire.adj

Generate a template for FBFM13 models 1–13::

    python3 farsite_adj_reader.py --generate --fuel-system 13 \\
        --output template_13.adj

Generate a template for specific models::

    python3 farsite_adj_reader.py --generate --models 4 9 10 \\
        --output my_adj.adj

Apply adjustments and write a modified CSV (for inspection / calibration)::

    python3 farsite_adj_reader.py --adj fire.adj \\
        --apply-csv rothermel_params.csv \\
        --output rothermel_adjusted.csv
"""

import argparse
import csv
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Core parsing
# ---------------------------------------------------------------------------

def parse_adj_file(filepath: str) -> dict[int, float]:
    """
    Parse a FARSITE .adj file.

    Returns a dict mapping *fuel_model_code* (int) to *adj_factor* (float).
    Lines starting with ``#`` or ``!`` are treated as comments.
    A single-token line that appears first is treated as the optional
    *num_models* declaration and is skipped.

    Parameters
    ----------
    filepath : str
        Path to the .adj file.

    Returns
    -------
    dict[int, float]
        ``{fuel_model_code: adj_factor}``
    """
    path = Path(filepath)
    if not path.exists():
        sys.exit(f"ERROR: file not found: {filepath}")

    adjustments: dict[int, float] = {}
    count_line_consumed = False

    with open(path) as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line or line[0] in ("#", "!"):
                continue
            tokens = line.split()
            if len(tokens) == 1:
                if not count_line_consumed:
                    count_line_consumed = True  # optional num_models
                continue
            if len(tokens) >= 2:
                try:
                    model = int(tokens[0])
                    factor = float(tokens[1])
                except ValueError:
                    print(f"WARNING: cannot parse line {lineno}: {raw.rstrip()}",
                          file=sys.stderr)
                    continue
                if factor < 0:
                    print(f"WARNING: negative adj_factor {factor} for model "
                          f"{model} on line {lineno} – clamped to 0",
                          file=sys.stderr)
                    factor = 0.0
                adjustments[model] = factor

    return adjustments


def print_adj_table(adjustments: dict[int, float]) -> None:
    """Pretty-print an adjustment table."""
    print(f"{'Fuel model':>12}  {'adj_factor':>12}")
    print("-" * 26)
    for model in sorted(adjustments):
        marker = "  ← modified" if adjustments[model] != 1.0 else ""
        print(f"{model:>12}  {adjustments[model]:>12.4f}{marker}")
    print(f"\nTotal entries: {len(adjustments)}")


# ---------------------------------------------------------------------------
# Template generation
# ---------------------------------------------------------------------------

# Anderson (1982) FBFM13 model numbers
FBFM13_CODES = list(range(1, 14))

# Scott & Burgan (2005) FBFM40 codes (NB, GR, GS, SH, TU, TL, SB series)
FBFM40_CODES = (
    list(range(101, 110)) +   # NB1–NB9
    list(range(121, 130)) +   # GR1–GR9
    list(range(141, 146)) +   # GS1–GS4 (141–144) + GS5 placeholder
    list(range(161, 166)) +   # SH1–SH9 (partial)
    [181, 182, 183, 184, 185, 186, 187, 188, 189] +  # TU1–TU5 (partial)
    [201, 202, 203, 204, 205, 206, 207, 208, 209] +  # TL1–TL9
    [221, 222, 223, 224]                               # SB1–SB4
)


def generate_adj_template(models: list[int], filepath: str,
                           comment: str = "") -> None:
    """
    Write a template .adj file with all adj_factors set to 1.0.

    Parameters
    ----------
    models : list[int]
        Fuel model codes to include.
    filepath : str
        Output file path.
    comment : str, optional
        Extra comment line to include in the header.
    """
    path = Path(filepath)
    with open(path, "w") as f:
        f.write("# FARSITE fuel adjustment file (.adj)\n")
        f.write("# Format: <fuel_model_number>  <adj_factor>\n")
        f.write("# adj_factor = multiplicative scale on Rothermel R0 (1.0 = no change)\n")
        if comment:
            f.write(f"# {comment}\n")
        f.write(f"{len(models)}\n")
        for m in sorted(models):
            f.write(f"{m:4d}  1.0000\n")
    print(f"Wrote template adj file with {len(models)} models → {path}")


# ---------------------------------------------------------------------------
# Apply to CSV (optional calibration helper)
# ---------------------------------------------------------------------------

def apply_adj_to_csv(adj: dict[int, float], in_csv: str, out_csv: str) -> None:
    """
    Read a Rothermel parameter CSV (with a 'fuel_model' and 'w0' column) and
    multiply w0 (total oven-dry fuel load) by the adj_factor for each row.

    This mirrors what the C++ solver does internally:  scaling w0 propagates
    through reaction intensity and R0 consistently.

    Parameters
    ----------
    adj : dict[int, float]
        Fuel model → adj_factor mapping from :func:`parse_adj_file`.
    in_csv : str
        Input CSV path.
    out_csv : str
        Output CSV path.
    """
    in_path  = Path(in_csv)
    out_path = Path(out_csv)

    if not in_path.exists():
        sys.exit(f"ERROR: input CSV not found: {in_csv}")

    with open(in_path, newline="") as fin, \
         open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            sys.exit("ERROR: empty CSV")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        applied = 0
        for row in reader:
            try:
                model = int(row["fuel_model"])
                if model in adj and adj[model] != 1.0:
                    for col in ("w0", "w_d1", "w_d10", "w_d100", "w_lh", "w_lw"):
                        if col in row and row[col]:
                            try:
                                row[col] = str(float(row[col]) * adj[model])
                            except ValueError:
                                pass
                    applied += 1
            except (KeyError, ValueError):
                pass
            writer.writerow(row)
    print(f"Applied adjustments to {applied} rows → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--adj", metavar="FILE",
                   help="Path to .adj file to inspect or apply")
    p.add_argument("--generate", action="store_true",
                   help="Generate a template .adj file (requires --output)")
    p.add_argument("--fuel-system", choices=["13", "40"], default="13",
                   metavar="{13|40}",
                   help="Fuel model system for template generation (default: 13)")
    p.add_argument("--models", nargs="+", type=int, metavar="MODEL",
                   help="Specific fuel model numbers to include in template")
    p.add_argument("--apply-csv", metavar="FILE",
                   help="Rothermel CSV to apply adjustments to (optional)")
    p.add_argument("--output", "-o", metavar="FILE",
                   help="Output file path (adj template or adjusted CSV)")
    return p


def main() -> None:
    args = _build_parser().parse_args()

    if args.generate:
        if not args.output:
            sys.exit("ERROR: --generate requires --output")
        models = (args.models if args.models
                  else (FBFM13_CODES if args.fuel_system == "13"
                        else FBFM40_CODES))
        generate_adj_template(
            models, args.output,
            comment=f"FBFM{args.fuel_system} template — set adj_factor to tune ROS per model",
        )
        return

    if args.adj:
        adj = parse_adj_file(args.adj)
        print(f"\nFuel adjustments from: {args.adj}")
        print_adj_table(adj)

        if args.apply_csv:
            out = args.output or "rothermel_adjusted.csv"
            apply_adj_to_csv(adj, args.apply_csv, out)
        return

    # No action specified
    _build_parser().print_help()


if __name__ == "__main__":
    main()
