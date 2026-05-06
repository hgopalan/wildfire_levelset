#!/usr/bin/env python3
"""
farsite_fms_reader.py  -  Parse FARSITE Fuel Moisture Scenario (.fms) files
and write solver-compatible moisture input values for wildfire_levelset.

A FARSITE fuel moisture scenario (.fms) file specifies per-fuel-model dead and
live fuel moisture contents:

    <num_fuel_models>
    <fuel_model_num>  <1hr%>  <10hr%>  <100hr%>  <live_herb%>  <live_wood%>
    ...

Moisture values are given as integer percentages (e.g., 8 means 8%).

Output modes
------------
1. **Print** (default): writes a summary table and a ready-to-paste block
   of wildfire_levelset solver input parameters for each fuel model.

2. **--write-inputs FILE**: writes an inputs.i stub with the FIRST fuel
   model's moisture values (or the model specified with --fuel-model).

3. **--write-csv FILE**: writes a CSV table with all fuel models and
   their moisture values for further processing.

4. **--query FUEL_NUM**: print only the entry for a specific fuel model.

Usage examples
--------------
  # Print moisture summary from an .fms file
  python3 farsite_fms_reader.py fire.fms

  # Write solver input stub using fuel model 4
  python3 farsite_fms_reader.py fire.fms --fuel-model 4 --write-inputs inputs_moisture.i

  # Export all moisture values as CSV
  python3 farsite_fms_reader.py fire.fms --write-csv moisture.csv

  # Query a single fuel model
  python3 farsite_fms_reader.py fire.fms --query 8

File format notes
-----------------
  • Lines beginning with # or ! are treated as comments.
  • Moisture values are percentages; the solver uses fractions (percent / 100).
  • FARSITE uses 13 (FBFM13) or 40 (FBFM40) standard fuel model numbers.
  • Fuel models 91-99 (non-burnable) are silently skipped.
  • Time-varying FMS data (multiple scenario blocks) is not supported by
    this format; use the FARSITE .wtr reader for time-varying moisture.

References
----------
  Finney, M.A. (2004). FARSITE: Fire Area Simulator — Model Development
    and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
  Scott, J.H. & Burgan, R.E. (2005). Standard Fire Behavior Fuel Models:
    A Comprehensive Set for Use with Rothermel's Surface Fire Spread Model.
    USDA Forest Service RMRS-GTR-153.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


# ---------------------------------------------------------------------------
# Non-burnable fuel model codes (FBFM13 and FBFM40)
# ---------------------------------------------------------------------------
NON_BURNABLE = {91, 92, 93, 98, 99}


# ---------------------------------------------------------------------------
# Data structure
# ---------------------------------------------------------------------------

@dataclass
class FuelMoisture:
    """Moisture values for a single fuel model."""
    fuel_model: int
    m_1hr: float    # 1-hr dead fuel moisture  [fraction]
    m_10hr: float   # 10-hr dead fuel moisture [fraction]
    m_100hr: float  # 100-hr dead fuel moisture [fraction]
    m_lh: float     # live herbaceous moisture  [fraction]
    m_lw: float     # live woody moisture       [fraction]


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_fms(path: str) -> List[FuelMoisture]:
    """
    Parse a FARSITE .fms fuel moisture scenario file.

    Parameters
    ----------
    path : str
        Path to the .fms file.

    Returns
    -------
    List[FuelMoisture]
        One entry per fuel model found in the file.  Fuel model codes in
        ``NON_BURNABLE`` are silently skipped.

    Raises
    ------
    ValueError
        If the file is malformed (missing header count, wrong column count).
    """
    records: List[FuelMoisture] = []

    with open(path, "r") as fh:
        lines = [ln.strip() for ln in fh
                 if ln.strip() and not ln.strip().startswith(("#", "!"))]

    if not lines:
        raise ValueError(f"No data found in {path}")

    # First non-comment line is the number of fuel models
    try:
        n_models = int(lines[0])
    except ValueError:
        raise ValueError(
            f"Expected integer fuel model count on first line, got: {lines[0]!r}"
        )

    data_lines = lines[1:]
    if len(data_lines) < n_models:
        # Tolerate files that list fewer rows than the declared count
        import warnings
        warnings.warn(
            f"File declares {n_models} fuel models but only {len(data_lines)} "
            f"data rows found.  Processing available rows."
        )

    for ln in data_lines:
        parts = ln.split()
        if len(parts) < 6:
            continue  # skip malformed or short lines
        try:
            fm  = int(parts[0])
            m1  = float(parts[1]) / 100.0
            m10 = float(parts[2]) / 100.0
            m100 = float(parts[3]) / 100.0
            mlh  = float(parts[4]) / 100.0
            mlw  = float(parts[5]) / 100.0
        except ValueError:
            continue

        if fm in NON_BURNABLE:
            continue

        records.append(FuelMoisture(
            fuel_model=fm,
            m_1hr=m1,
            m_10hr=m10,
            m_100hr=m100,
            m_lh=mlh,
            m_lw=mlw,
        ))

    return records


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _solver_block(fm: FuelMoisture, comment: str = "") -> str:
    """Return a solver inputs.i block for one fuel model."""
    lines = []
    if comment:
        lines.append(f"# {comment}")
    lines += [
        f"# Fuel model {fm.fuel_model} moisture (from .fms scenario)",
        f"rothermel.M_d1   = {fm.m_1hr:.4f}   # 1-hr  dead moisture [fraction]",
        f"rothermel.M_d10  = {fm.m_10hr:.4f}   # 10-hr dead moisture [fraction]",
        f"rothermel.M_d100 = {fm.m_100hr:.4f}   # 100-hr dead moisture [fraction]",
        f"rothermel.M_lh   = {fm.m_lh:.4f}   # live herbaceous moisture [fraction]",
        f"rothermel.M_lw   = {fm.m_lw:.4f}   # live woody moisture [fraction]",
    ]
    return "\n".join(lines)


def _table_row(fm: FuelMoisture) -> str:
    return (
        f"  FM {fm.fuel_model:3d}  "
        f"1hr={fm.m_1hr*100:5.1f}%  "
        f"10hr={fm.m_10hr*100:5.1f}%  "
        f"100hr={fm.m_100hr*100:5.1f}%  "
        f"lh={fm.m_lh*100:5.1f}%  "
        f"lw={fm.m_lw*100:5.1f}%"
    )


# ---------------------------------------------------------------------------
# Output functions
# ---------------------------------------------------------------------------

def print_summary(records: List[FuelMoisture]) -> None:
    print(f"Parsed {len(records)} fuel model moisture record(s):\n")
    print(f"  {'FM':>4}  {'1hr%':>7}  {'10hr%':>7}  {'100hr%':>8}  "
          f"{'LH%':>7}  {'LW%':>7}")
    print("  " + "-" * 55)
    for fm in records:
        print(_table_row(fm))
    print()
    print("Solver input blocks (paste into your inputs.i file):\n")
    for fm in records:
        print(_solver_block(fm))
        print()


def write_inputs_stub(path: str, fm: FuelMoisture) -> None:
    """Write a minimal solver inputs.i stub for one fuel model."""
    content = (
        "# Auto-generated by farsite_fms_reader.py\n"
        + _solver_block(fm, f"Generated from .fms file for fuel model {fm.fuel_model}")
        + "\n"
    )
    with open(path, "w") as fh:
        fh.write(content)
    print(f"Wrote solver moisture stub: {path}")


def write_csv(path: str, records: List[FuelMoisture]) -> None:
    """Write all fuel model moisture values as CSV."""
    fieldnames = [
        "fuel_model", "m_1hr_frac", "m_10hr_frac", "m_100hr_frac",
        "m_lh_frac", "m_lw_frac",
        "m_1hr_pct", "m_10hr_pct", "m_100hr_pct",
        "m_lh_pct", "m_lw_pct",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for fm in records:
            writer.writerow({
                "fuel_model":   fm.fuel_model,
                "m_1hr_frac":   f"{fm.m_1hr:.4f}",
                "m_10hr_frac":  f"{fm.m_10hr:.4f}",
                "m_100hr_frac": f"{fm.m_100hr:.4f}",
                "m_lh_frac":    f"{fm.m_lh:.4f}",
                "m_lw_frac":    f"{fm.m_lw:.4f}",
                "m_1hr_pct":    f"{fm.m_1hr*100:.1f}",
                "m_10hr_pct":   f"{fm.m_10hr*100:.1f}",
                "m_100hr_pct":  f"{fm.m_100hr*100:.1f}",
                "m_lh_pct":     f"{fm.m_lh*100:.1f}",
                "m_lw_pct":     f"{fm.m_lw*100:.1f}",
            })
    print(f"Wrote CSV: {path}  ({len(records)} rows)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Parse a FARSITE .fms fuel moisture scenario file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "fms_file",
        metavar="FILE",
        help="Input FARSITE .fms file.",
    )
    parser.add_argument(
        "--fuel-model", "-f",
        type=int, default=None, metavar="NUM",
        help="Fuel model number to use for --write-inputs (default: first in file).",
    )
    parser.add_argument(
        "--write-inputs",
        default=None, metavar="FILE",
        help="Write a solver inputs.i moisture stub to FILE.",
    )
    parser.add_argument(
        "--write-csv",
        default=None, metavar="FILE",
        help="Write all moisture values as a CSV table to FILE.",
    )
    parser.add_argument(
        "--query", "-q",
        type=int, default=None, metavar="NUM",
        help="Print solver block for a specific fuel model number.",
    )

    args = parser.parse_args(argv)

    # Parse the FMS file
    try:
        records = parse_fms(args.fms_file)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if not records:
        print("ERROR: no burnable fuel model records found in the file.",
              file=sys.stderr)
        sys.exit(1)

    # Build lookup dict
    by_fm: Dict[int, FuelMoisture] = {r.fuel_model: r for r in records}

    # --query: print a single fuel model's block
    if args.query is not None:
        if args.query not in by_fm:
            print(f"ERROR: fuel model {args.query} not found in file.", file=sys.stderr)
            sys.exit(1)
        print(_solver_block(by_fm[args.query]))
        return

    # Default: print full summary
    print_summary(records)

    # --write-inputs
    if args.write_inputs is not None:
        if args.fuel_model is not None:
            if args.fuel_model not in by_fm:
                print(f"ERROR: fuel model {args.fuel_model} not in file.",
                      file=sys.stderr)
                sys.exit(1)
            target = by_fm[args.fuel_model]
        else:
            target = records[0]
        write_inputs_stub(args.write_inputs, target)

    # --write-csv
    if args.write_csv is not None:
        write_csv(args.write_csv, records)


if __name__ == "__main__":
    main()
