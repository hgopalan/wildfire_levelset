#!/usr/bin/env python3
"""
farsite_fmd_reader.py – Read, inspect, and convert FARSITE fuel moisture
schedule (`.fmd`) files.

A FARSITE ``.fmd`` file carries a time-stamped schedule of fuel moisture
contents by fuel model.  This script can:

  1. **Inspect** an existing ``.fmd`` file (print the schedule as a table or
     CSV to stdout).
  2. **Convert** the schedule to a flat CSV that the solver can use for
     time-varying moisture (one row per timestamp × fuel model).
  3. **Interpolate** and **query** the schedule at an arbitrary simulation
     time so you can verify what moisture values the solver will use at any
     given time.
  4. **Generate** a constant-moisture ``.fmd`` template for a given time
     range and set of fuel models.

FMD file format
---------------
Each *timestamp block* has a header line followed by one data row per fuel
model::

    MONTH  DAY  HOUR  PRECIP  NUM_MODELS
    MODEL  1HR%  10HR%  100HR%  LHERB%  LWOOD%
    ...

- ``HOUR`` is in HHMM format (e.g. ``1400`` = 14:00 local time).
- ``PRECIP`` is precipitation in 100ths of an inch (stored but not used by
  this script).
- Moisture values are **integer percentages** (e.g. ``8`` → 8 %, which is
  0.08 in fraction form).
- A ``MODEL`` value of ``0`` means "apply to all fuel models".
- Lines starting with ``#`` or ``!`` are comments.

Output CSV columns
------------------
``month,day,hour,model,M_d1_pct,M_d10_pct,M_d100_pct,M_lh_pct,M_lw_pct``

**References**
    Finney, M.A. (2004). *FARSITE: Fire Area Simulator—Model Development and
    Evaluation*. USDA Forest Service Research Paper RMRS-RP-4.

Usage examples
--------------
Inspect an FMD file::

    python3 farsite_fmd_reader.py --fmd fire.fmd

Convert to CSV (for use with the solver or external analysis)::

    python3 farsite_fmd_reader.py --fmd fire.fmd --output moisture_schedule.csv

Query moisture at a specific simulation time (seconds since first record)::

    python3 farsite_fmd_reader.py --fmd fire.fmd \\
        --query-time 3600 --fuel-model 4

Generate a constant-moisture template covering 24 h at 1-h intervals::

    python3 farsite_fmd_reader.py --generate \\
        --models 4 9 10 \\
        --start-month 7 --start-day 15 --hours 24 \\
        --M-d1 8 --M-d10 10 --M-d100 15 \\
        --M-lh 90 --M-lw 120 \\
        --output template.fmd
"""

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MoistureRecord:
    """Fuel moisture values at one (timestamp, model) pair."""
    month: int
    day: int
    hour: int        # HHMM integer
    precip: int      # 100ths of an inch
    model: int       # fuel model code; 0 = "all models"
    M_d1: float      # 1-hr   dead fuel moisture [%]
    M_d10: float     # 10-hr  dead fuel moisture [%]
    M_d100: float    # 100-hr dead fuel moisture [%]
    M_lh: float      # live herbaceous moisture [%]
    M_lw: float      # live woody moisture [%]


@dataclass
class FmdSchedule:
    """Parsed FMD schedule."""
    records: list[MoistureRecord] = field(default_factory=list)
    source: str = ""

    @property
    def timestamps(self) -> list[tuple[int, int, int]]:
        """Unique (month, day, hour) tuples in order of appearance."""
        seen: list[tuple[int, int, int]] = []
        for r in self.records:
            t = (r.month, r.day, r.hour)
            if t not in seen:
                seen.append(t)
        return seen

    @property
    def fuel_models(self) -> list[int]:
        """Sorted list of unique fuel model codes."""
        return sorted(set(r.model for r in self.records))


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_fmd_file(filepath: str) -> FmdSchedule:
    """
    Parse a FARSITE .fmd file.

    Parameters
    ----------
    filepath : str
        Path to the ``.fmd`` file.

    Returns
    -------
    FmdSchedule
        Parsed schedule with all records.

    Raises
    ------
    SystemExit
        If the file cannot be opened or contains no valid records.
    """
    path = Path(filepath)
    if not path.exists():
        sys.exit(f"ERROR: FMD file not found: {filepath}")

    sched = FmdSchedule(source=str(path))
    lines = path.read_text().splitlines()
    idx = 0
    n = len(lines)

    def _next_data_line() -> Optional[str]:
        nonlocal idx
        while idx < n:
            line = lines[idx].strip()
            idx += 1
            if line and line[0] not in ("#", "!"):
                return line
        return None

    while idx < n:
        header = _next_data_line()
        if header is None:
            break
        tokens = header.split()
        if len(tokens) < 5:
            continue
        try:
            month  = int(tokens[0])
            day    = int(tokens[1])
            hour   = int(tokens[2])
            precip = int(tokens[3])
            nmod   = int(tokens[4])
        except ValueError:
            continue

        if month < 1 or month > 12 or day < 1 or day > 31 or nmod < 1:
            continue

        for _ in range(nmod):
            data_line = _next_data_line()
            if data_line is None:
                break
            dtok = data_line.split()
            if len(dtok) < 6:
                continue
            try:
                model  = int(dtok[0])
                m1     = float(dtok[1])
                m10    = float(dtok[2])
                m100   = float(dtok[3])
                mlh    = float(dtok[4])
                mlw    = float(dtok[5])
            except ValueError:
                continue
            sched.records.append(MoistureRecord(
                month=month, day=day, hour=hour, precip=precip,
                model=model,
                M_d1=m1, M_d10=m10, M_d100=m100, M_lh=mlh, M_lw=mlw,
            ))

    if not sched.records:
        sys.exit(f"ERROR: no valid records found in {filepath}")
    return sched


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def print_schedule(sched: FmdSchedule) -> None:
    """Print a human-readable summary of the schedule."""
    print(f"FMD schedule: {sched.source}")
    print(f"  Timestamps : {len(sched.timestamps)}")
    print(f"  Fuel models: {sched.fuel_models}")
    print(f"  Total records: {len(sched.records)}")
    print()
    # Header
    print(f"{'Month':>6} {'Day':>4} {'Hour':>6} {'Model':>6} "
          f"{'1hr%':>6} {'10hr%':>6} {'100hr%':>7} {'LHrb%':>6} {'LWdy%':>6}")
    print("-" * 57)
    for r in sched.records:
        print(f"{r.month:6d} {r.day:4d} {r.hour:6d} {r.model:6d} "
              f"{r.M_d1:6.1f} {r.M_d10:6.1f} {r.M_d100:7.1f} "
              f"{r.M_lh:6.1f} {r.M_lw:6.1f}")


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def schedule_to_csv(sched: FmdSchedule, out_path: str) -> None:
    """
    Write the schedule to a flat CSV.

    Columns: ``month, day, hour, model, M_d1_pct, M_d10_pct, M_d100_pct,
    M_lh_pct, M_lw_pct``.

    Parameters
    ----------
    sched : FmdSchedule
        Parsed schedule.
    out_path : str
        Output CSV file path.
    """
    fieldnames = [
        "month", "day", "hour", "model",
        "M_d1_pct", "M_d10_pct", "M_d100_pct", "M_lh_pct", "M_lw_pct",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in sched.records:
            w.writerow({
                "month": r.month, "day": r.day, "hour": r.hour,
                "model": r.model,
                "M_d1_pct":   r.M_d1,
                "M_d10_pct":  r.M_d10,
                "M_d100_pct": r.M_d100,
                "M_lh_pct":   r.M_lh,
                "M_lw_pct":   r.M_lw,
            })
    print(f"Wrote {len(sched.records)} records → {out_path}")


# ---------------------------------------------------------------------------
# Query / interpolation
# ---------------------------------------------------------------------------

def _hhmm_to_seconds(hhmm: int) -> float:
    """Convert HHMM integer to seconds within a day."""
    return (hhmm // 100) * 3600.0 + (hhmm % 100) * 60.0


def query_at_time(sched: FmdSchedule,
                  sim_time_s: float,
                  fuel_model: int = 0) -> MoistureRecord:
    """
    Return linearly interpolated moisture values at *sim_time_s* seconds
    since the first record's timestamp.

    Searches for records matching *fuel_model* (or model 0 as a fallback),
    brackets by simulation time, and linearly interpolates.

    Parameters
    ----------
    sched : FmdSchedule
        Parsed schedule.
    sim_time_s : float
        Simulation time in seconds (0 = first timestamp).
    fuel_model : int
        Fuel model code to look up; 0 = "all models" global entry.

    Returns
    -------
    MoistureRecord
        Interpolated record (month/day/hour/model reflect the lower bracket).
    """
    # Build list of (t_s, record) for the requested model or model=0
    candidates = [
        r for r in sched.records
        if r.model == fuel_model or r.model == 0
    ]
    if not candidates:
        sys.exit(f"ERROR: no records found for fuel model {fuel_model}")

    # Compute relative time from first candidate
    t0 = (_hhmm_to_seconds(candidates[0].hour)
          + candidates[0].day * 86400
          + candidates[0].month * 30 * 86400)
    timed = [((_hhmm_to_seconds(r.hour)
               + r.day * 86400
               + r.month * 30 * 86400) - t0, r)
             for r in candidates]
    timed.sort(key=lambda x: x[0])

    # Find brackets
    lo_t, lo_r = timed[0]
    hi_t, hi_r = timed[-1]
    for (t, r) in timed:
        if t <= sim_time_s:
            lo_t, lo_r = t, r
        elif t > sim_time_s and t < hi_t:
            hi_t, hi_r = t, r
            break

    if lo_t >= hi_t or lo_r is hi_r:
        return lo_r  # clamp

    alpha = max(0.0, min(1.0, (sim_time_s - lo_t) / (hi_t - lo_t)))
    return MoistureRecord(
        month=lo_r.month, day=lo_r.day, hour=lo_r.hour,
        precip=lo_r.precip, model=lo_r.model,
        M_d1   = lo_r.M_d1   + alpha * (hi_r.M_d1   - lo_r.M_d1),
        M_d10  = lo_r.M_d10  + alpha * (hi_r.M_d10  - lo_r.M_d10),
        M_d100 = lo_r.M_d100 + alpha * (hi_r.M_d100 - lo_r.M_d100),
        M_lh   = lo_r.M_lh   + alpha * (hi_r.M_lh   - lo_r.M_lh),
        M_lw   = lo_r.M_lw   + alpha * (hi_r.M_lw   - lo_r.M_lw),
    )


# ---------------------------------------------------------------------------
# Template generation
# ---------------------------------------------------------------------------

def generate_fmd_template(
    models: list[int],
    start_month: int,
    start_day: int,
    n_hours: int,
    M_d1: float,
    M_d10: float,
    M_d100: float,
    M_lh: float,
    M_lw: float,
    precip: int,
    filepath: str,
) -> None:
    """
    Write a constant-moisture FMD template covering *n_hours* hours at
    hourly intervals starting from (start_month, start_day, 0000).

    Parameters
    ----------
    models : list[int]
        Fuel model codes to include.
    start_month, start_day : int
        Reference start date.
    n_hours : int
        Number of hourly snapshots to write.
    M_d1, M_d10, M_d100, M_lh, M_lw : float
        Constant moisture values [%].
    precip : int
        Precipitation (100ths of an inch; default 0).
    filepath : str
        Output ``.fmd`` file path.
    """
    path = Path(filepath)
    nmod = len(models)
    with open(path, "w") as f:
        f.write("# FARSITE fuel moisture schedule (.fmd)\n")
        f.write(f"# Constant moisture template: {nmod} models, {n_hours} hours\n")
        for hour_idx in range(n_hours):
            hh = hour_idx % 24
            dd_offset = hour_idx // 24
            day   = start_day + dd_offset
            month = start_month
            # Naive day-rollover (no true calendar)
            while day > 28:
                day   -= 28
                month += 1
            if month > 12:
                month = 1
            hhmm = hh * 100
            f.write(f"{month} {day} {hhmm:04d} {precip} {nmod}\n")
            for m in models:
                f.write(f"{m} {M_d1:.1f} {M_d10:.1f} "
                        f"{M_d100:.1f} {M_lh:.1f} {M_lw:.1f}\n")
    print(f"Wrote {n_hours} × {nmod}-model FMD template → {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--fmd", metavar="FILE",
                   help="Path to .fmd file to parse/inspect/convert")
    p.add_argument("--output", "-o", metavar="FILE",
                   help="Output CSV or .fmd file path")
    p.add_argument("--query-time", type=float, metavar="SECONDS",
                   help="Query interpolated moisture at this simulation time [s]")
    p.add_argument("--fuel-model", type=int, default=0, metavar="N",
                   help="Fuel model code for --query-time (0 = global, default: 0)")

    gen = p.add_argument_group("template generation (requires --output)")
    gen.add_argument("--generate", action="store_true",
                     help="Generate a constant-moisture .fmd template")
    gen.add_argument("--models", nargs="+", type=int, metavar="N",
                     default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                     help="Fuel model codes (default: FBFM13 1–13)")
    gen.add_argument("--start-month", type=int, default=7, metavar="M",
                     help="Start month (1–12, default: 7)")
    gen.add_argument("--start-day", type=int, default=1, metavar="D",
                     help="Start day (1–31, default: 1)")
    gen.add_argument("--hours", type=int, default=24, metavar="N",
                     help="Number of hourly snapshots (default: 24)")
    gen.add_argument("--M-d1",   type=float, default=8.0)
    gen.add_argument("--M-d10",  type=float, default=10.0)
    gen.add_argument("--M-d100", type=float, default=15.0)
    gen.add_argument("--M-lh",   type=float, default=90.0)
    gen.add_argument("--M-lw",   type=float, default=120.0)
    gen.add_argument("--precip", type=int, default=0,
                     help="Precipitation in 100ths of an inch (default: 0)")
    return p


def main() -> None:
    args = _build_parser().parse_args()

    if args.generate:
        if not args.output:
            sys.exit("ERROR: --generate requires --output")
        generate_fmd_template(
            models=args.models,
            start_month=args.start_month,
            start_day=args.start_day,
            n_hours=args.hours,
            M_d1=args.M_d1,
            M_d10=args.M_d10,
            M_d100=args.M_d100,
            M_lh=args.M_lh,
            M_lw=args.M_lw,
            precip=args.precip,
            filepath=args.output,
        )
        return

    if args.fmd:
        sched = parse_fmd_file(args.fmd)

        if args.query_time is not None:
            r = query_at_time(sched, args.query_time, args.fuel_model)
            print(f"\nInterpolated moisture at t = {args.query_time:.1f} s"
                  f"  (model {args.fuel_model}):")
            print(f"  1-hr dead:          {r.M_d1:.2f} %  ({r.M_d1/100:.4f} fraction)")
            print(f"  10-hr dead:         {r.M_d10:.2f} %")
            print(f"  100-hr dead:        {r.M_d100:.2f} %")
            print(f"  Live herbaceous:    {r.M_lh:.2f} %")
            print(f"  Live woody:         {r.M_lw:.2f} %")
            return

        print_schedule(sched)

        if args.output:
            schedule_to_csv(sched, args.output)
        return

    _build_parser().print_help()


if __name__ == "__main__":
    main()
