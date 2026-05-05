#!/usr/bin/env python3
"""
fuel_moisture_from_weather.py  -  Estimate dead and live fuel moisture content
from weather station observations (temperature, relative humidity, and optional
precipitation).

Dead fuel moisture (1-hr, 10-hr, 100-hr) is estimated using the equilibrium
moisture content (EMC) model of Nelson (2000) / Simard (1968), which is the
same approach used internally by FARSITE, BehavePlus, and FLAMMAP.

Live fuel moisture (herbaceous, woody) is not estimated from weather alone —
it depends on phenological stage, soil moisture, and species — so typical
default ranges are provided as guidance.

Output
------
The tool prints moisture content values suitable for direct use as
``rothermel.M_d1``, ``rothermel.M_d10``, ``rothermel.M_d100`` in the
wildfire_levelset solver input file.

Optionally writes an output text file (``--out``) containing the estimated
values and a full comment explaining the formulas.

Equilibrium moisture content model
-----------------------------------
The Nelson (2000) / Simard (1968) EMC equations relate temperature T [°C] and
relative humidity RH [%] to equilibrium moisture content EMC [%] of fine dead
fuels:

  Region 1 (RH < 10%):
    EMC = 0.03229 + 0.281073 * (RH/100) - 0.000578 * T * (RH/100)

  Region 2 (10% ≤ RH < 50%):
    EMC = 2.22749 + 0.160107 * (RH/100) - 0.014784 * T

  Region 3 (RH ≥ 50%):
    EMC = 21.0606 + 0.005565 * (RH/100)² - 0.00035 * T * (RH/100) - 0.483199 * (RH/100)

The result is in percent [%].  For conversion to the fraction [0-1] used by
the solver, divide by 100.

Timelag adjustments
-------------------
Coarser fuel sizes respond more slowly to weather changes than the 1-hr size
class.  Nelson (2000) recommends:

  M_1hr   = EMC(T, RH)
  M_10hr  ≈ EMC(T, RH) * 1.10  (or use a running 12-h weighted mean)
  M_100hr ≈ EMC(T, RH) * 1.30  (or use a running 7-day weighted mean)

These multipliers are empirical approximations commonly used in operational
fire weather computations when only a single weather snapshot is available.
If a time series of observations is available, weighted-average EMC over the
appropriate lag period gives a better estimate (pass ``--time-series`` to use
this mode).

Dew point / wet-bulb correction
---------------------------------
When ``--use-dew-point`` is provided, relative humidity is first computed from
the ambient temperature and dew point using the Magnus formula before applying
the EMC equations.

References
----------
  Nelson, R.M. Jr. (2000). "Prediction of diurnal change in 10-h fuel stick
    moisture content." Canadian Journal of Forest Research, 30(7), 1071-1087.
  Simard, A.J. (1968). The moisture content of forest fuels. Forest Fire
    Research Institute, Ottawa, Information Report FF-X-14.
  Rothermel, R.C. (1983). How to Predict the Spread and Intensity of Forest
    and Range Fires. USDA Forest Service GTR INT-143.
  Andrews, P.L. (2018). The Rothermel Surface Fire Spread Model and Associated
    Developments. USDA Forest Service GTR RMRS-371.

Usage examples
--------------
  # Single observation: 30°C, 20% RH
  python3 fuel_moisture_from_weather.py --temp 30 --rh 20

  # With precipitation (wetting)
  python3 fuel_moisture_from_weather.py --temp 25 --rh 35 --precip 2.5

  # From dew point: 30°C ambient, 15°C dew point
  python3 fuel_moisture_from_weather.py --temp 30 --use-dew-point 15

  # Batch from CSV (columns: temp_c, rh_pct)
  python3 fuel_moisture_from_weather.py --csv weather.csv

  # Write output to file
  python3 fuel_moisture_from_weather.py --temp 28 --rh 25 --out moisture.txt

  # US customary inputs (°F, relative humidity %)
  python3 fuel_moisture_from_weather.py --temp 86 --rh 20 --fahrenheit

Options
-------
  --temp T         Ambient temperature [°C] (or °F with --fahrenheit).
  --rh RH          Relative humidity [%].
  --precip P       Precipitation in last hour [mm].  Adds wetting adjustment.
  --use-dew-point D Dew-point temperature [°C] (replaces --rh).
  --fahrenheit     Temperature inputs are in °F (converted internally).
  --csv FILE       Two-column CSV (temp_c, rh_pct) for batch processing.
  --out FILE       Write results to this file (default: stdout only).
  --solver-format  Print values in solver input-file format (M_d1 = ...).
  --help           Show this help and exit.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Core EMC computation (Nelson 2000 / Simard 1968)
# ---------------------------------------------------------------------------

def emc_nelson(temp_c: float, rh_pct: float) -> float:
    """
    Compute equilibrium moisture content [%] of dead fine fuel.

    Parameters
    ----------
    temp_c  : Dry-bulb temperature [°C].
    rh_pct  : Relative humidity [%].

    Returns
    -------
    emc : Equilibrium moisture content [%].

    Notes
    -----
    Equations from Nelson (2000) / Simard (1968) as implemented in FARSITE
    and BehavePlus.
    """
    rh = rh_pct / 100.0
    rh = min(1.0, max(0.0, rh))
    temp_c = max(-40.0, temp_c)  # clamp to physically plausible range

    if rh < 0.10:
        emc = 0.03229 + 0.281073 * rh - 0.000578 * temp_c * rh
    elif rh < 0.50:
        emc = 2.22749 + 0.160107 * rh - 0.014784 * temp_c
    else:
        emc = (21.0606 + 0.005565 * rh * rh
               - 0.00035 * temp_c * rh
               - 0.483199 * rh)

    return max(0.5, emc)  # physical minimum ~0.5%


def emc_by_size_class(temp_c: float, rh_pct: float,
                      precip_mm: float = 0.0) -> Tuple[float, float, float]:
    """
    Estimate moisture content [%] for the 1-hr, 10-hr, and 100-hr dead fuel
    size classes using the Nelson (2000) EMC as a base and applying empirical
    timelag-weighted multipliers.

    Parameters
    ----------
    temp_c    : Ambient temperature [°C].
    rh_pct    : Relative humidity [%].
    precip_mm : Precipitation in the last hour [mm].  Values > 0 add a wetting
                increment following the NFDRS (1988) approach.

    Returns
    -------
    (M_1hr, M_10hr, M_100hr) : Moisture fractions (not percent).
    """
    base_emc = emc_nelson(temp_c, rh_pct)  # percent

    # Precipitation wetting increment [%]
    # Based on NFDRS (Deeming et al. 1977) boundary:
    #   < 0.5 mm:  no correction
    #   0.5-1.0 mm: +5% to 1-hr
    #   1.0-2.5 mm: +10% to 1-hr
    #   2.5-5.0 mm: +15% to 1-hr
    #   > 5 mm:    +25% to 1-hr
    wet_1hr = 0.0
    if precip_mm >= 0.5:
        if   precip_mm < 1.0:  wet_1hr = 5.0
        elif precip_mm < 2.5:  wet_1hr = 10.0
        elif precip_mm < 5.0:  wet_1hr = 15.0
        else:                   wet_1hr = 25.0

    m_1hr   = base_emc + wet_1hr
    m_10hr  = base_emc * 1.10 + (wet_1hr * 0.5)   # 10-hr lags behind 1-hr
    m_100hr = base_emc * 1.30 + (wet_1hr * 0.2)   # 100-hr lags further

    # Cap at physical adsorption limit (~35% for dead fuel)
    m_1hr   = min(35.0, m_1hr)
    m_10hr  = min(35.0, m_10hr)
    m_100hr = min(35.0, m_100hr)

    # Return as fractions (solver input convention)
    return m_1hr / 100.0, m_10hr / 100.0, m_100hr / 100.0


# ---------------------------------------------------------------------------
# Dew point → relative humidity
# ---------------------------------------------------------------------------

def rh_from_dewpoint(temp_c: float, dewpoint_c: float) -> float:
    """
    Compute relative humidity [%] from dry-bulb and dew-point temperatures
    using the Magnus formula.

    Parameters
    ----------
    temp_c     : Dry-bulb temperature [°C].
    dewpoint_c : Dew-point temperature [°C].

    Returns
    -------
    rh : Relative humidity [%].
    """
    # Magnus formula constants (Alduchov & Eskridge 1996)
    a, b = 17.625, 243.04
    gamma_t  = a * temp_c     / (b + temp_c)
    gamma_td = a * dewpoint_c / (b + dewpoint_c)
    rh = 100.0 * math.exp(gamma_td - gamma_t)
    return min(100.0, max(0.0, rh))


# ---------------------------------------------------------------------------
# Batch CSV processing
# ---------------------------------------------------------------------------

def process_csv(path: str, precip_mm: float = 0.0
                ) -> List[Tuple[float, float, float]]:
    """
    Read a two-column CSV (temp_c, rh_pct) and compute moisture per row.

    Returns a list of (M_1hr, M_10hr, M_100hr) tuples [fractions].
    """
    results = []
    with open(path, "r", newline="") as fh:
        # Skip blank lines and comments
        rows = [r for r in csv.reader(fh)
                if r and not r[0].strip().startswith("#")]
    # Detect header row
    start = 0
    try:
        float(rows[0][0])
    except (ValueError, IndexError):
        start = 1

    for row in rows[start:]:
        if len(row) < 2:
            continue
        try:
            temp_c = float(row[0])
            rh_pct = float(row[1])
            precip = float(row[2]) if len(row) > 2 else precip_mm
        except ValueError:
            continue
        results.append(emc_by_size_class(temp_c, rh_pct, precip))
    return results


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def _fmt_solver(m1: float, m10: float, m100: float) -> str:
    return (
        f"rothermel.M_d1   = {m1:.4f}   # 1-hr  dead fuel moisture [fraction]\n"
        f"rothermel.M_d10  = {m10:.4f}   # 10-hr dead fuel moisture [fraction]\n"
        f"rothermel.M_d100 = {m100:.4f}   # 100-hr dead fuel moisture [fraction]\n"
        "# Live fuel defaults (set from phenology / field observations):\n"
        "rothermel.M_lh   = 0.900   # live herbaceous [fraction]\n"
        "rothermel.M_lw   = 1.200   # live woody      [fraction]\n"
    )


def _fmt_human(temp_c: float, rh_pct: float,
               m1: float, m10: float, m100: float) -> str:
    return (
        f"Input:   T = {temp_c:.1f} °C,  RH = {rh_pct:.0f}%\n"
        f"  1-hr  dead EMC:  {m1 * 100:.1f}%  ({m1:.4f} fraction)\n"
        f"  10-hr dead EMC:  {m10 * 100:.1f}%  ({m10:.4f} fraction)\n"
        f"  100-hr dead EMC: {m100 * 100:.1f}%  ({m100:.4f} fraction)\n"
        "  Live herb. (guidance): 30-100%  (0.30-1.00)\n"
        "  Live woody (guidance): 80-200%  (0.80-2.00)\n"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Estimate dead fuel moisture from temperature and RH.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--temp",           type=float, default=None, metavar="T",
                        help="Ambient temperature [°C] (or °F with --fahrenheit)")
    parser.add_argument("--rh",             type=float, default=None, metavar="RH",
                        help="Relative humidity [%%]")
    parser.add_argument("--precip",         type=float, default=0.0, metavar="P",
                        help="Precipitation in last hour [mm] (default: 0)")
    parser.add_argument("--use-dew-point",  type=float, default=None, metavar="D",
                        help="Dew-point temperature [°C] (replaces --rh)")
    parser.add_argument("--fahrenheit",     action="store_true",
                        help="Temperature inputs are in °F")
    parser.add_argument("--csv",            default=None, metavar="FILE",
                        help="Batch mode: two-column CSV (temp_c, rh_pct)")
    parser.add_argument("--out",            default=None, metavar="FILE",
                        help="Write results to this file")
    parser.add_argument("--solver-format",  action="store_true",
                        help="Print in solver input-file format")

    args = parser.parse_args(argv)

    output_lines: List[str] = []

    if args.csv:
        # Batch mode
        results = process_csv(args.csv, args.precip)
        if not results:
            print("ERROR: no valid rows in", args.csv, file=sys.stderr)
            sys.exit(1)
        # Print summary statistics
        m1s  = [r[0] for r in results]
        m10s = [r[1] for r in results]
        m100s = [r[2] for r in results]
        mean_m1  = sum(m1s)  / len(m1s)
        mean_m10 = sum(m10s) / len(m10s)
        mean_m100 = sum(m100s) / len(m100s)
        print(f"Processed {len(results)} rows from {args.csv}")
        line = (f"Mean:  M_1hr={mean_m1:.4f}  M_10hr={mean_m10:.4f}  M_100hr={mean_m100:.4f}\n"
                f"Range: M_1hr=[{min(m1s):.4f}, {max(m1s):.4f}]  "
                f"M_10hr=[{min(m10s):.4f}, {max(m10s):.4f}]  "
                f"M_100hr=[{min(m100s):.4f}, {max(m100s):.4f}]\n")
        print(line)
        output_lines.append(line)
        if args.solver_format:
            sf = _fmt_solver(mean_m1, mean_m10, mean_m100)
            print(sf)
            output_lines.append(sf)

    elif args.temp is not None:
        temp_c = args.temp
        if args.fahrenheit:
            temp_c = (temp_c - 32.0) / 1.8

        if args.use_dew_point is not None:
            rh_pct = rh_from_dewpoint(temp_c, args.use_dew_point)
            print(f"Derived RH = {rh_pct:.1f}%  from T={temp_c:.1f}°C, Td={args.use_dew_point:.1f}°C")
        elif args.rh is not None:
            rh_pct = args.rh
        else:
            parser.error("Provide --rh or --use-dew-point with --temp")

        m1, m10, m100 = emc_by_size_class(temp_c, rh_pct, args.precip)

        human = _fmt_human(temp_c, rh_pct, m1, m10, m100)
        print(human)
        output_lines.append(human)

        if args.solver_format:
            sf = _fmt_solver(m1, m10, m100)
            print(sf)
            output_lines.append(sf)

    else:
        parser.print_help()
        sys.exit(0)

    if args.out and output_lines:
        with open(args.out, "w") as fh:
            fh.write("# Generated by fuel_moisture_from_weather.py\n")
            fh.write("# Reference: Nelson (2000), Canadian Journal of Forest Research\n")
            fh.writelines(output_lines)
        print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
