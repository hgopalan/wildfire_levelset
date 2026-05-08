#!/usr/bin/env python3
"""
ignition_probability_table.py – Anderson (1970) / Rothermel (1983) probability
of ignition (P_i) worksheet tool.

Computes the probability that a firebrand landing on fine dead fuel will cause
a sustained ignition, swept across a grid of fine-fuel temperatures and 1-hr
dead fuel moisture contents.  The result is printed as an ASCII table and
written to CSV — analogous to the *Probability of Ignition* worksheet in
BehavePlus.

Background
----------
The probability-of-ignition formula (Anderson 1970, as reported in
Rothermel 1983) is:

    P_i = min(1, max(0, 0.000048 × T_fuel_F^1.4 × exp(−0.07 × MC%)))

where:
    T_fuel_F  – fine-fuel temperature [°F] = T_ambient_F + solar_increment_F
    MC%       – 1-hr dead fuel moisture content [%]

The formula gives the fraction of ignition attempts (firebrands landing on
fine fuel) that result in sustained fire.

Sweeping both dimensions produces a full ignition probability matrix useful
for calibrating spotting models, scheduling prescribed burns, and assessing
wildfire risk under various weather and fuel moisture conditions.

Requirements
------------
  None (pure Python).  Optional ``matplotlib`` for heatmap plots.

Usage
-----
  # Default: ambient 20–45°C (10 °F steps), moisture 2–20% (9 steps)
  python3 tools/ignition_probability_table.py

  # Custom temperature range
  python3 tools/ignition_probability_table.py \\
      --temp-min 20 --temp-max 50 --temp-steps 7 \\
      --moisture-min 2 --moisture-max 25 --moisture-steps 10

  # Add a solar heating increment (default: 25 °F)
  python3 tools/ignition_probability_table.py --solar-heating 30

  # Save as CSV
  python3 tools/ignition_probability_table.py --out ignition_prob.csv

  # Save heatmap PNG
  python3 tools/ignition_probability_table.py --plot

References
----------
  Anderson, H.E. (1970). Forest fuel ignitibility.
    Fire Technology, 6(4), 312–319.
  Rothermel, R.C. (1983). How to predict the spread and intensity of forest
    and range fires.  USDA Forest Service General Technical Report INT-143.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional


# ---------------------------------------------------------------------------
# Anderson (1970) probability-of-ignition constants
# ---------------------------------------------------------------------------
# Anderson, H.E. (1970). Forest fuel ignitibility.
#   Fire Technology, 6(4), 312–319.  Eq. 4 (as reported in Rothermel 1983).
#
#   P_i = min(1, max(0, _PI_COEFF * T_fuel_F^_PI_TEMP_EXP * exp(_PI_MC_COEFF * MC%)))
_PI_COEFF      = 0.000048   # empirical coefficient  (Anderson 1970, Eq. 4)
_PI_TEMP_EXP   = 1.4        # temperature exponent   (Anderson 1970, Eq. 4)
_PI_MC_COEFF   = -0.07      # moisture coefficient   (Anderson 1970, Eq. 4)


# ---------------------------------------------------------------------------
# Core P_i calculation
# ---------------------------------------------------------------------------

def prob_ignition(T_fuel_F: float, MC_pct: float) -> float:
    """Compute Anderson (1970) probability of ignition.

    Uses the empirical formula reported in Rothermel (1983):

        P_i = min(1, max(0, _PI_COEFF × T_fuel_F^_PI_TEMP_EXP × exp(_PI_MC_COEFF × MC%)))

    Parameters
    ----------
    T_fuel_F : float
        Fine-fuel temperature [°F].  Typically ambient + solar increment.
    MC_pct : float
        1-hr dead fuel moisture content [% of oven-dry weight].

    Returns
    -------
    float
        Probability of ignition in [0, 1].  Returns 0 when T_fuel_F ≤ 0.
    """
    if T_fuel_F <= 0.0:
        return 0.0
    raw = _PI_COEFF * (T_fuel_F ** _PI_TEMP_EXP) * math.exp(_PI_MC_COEFF * MC_pct)
    return min(1.0, max(0.0, raw))


# ---------------------------------------------------------------------------
# Table computation
# ---------------------------------------------------------------------------

def compute_pi_table(
    temps_C: List[float],
    moistures_pct: List[float],
    solar_heating_F: float = 25.0,
) -> List[Dict]:
    """Compute probability-of-ignition matrix.

    Sweeps ambient temperature (°C) and 1-hr dead fuel moisture (%).
    The fine-fuel temperature is computed as:

        T_fuel_F = T_ambient_C × 9/5 + 32 + solar_heating_F

    Parameters
    ----------
    temps_C : list[float]
        Ambient air temperatures [°C].
    moistures_pct : list[float]
        1-hr dead fuel moisture values [%].
    solar_heating_F : float
        Solar temperature increment added to the ambient to get fuel
        temperature [°F].  Rothermel (1983) recommends 25 °F as a
        typical clear-sky default.

    Returns
    -------
    list[dict]
        One dict per (temperature, moisture) pair with keys:
        ``T_ambient_C``, ``T_fuel_F``, ``MC_pct``, ``P_ignition``.
    """
    rows = []
    for T_C in temps_C:
        T_a_F    = T_C * 9.0 / 5.0 + 32.0
        T_fuel_F = T_a_F + solar_heating_F
        for MC in moistures_pct:
            P_i = prob_ignition(T_fuel_F, MC)
            rows.append({
                "T_ambient_C": round(T_C, 1),
                "T_fuel_F":    round(T_fuel_F, 1),
                "MC_pct":      round(MC, 1),
                "P_ignition":  round(P_i, 4),
            })
    return rows


# ---------------------------------------------------------------------------
# ASCII table
# ---------------------------------------------------------------------------

def print_pi_table(rows: List[Dict]) -> None:
    """Print probability of ignition as a pivot table.

    Rows are ambient temperatures; columns are moisture contents.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_pi_table`.
    """
    if not rows:
        return

    temps  = sorted({r["T_ambient_C"] for r in rows})
    moist  = sorted({r["MC_pct"] for r in rows})
    lookup = {(r["T_ambient_C"], r["MC_pct"]): r["P_ignition"] for r in rows}

    col_w  = 8
    header = f"{'T_amb [°C]':>12s}" + "".join(f"  MC={m:4.1f}%" for m in moist)
    print("\nProbability of Ignition  P_i  (Anderson 1970 / Rothermel 1983)")
    print("-" * len(header))
    print(header)
    print("-" * len(header))
    for T in temps:
        row_str = f"{T:12.1f}"
        for m in moist:
            val = lookup.get((T, m), float("nan"))
            row_str += f"  {val:>8.3f}"
        print(row_str)
    print()


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def write_pi_csv(rows: List[Dict], out_path: str) -> None:
    """Write the P_i table rows to a CSV file.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_pi_table`.
    out_path : str
        Destination file path.
    """
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows → {out_path}")


# ---------------------------------------------------------------------------
# Optional heatmap
# ---------------------------------------------------------------------------

def make_pi_heatmap(rows: List[Dict], out_path: str) -> None:
    """Save a heatmap of P_i as a PNG figure.

    Requires ``matplotlib`` and ``numpy``.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_pi_table`.
    out_path : str
        Destination PNG file path.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib/numpy not available – skipping heatmap.",
              file=sys.stderr)
        return

    temps = sorted({r["T_ambient_C"] for r in rows})
    moist = sorted({r["MC_pct"] for r in rows})
    lookup = {(r["T_ambient_C"], r["MC_pct"]): r["P_ignition"] for r in rows}

    data = np.array([[lookup.get((T, m), 0.0) for m in moist] for T in temps])

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(data, aspect="auto", origin="lower", cmap="YlOrRd",
                   vmin=0.0, vmax=1.0,
                   extent=[min(moist) - 0.5, max(moist) + 0.5,
                           min(temps) - 1.0,  max(temps) + 1.0])
    ax.set_xlabel("1-hr Dead Fuel Moisture [%]")
    ax.set_ylabel("Ambient Temperature [°C]")
    ax.set_title("Probability of Ignition (Anderson 1970 / Rothermel 1983)")
    fig.colorbar(im, ax=ax, label="P_i")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved heatmap → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Anderson (1970) / Rothermel (1983) probability-of-ignition worksheet.\n\n"
            "Sweeps ambient temperature [°C] and 1-hr dead fuel moisture [%] to\n"
            "produce a P_i table indicating the likelihood that a firebrand landing\n"
            "on fine fuel will cause a sustained ignition."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--temp-min",    type=float, default=20.0,
                        help="Minimum ambient temperature [°C] (default: 20)")
    parser.add_argument("--temp-max",    type=float, default=45.0,
                        help="Maximum ambient temperature [°C] (default: 45)")
    parser.add_argument("--temp-steps",  type=int,   default=6,
                        help="Number of temperature steps (default: 6)")
    parser.add_argument("--moisture-min", type=float, default=2.0,
                        help="Minimum 1-hr dead fuel moisture [%%] (default: 2)")
    parser.add_argument("--moisture-max", type=float, default=20.0,
                        help="Maximum 1-hr dead fuel moisture [%%] (default: 20)")
    parser.add_argument("--moisture-steps", type=int, default=9,
                        help="Number of moisture steps (default: 9)")
    parser.add_argument("--solar-heating", type=float, default=25.0,
                        help="Solar temperature increment added to ambient [°F] "
                             "(default: 25 – Rothermel 1983 recommendation)")
    parser.add_argument("--out",         default="ignition_probability.csv",
                        help="Output CSV path (default: ignition_probability.csv)")
    parser.add_argument("--plot",        action="store_true",
                        help="Save heatmap PNG (requires matplotlib + numpy)")
    parser.add_argument("--plot-out",    default="ignition_probability_heatmap.png",
                        help="Heatmap PNG path (default: ignition_probability_heatmap.png)")

    args = parser.parse_args(argv)

    def linspace(lo, hi, n):
        if n == 1:
            return [lo]
        return [lo + (hi - lo) * i / (n - 1) for i in range(n)]

    temps  = linspace(args.temp_min,     args.temp_max,     args.temp_steps)
    moist  = linspace(args.moisture_min, args.moisture_max, args.moisture_steps)

    print(f"Temperature:   {len(temps)} steps  [{args.temp_min:.1f} – {args.temp_max:.1f} °C]")
    print(f"Moisture:      {len(moist)} steps  [{args.moisture_min:.1f} – {args.moisture_max:.1f} %]")
    print(f"Solar heating: {args.solar_heating:.1f} °F increment")

    rows = compute_pi_table(temps, moist, args.solar_heating)

    print_pi_table(rows)
    write_pi_csv(rows, args.out)

    if args.plot:
        make_pi_heatmap(rows, args.plot_out)


if __name__ == "__main__":
    main()
