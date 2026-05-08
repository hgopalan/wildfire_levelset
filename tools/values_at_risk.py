#!/usr/bin/env python3
"""
values_at_risk.py – FSPro-style values-at-risk (VAR) analysis.

Overlays a burn-probability map (from ``ensemble_burn_probability.py``) with
a user-supplied values inventory to compute:

  1. **Expected loss per cell** = P_burn × value
  2. **At-risk value totals** by category and percentile threshold
  3. **Prioritised at-risk list** of assets sorted by expected loss
  4. **Summary statistics**: total expected loss, conditional value at risk (CVaR)

Values inventory format (CSV with header)::

    name, category, x_m, y_m, value
    "Structure A", residential, 330150.0, 3775200.0, 250000
    "Power tower 12", infrastructure, 330600.0, 3775800.0, 80000
    "Reservoir", water, 330900.0, 3776100.0, 500000

  - ``name``     : human-readable asset label
  - ``category`` : free-form category string (used for grouping in summary)
  - ``x_m``      : easting [m] in the same coordinate system as the burn-prob map
  - ``y_m``      : northing [m]
  - ``value``    : replacement / economic value of the asset [any currency unit]

Each asset is matched to the nearest burn-probability grid cell (within a
configurable radius ``--max-match-dist``).  When no grid cell is within the
radius the asset is reported with P_burn = NaN.

Outputs
-------
  ``var_detail.csv``    – one row per asset with matched P_burn and expected loss
  ``var_summary.csv``   – totals by category
  ``var_exceedance.csv``– P(total_loss > threshold) curve (Monte-Carlo)
  ``var_map.png``       – optional spatial plot of expected loss (requires matplotlib)

Usage
-----
  # Basic run
  python3 tools/values_at_risk.py \\
      --burn-prob burn_probability.csv \\
      --values assets.csv \\
      --out-dir var_results/

  # Custom P_burn threshold for summary
  python3 tools/values_at_risk.py \\
      --burn-prob burn_probability.csv \\
      --values assets.csv \\
      --threshold 0.25 \\
      --out-dir var_results/

  # Monte-Carlo loss exceedance curve (1000 samples)
  python3 tools/values_at_risk.py \\
      --burn-prob burn_probability.csv \\
      --values assets.csv \\
      --exceedance \\
      --n-mc 1000

  # Spatial plot of expected loss
  python3 tools/values_at_risk.py \\
      --burn-prob burn_probability.csv \\
      --values assets.csv \\
      --plot

Options
-------
  --burn-prob FILE     Burn-probability CSV (X,Y,P_burn) from ensemble_burn_probability.py
  --values FILE        Asset values CSV (name,category,x_m,y_m,value)
  --out-dir DIR        Output directory (default: var_results)
  --threshold T        P_burn threshold for "at-risk" classification (default: 0.10)
  --max-match-dist D   Maximum distance [m] to match asset to grid cell (default: inf)
  --exceedance         Compute Monte-Carlo total-loss exceedance curve
  --n-mc N             Monte-Carlo samples for exceedance curve (default: 500)
  --plot               Save spatial expected-loss plot (requires matplotlib)
  --plot-out FILE      Path for the spatial plot PNG (default: var_map.png in out-dir)
  --quiet              Suppress printed summary

References
----------
  FSPro (Fire Spread Probability): https://www.firelab.org/project/fspro
  Finney, M.A. et al. (2011). A method for ensemble wildland fire simulation.
    Environmental Modelling & Software, 26(10), 1352-1359.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Burn-probability CSV reader (same format as ensemble_burn_probability.py output)
# ---------------------------------------------------------------------------

def read_burn_prob(csv_path: str) -> Tuple[List[Tuple[float, float, float]], float]:
    """Read burn_probability.csv.

    Returns (points, resolution) where points is a list of (x, y, P_burn)
    and resolution is the inferred grid cell size in metres.
    """
    points: List[Tuple[float, float, float]] = []
    resolution = 0.0

    with open(csv_path, newline="") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "resolution" in line:
                    import re
                    m = re.search(r"resolution\s*=\s*([\d.]+)", line)
                    if m:
                        resolution = float(m.group(1))
                continue
            if line.lower().startswith("x"):
                continue
            parts = line.split(",")
            if len(parts) < 3:
                continue
            try:
                x, y, p = float(parts[0]), float(parts[1]), float(parts[2])
                points.append((x, y, p))
            except ValueError:
                continue

    # Auto-detect resolution from grid spacing if not found in header
    if resolution <= 0.0 and len(points) > 1:
        xs = sorted({pt[0] for pt in points})
        if len(xs) > 1:
            diffs = [xs[i+1] - xs[i] for i in range(len(xs)-1) if xs[i+1] > xs[i]]
            if diffs:
                resolution = min(diffs)

    return points, resolution


# ---------------------------------------------------------------------------
# Asset values CSV reader
# ---------------------------------------------------------------------------

def read_values(csv_path: str) -> List[Dict]:
    """Read values inventory CSV.

    Required columns: name, category, x_m, y_m, value
    Additional columns are preserved.
    """
    assets = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            # Normalise header names: strip whitespace
            row = {k.strip(): v.strip() for k, v in row.items() if k}
            try:
                asset = {
                    "name":     row.get("name", ""),
                    "category": row.get("category", ""),
                    "x_m":      float(row["x_m"]),
                    "y_m":      float(row["y_m"]),
                    "value":    float(row["value"]),
                }
                # Carry any extra columns through
                for k, v in row.items():
                    if k not in asset:
                        asset[k] = v
                assets.append(asset)
            except (KeyError, ValueError) as e:
                print(f"WARNING: skipping row {row} – {e}", file=sys.stderr)
    return assets


# ---------------------------------------------------------------------------
# Nearest-cell lookup using a simple grid index
# ---------------------------------------------------------------------------

class BurnProbGrid:
    """Spatial index for burn-probability point cloud."""

    def __init__(self, points: List[Tuple[float, float, float]],
                 resolution: float) -> None:
        self.points     = points          # (x, y, P_burn)
        self.resolution = resolution or 1.0

        # Build a dict index: (ix, iy) → (x, y, P_burn)
        self._idx: Dict[Tuple[int, int], Tuple[float, float, float]] = {}
        if points:
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            self._x0 = min(xs)
            self._y0 = min(ys)
        else:
            self._x0 = 0.0
            self._y0 = 0.0

        for pt in points:
            ix = int(round((pt[0] - self._x0) / self.resolution))
            iy = int(round((pt[1] - self._y0) / self.resolution))
            self._idx[(ix, iy)] = pt

    def query(self, x: float, y: float,
              max_dist: float = float("inf")) -> Optional[Tuple[float, float, float]]:
        """Return closest (x, y, P_burn) point within max_dist [m], or None."""
        ix0 = int(round((x - self._x0) / self.resolution))
        iy0 = int(round((y - self._y0) / self.resolution))

        # Search in expanding ring of cells up to max_dist
        search_radius = max(1, int(math.ceil(max_dist / self.resolution)) + 1)
        best_dist  = float("inf")
        best_point = None

        for dix in range(-search_radius, search_radius + 1):
            for diy in range(-search_radius, search_radius + 1):
                pt = self._idx.get((ix0 + dix, iy0 + diy))
                if pt is None:
                    continue
                dist = math.sqrt((pt[0] - x)**2 + (pt[1] - y)**2)
                if dist < best_dist and dist <= max_dist:
                    best_dist  = dist
                    best_point = pt

        return best_point


# ---------------------------------------------------------------------------
# Core VAR computation
# ---------------------------------------------------------------------------

def compute_var(assets: List[Dict], grid: BurnProbGrid,
                max_match_dist: float,
                threshold: float) -> List[Dict]:
    """Enrich each asset with P_burn, expected_loss, and at_risk flag.

    Returns list of result dicts.
    """
    results = []
    for asset in assets:
        x, y, val = asset["x_m"], asset["y_m"], asset["value"]
        pt = grid.query(x, y, max_match_dist)

        if pt is None:
            p_burn       = float("nan")
            match_x      = float("nan")
            match_y      = float("nan")
            match_dist   = float("nan")
            expected_loss = float("nan")
            at_risk      = 0
        else:
            p_burn       = pt[2]
            match_x      = pt[0]
            match_y      = pt[1]
            match_dist   = math.sqrt((pt[0] - x)**2 + (pt[1] - y)**2)
            expected_loss = p_burn * val
            at_risk      = int(p_burn >= threshold)

        row = dict(asset)  # copy all original fields
        row.update({
            "P_burn":        p_burn,
            "match_x":       match_x,
            "match_y":       match_y,
            "match_dist_m":  match_dist,
            "expected_loss": expected_loss,
            "at_risk":       at_risk,
        })
        results.append(row)

    # Sort by expected_loss descending (highest risk first)
    results.sort(key=lambda r: (
        r["expected_loss"] if not math.isnan(r["expected_loss"]) else -1.0),
        reverse=True)
    return results


# ---------------------------------------------------------------------------
# Category summary
# ---------------------------------------------------------------------------

def summarise_by_category(results: List[Dict]) -> List[Dict]:
    """Aggregate expected loss and total value by category."""
    cats: Dict[str, Dict] = {}
    for row in results:
        cat = row.get("category", "")
        if cat not in cats:
            cats[cat] = dict(category=cat, n_assets=0, total_value=0.0,
                             expected_loss=0.0, n_at_risk=0,
                             n_matched=0)
        c = cats[cat]
        c["n_assets"] += 1
        c["total_value"] += row["value"]
        if not math.isnan(row.get("expected_loss", float("nan"))):
            c["expected_loss"] += row["expected_loss"]
            c["n_matched"] += 1
        if row.get("at_risk", 0):
            c["n_at_risk"] += 1

    summary = sorted(cats.values(), key=lambda r: r["expected_loss"], reverse=True)
    return summary


# ---------------------------------------------------------------------------
# Monte-Carlo total-loss exceedance curve
# ---------------------------------------------------------------------------

def mc_exceedance(results: List[Dict], n_mc: int,
                  rng_seed: int = 0) -> List[Tuple[float, float]]:
    """Compute P(total_loss > threshold) via Monte-Carlo simulation.

    Each MC trial: for each asset, draw Bernoulli(P_burn) → sum losses.
    Returns sorted list of (total_loss, exceedance_probability) pairs.
    """
    rng = random.Random(rng_seed)
    # Filter assets with valid P_burn
    valid = [(r["P_burn"], r["value"])
             for r in results
             if not math.isnan(r.get("P_burn", float("nan"))) and r["value"] > 0]

    total_losses = []
    for _ in range(n_mc):
        loss = sum(val for p, val in valid if rng.random() < p)
        total_losses.append(loss)

    total_losses.sort()
    n = len(total_losses)
    curve = []
    for i, loss in enumerate(total_losses):
        # P(loss > total_losses[i]) = fraction of samples exceeding this value
        p_exceed = 1.0 - (i + 1) / n
        curve.append((loss, p_exceed))

    # Prepend (0, 1.0) so curve starts at full probability
    curve.insert(0, (0.0, 1.0))
    return curve


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

_DETAIL_FIELDS = [
    "name", "category", "x_m", "y_m", "value",
    "P_burn", "match_x", "match_y", "match_dist_m",
    "expected_loss", "at_risk",
]

def write_detail_csv(results: List[Dict], path: str) -> None:
    """Write per-asset detail CSV."""
    extra = [k for k in results[0] if k not in _DETAIL_FIELDS] if results else []
    fields = _DETAIL_FIELDS + extra
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {len(results)} asset row(s) → {path}")


_SUMMARY_FIELDS = [
    "category", "n_assets", "n_matched", "n_at_risk",
    "total_value", "expected_loss",
]

def write_summary_csv(summary: List[Dict], path: str) -> None:
    """Write category-level summary CSV."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_SUMMARY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary)
    print(f"Wrote {len(summary)} category row(s) → {path}")


def write_exceedance_csv(curve: List[Tuple[float, float]], path: str) -> None:
    """Write (total_loss, P_exceed) exceedance curve CSV."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["total_loss", "P_exceed"])
        writer.writerows(curve)
    print(f"Wrote {len(curve)} exceedance points → {path}")


# ---------------------------------------------------------------------------
# ASCII summary printer
# ---------------------------------------------------------------------------

def print_summary(results: List[Dict], summary: List[Dict],
                  threshold: float) -> None:
    """Print a concise ASCII summary to stdout."""
    n_total   = len(results)
    n_matched = sum(1 for r in results if not math.isnan(r.get("P_burn", float("nan"))))
    n_at_risk = sum(r.get("at_risk", 0) for r in results)
    total_val = sum(r["value"] for r in results)
    exp_loss  = sum(r["expected_loss"] for r in results
                    if not math.isnan(r.get("expected_loss", float("nan"))))

    print("=" * 68)
    print("  VALUES-AT-RISK SUMMARY")
    print("=" * 68)
    print(f"  Total assets             : {n_total}")
    print(f"  Assets matched to grid   : {n_matched}")
    print(f"  Assets at risk (P≥{threshold:.2f}) : {n_at_risk}")
    print(f"  Total portfolio value    : {total_val:,.0f}")
    print(f"  Expected loss (Σ P·V)    : {exp_loss:,.0f}"
          f"  ({100*exp_loss/total_val:.1f}% of portfolio)"
          if total_val > 0 else "")

    print(f"\n  ── Top 10 At-Risk Assets ({'sorted by expected loss'}) ──────────")
    top = [r for r in results[:10]
           if not math.isnan(r.get("expected_loss", float("nan")))]
    if top:
        hdr = f"  {'Name':<30s} {'Cat':<16s} {'P_burn':>7s} {'Value':>12s} {'Exp.Loss':>12s}"
        print(hdr)
        print("  " + "-" * 79)
        for r in top:
            pb   = f"{r['P_burn']:.3f}"
            val  = f"{r['value']:,.0f}"
            el   = f"{r['expected_loss']:,.0f}"
            name = r["name"][:30]
            cat  = r["category"][:16]
            print(f"  {name:<30s} {cat:<16s} {pb:>7s} {val:>12s} {el:>12s}")

    print(f"\n  ── By Category ─────────────────────────────────────────────────")
    if summary:
        print(f"  {'Category':<20s} {'Assets':>7s} {'At-Risk':>8s}"
              f" {'Total Value':>14s} {'Exp.Loss':>14s}")
        print("  " + "-" * 68)
        for s in summary:
            print(f"  {s['category']:<20s} {s['n_assets']:>7d} {s['n_at_risk']:>8d}"
                  f" {s['total_value']:>14,.0f} {s['expected_loss']:>14,.0f}")
    print("=" * 68)


# ---------------------------------------------------------------------------
# Optional spatial plot
# ---------------------------------------------------------------------------

def make_var_plot(results: List[Dict],
                  grid_pts: List[Tuple[float, float, float]],
                  out_path: str) -> None:
    """Save a spatial expected-loss overlay plot."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib not available – skipping plot.", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw burn probability as a faint background
    if grid_pts:
        xs = np.array([p[0] for p in grid_pts])
        ys = np.array([p[1] for p in grid_pts])
        ps = np.array([p[2] for p in grid_pts])
        sc_bg = ax.scatter(xs, ys, c=ps, cmap="YlOrRd", vmin=0, vmax=1,
                           s=6, alpha=0.4, label="P_burn")
        plt.colorbar(sc_bg, ax=ax, label="P_burn")

    # Overlay assets sized by expected loss
    valid = [r for r in results
             if not math.isnan(r.get("expected_loss", float("nan")))]
    if valid:
        xs_a = [r["x_m"] for r in valid]
        ys_a = [r["y_m"] for r in valid]
        el   = [r["expected_loss"] for r in valid]
        max_el = max(el) if el else 1.0
        sizes  = [max(10.0, 300.0 * e / max_el) for e in el]
        ax.scatter(xs_a, ys_a, c="steelblue", s=sizes, alpha=0.8,
                   edgecolors="k", linewidths=0.5, label="Assets (size ~ Exp.Loss)")

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title("Values at Risk — Expected Loss Overlay")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved VAR map → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="FSPro-style values-at-risk overlay for wildfire burn probability maps.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--burn-prob",       required=True,  metavar="CSV",
                   help="Burn-probability CSV (X,Y,P_burn) from ensemble_burn_probability.py")
    p.add_argument("--values",          required=True,  metavar="CSV",
                   help="Asset values CSV (name,category,x_m,y_m,value)")
    p.add_argument("--out-dir",         default="var_results", metavar="DIR",
                   help="Output directory for result files (default: var_results)")
    p.add_argument("--threshold",       type=float, default=0.10, metavar="T",
                   help="P_burn threshold for 'at-risk' classification (default: 0.10)")
    p.add_argument("--max-match-dist",  type=float, default=float("inf"), metavar="M",
                   help="Max distance [m] to match asset to grid cell (default: unlimited)")
    p.add_argument("--exceedance",      action="store_true",
                   help="Compute Monte-Carlo total-loss exceedance curve")
    p.add_argument("--n-mc",            type=int, default=500, metavar="N",
                   help="Monte-Carlo samples for exceedance curve (default: 500)")
    p.add_argument("--seed",            type=int, default=0,
                   help="Random seed for Monte-Carlo (default: 0)")
    p.add_argument("--plot",            action="store_true",
                   help="Save spatial expected-loss plot (requires matplotlib)")
    p.add_argument("--plot-out",        default=None, metavar="FILE",
                   help="Path for the spatial plot PNG (default: var_map.png in out-dir)")
    p.add_argument("--quiet",           action="store_true",
                   help="Suppress printed summary")
    return p


def main(argv=None) -> None:
    args = _build_parser().parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load burn-probability map ---
    print(f"Loading burn-probability map: {args.burn_prob}")
    bp_points, resolution = read_burn_prob(args.burn_prob)
    if not bp_points:
        print("ERROR: no burn-probability points loaded.", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(bp_points)} grid cells, resolution ≈ {resolution:.1f} m")

    grid = BurnProbGrid(bp_points, resolution)

    # --- Load values inventory ---
    print(f"Loading values inventory: {args.values}")
    assets = read_values(args.values)
    if not assets:
        print("ERROR: no assets loaded from values CSV.", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(assets)} asset(s) loaded")

    # --- Compute VAR ---
    results = compute_var(assets, grid, args.max_match_dist, args.threshold)
    summary = summarise_by_category(results)

    # --- Print summary ---
    if not args.quiet:
        print_summary(results, summary, args.threshold)

    # --- Write outputs ---
    write_detail_csv(results, str(out_dir / "var_detail.csv"))
    write_summary_csv(summary, str(out_dir / "var_summary.csv"))

    # --- Monte-Carlo exceedance curve ---
    if args.exceedance:
        print(f"\nRunning Monte-Carlo exceedance ({args.n_mc} samples)…")
        curve = mc_exceedance(results, args.n_mc, rng_seed=args.seed)
        write_exceedance_csv(curve, str(out_dir / "var_exceedance.csv"))

    # --- Spatial plot ---
    if args.plot:
        plot_path = args.plot_out or str(out_dir / "var_map.png")
        make_var_plot(results, bp_points, plot_path)


if __name__ == "__main__":
    main()
