// spatial_grid.H – Cell-list spatial acceleration for scattered-data interpolation.
//
// Provides a CPU-side builder and a GPU-callable query function that together
// replace the O(N) full-scan Gaussian RBF loop with an O(K) local search,
// where K is the (small) number of data points within the support radius.
//
// Kernel choice – Wendland C² (compactly supported):
//
//   phi(r) = (1 - r/R)^4 * (4*r/R + 1),  r <= R
//   phi(r) = 0,                            r > R
//
// with R = 3*h (three times the characteristic data spacing h = sqrt(area/N)).
//
// This is the "compromise between IDW and RBF":
//   • Smooth everywhere (C², no IDW bullseye artefact)
//   • Zero beyond R    (compact support → fast spatial lookup, like IDW)
//   • Cheap polynomial (no exp/transcendental calls)
//
// Usage:
//   1. Call build_cell_list() on the CPU to sort data into uniform bins and
//      copy them to device memory.
//   2. Set the sx/sy/bstart/bcount pointers in the returned CellListParams from
//      the corresponding Gpu::DeviceVector::data() calls.
//   3. In a GPU kernel (AMREX_GPU_DEVICE lambda), call cell_list_interp() to
//      evaluate the weighted average at a query point.

#ifndef SPATIAL_GRID_H
#define SPATIAL_GRID_H

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace amrex;

// ---------------------------------------------------------------------------
// GPU-visible descriptor (plain-old-data, safe to capture by value in lambdas)
// ---------------------------------------------------------------------------
struct CellListParams {
    Real const* sx;      // sorted x coordinates (device pointer)
    Real const* sy;      // sorted y coordinates (device pointer)
    int  const* bstart;  // start index in sorted arrays for each bin
    int  const* bcount;  // number of points in each bin
    Real xmin, ymin;     // domain origin (lower-left corner of bin grid)
    Real cell_sz;        // bin edge length (== cutoff radius R)
    int  nbx, nby;       // number of bins along x and y
    int  n_pts;          // total number of data points
    Real R;              // cutoff radius (== cell_sz)
};

// ---------------------------------------------------------------------------
// CPU helper – build the cell-list structure from host x/y vectors.
//
// Parameters
//   hx, hy      – host-side x and y coordinates of the data points
//   h_spacing   – characteristic data spacing (sqrt(area / N)); callers
//                 should pass the same value they use for their kernel bandwidth
//   d_sx/d_sy   – output: device vectors of sorted coordinates
//   d_bstart    – output: device vector of bin start indices (size nbx*nby)
//   d_bcount    – output: device vector of bin point counts   (size nbx*nby)
//   params      – output: filled with all scalar fields; CALLER must set the
//                 four pointer fields (sx, sy, bstart, bcount) after this call
//                 using the data() of the corresponding DeviceVectors.
//
// Returns
//   A host-side vector sorted_idx where sorted_idx[new_pos] == original_pos.
//   Use it to reorder any per-point value arrays to match the sorted layout.
// ---------------------------------------------------------------------------
static std::vector<int> build_cell_list(
        const std::vector<Real>& hx,
        const std::vector<Real>& hy,
        Real                     h_spacing,
        Gpu::DeviceVector<Real>& d_sx,
        Gpu::DeviceVector<Real>& d_sy,
        Gpu::DeviceVector<int>&  d_bstart,
        Gpu::DeviceVector<int>&  d_bcount,
        CellListParams&          params)
{
    const int n = static_cast<int>(hx.size());

    // Cutoff = 3 * characteristic spacing  (Wendland C2 support radius)
    const Real R        = Real(3.0) * h_spacing;
    const Real cell_sz  = std::max(R, Real(1.0e-6));  // guard against degenerate spacing

    const Real xmin = *std::min_element(hx.begin(), hx.end());
    const Real xmax = *std::max_element(hx.begin(), hx.end());
    const Real ymin = *std::min_element(hy.begin(), hy.end());
    const Real ymax = *std::max_element(hy.begin(), hy.end());

    const int nbx = static_cast<int>(std::ceil((xmax - xmin) / cell_sz)) + 1;
    const int nby = static_cast<int>(std::ceil((ymax - ymin) / cell_sz)) + 1;
    const int n_bins = nbx * nby;

    // --- assign each point to a bin ---
    std::vector<int> bin_of(n);
    std::vector<int> count(n_bins, 0);
    for (int i = 0; i < n; ++i) {
        int bx = static_cast<int>((hx[i] - xmin) / cell_sz);
        int by = static_cast<int>((hy[i] - ymin) / cell_sz);
        bx = std::max(0, std::min(nbx - 1, bx));
        by = std::max(0, std::min(nby - 1, by));
        bin_of[i] = by * nbx + bx;
        ++count[bin_of[i]];
    }

    // --- exclusive prefix-sum → bin start indices ---
    std::vector<int> h_bstart(n_bins, 0);
    for (int b = 1; b < n_bins; ++b)
        h_bstart[b] = h_bstart[b - 1] + count[b - 1];

    // --- scatter points into sorted arrays ---
    std::vector<Real> sx(n), sy(n);
    std::vector<int>  sorted_idx(n);
    std::vector<int>  fill(h_bstart);   // working copy of starts as fill cursors

    for (int i = 0; i < n; ++i) {
        const int b   = bin_of[i];
        const int pos = fill[b]++;
        sx[pos]          = hx[i];
        sy[pos]          = hy[i];
        sorted_idx[pos]  = i;
    }

    // --- copy to device ---
    d_sx.resize(n);     d_sy.resize(n);
    d_bstart.resize(n_bins);  d_bcount.resize(n_bins);

    Gpu::copy(Gpu::hostToDevice, sx.begin(),       sx.end(),       d_sx.begin());
    Gpu::copy(Gpu::hostToDevice, sy.begin(),       sy.end(),       d_sy.begin());
    Gpu::copy(Gpu::hostToDevice, h_bstart.begin(), h_bstart.end(), d_bstart.begin());
    Gpu::copy(Gpu::hostToDevice, count.begin(),    count.end(),    d_bcount.begin());

    // --- fill scalar params (caller sets pointer fields after this) ---
    params.xmin    = xmin;
    params.ymin    = ymin;
    params.cell_sz = cell_sz;
    params.nbx     = nbx;
    params.nby     = nby;
    params.n_pts   = n;
    params.R       = R;

    return sorted_idx;   // caller uses this to reorder value arrays
}

// ---------------------------------------------------------------------------
// Reorder a host value array to match the sorted layout produced by
// build_cell_list().  sorted_idx is the return value of build_cell_list().
// ---------------------------------------------------------------------------
static std::vector<Real> reorder_values(const std::vector<Real>& vals,
                                        const std::vector<int>&  sorted_idx)
{
    const int n = static_cast<int>(vals.size());
    std::vector<Real> out(n);
    for (int i = 0; i < n; ++i)
        out[i] = vals[sorted_idx[i]];
    return out;
}

// ---------------------------------------------------------------------------
// GPU kernel helper – Wendland C² weighted interpolation via cell list.
//
// Evaluates the normalised weighted average
//
//   f(p) = sum_i phi(|p - x_i| / R) * vals[i]
//          ----------------------------------------
//          sum_i phi(|p - x_i| / R)
//
// with phi(t) = (1-t)^4 (4t+1) for t in [0,1], 0 otherwise.
//
// Only the 3×3 = 9 bins centred on the query bin are searched.
//
// Falls back to nearest-neighbour value if the weighted sum is zero (i.e.
// the query point has no data within radius R – possible near boundaries).
// ---------------------------------------------------------------------------
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real cell_list_interp(const CellListParams& cl,
                      const Real* AMREX_RESTRICT vals,
                      Real px, Real py) noexcept
{
    // Bin containing the query point (clamped to valid range)
    int bx = static_cast<int>((px - cl.xmin) / cl.cell_sz);
    int by = static_cast<int>((py - cl.ymin) / cl.cell_sz);
    bx = max(0, min(cl.nbx - 1, bx));
    by = max(0, min(cl.nby - 1, by));

    Real R_inv = Real(1.0) / cl.R;
    Real R2    = cl.R * cl.R;
    Real sw = Real(0.0);
    Real sv = Real(0.0);

    // Nearest-neighbour fallback fields
    Real best_d2 = Real(1.0e30);
    Real best_v  = Real(0.0);

    // Search the 3×3 neighbourhood of bins (radius 1 bin = 1 cutoff)
    for (int dby = max(0, by - 1); dby <= min(cl.nby - 1, by + 1); ++dby) {
        for (int dbx = max(0, bx - 1); dbx <= min(cl.nbx - 1, bx + 1); ++dbx) {
            const int b   = dby * cl.nbx + dbx;
            const int end = cl.bstart[b] + cl.bcount[b];
            for (int k = cl.bstart[b]; k < end; ++k) {
                const Real ddx = px - cl.sx[k];
                const Real ddy = py - cl.sy[k];
                const Real d2  = ddx * ddx + ddy * ddy;

                // Wendland C2: phi(t) = (1-t)^4*(4t+1), t = sqrt(d2)/R
                if (d2 <= R2) {
                    const Real t   = std::sqrt(d2) * R_inv;
                    const Real om  = Real(1.0) - t;           // (1 - t)
                    const Real om2 = om * om;
                    const Real w   = om2 * om2 * (Real(4.0) * t + Real(1.0));
                    sw += w;
                    sv += w * vals[k];
                }

                // Track nearest point for fallback
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_v  = vals[k];
                }
            }
        }
    }

    if (sw > Real(0.0)) return sv / sw;

    // Fallback: if no point is within R (boundary / extrapolation), return the
    // nearest-neighbour value found in the searched bins.  If nothing was found
    // at all (empty neighbourhood), expand search one extra ring.
    if (best_d2 < Real(1.0e30)) return best_v;

    // Last resort: wider search (2-bin radius) to handle very sparse data
    for (int dby = max(0, by - 2); dby <= min(cl.nby - 1, by + 2); ++dby) {
        for (int dbx = max(0, bx - 2); dbx <= min(cl.nbx - 1, bx + 2); ++dbx) {
            const int b   = dby * cl.nbx + dbx;
            const int end = cl.bstart[b] + cl.bcount[b];
            for (int k = cl.bstart[b]; k < end; ++k) {
                const Real ddx = px - cl.sx[k];
                const Real ddy = py - cl.sy[k];
                const Real d2  = ddx * ddx + ddy * ddy;
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_v  = vals[k];
                }
            }
        }
    }
    return best_v;
}

#endif // SPATIAL_GRID_H
