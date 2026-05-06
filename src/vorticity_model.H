// ============================================================================
// Vorticity Diagnostic Field
// ============================================================================
// Computes the vertical component of vorticity (ω_z) from the 2-D (or 3-D
// projected) velocity field.  This is a standard diagnostic used in coupled
// fire-atmosphere codes (WRF-SFIRE, FireTEC, QUIC-Fire) to identify fire
// whirls, horizontal roll vortices, and regions of strong wind shear near
// the fire perimeter.
//
// Definition (2-D):
//   ω_z = ∂v/∂x − ∂u/∂y   [s⁻¹]
//
// where u = Ux (east-west component) and v = Uy (north-south component).
// Positive ω_z corresponds to counter-clockwise rotation (cyclonic in the
// northern hemisphere).
//
// Numerical approximation:
//   Centred finite differences on the cell-centred velocity MultiFab.
//   The velocity MultiFab must have at least 1 ghost cell filled (the
//   solver's "vel" field is declared with ng = 1 throughout).
//   A FillBoundary call is performed internally to ensure ghost cells are
//   consistent across MPI subdomain boundaries before the kernel is launched.
//
// References:
//   Weise, D.R. & Biging, G.S. (1996). Effects of wind velocity and
//     slope on flame properties. Canadian Journal of Forest Research.
//   Cunningham, P. & Linn, R.R. (2007). Numerical simulations of grass
//     fires using a coupled atmosphere-fire model. J. Geophys. Res.
// ============================================================================

#ifndef VORTICITY_MODEL_H
#define VORTICITY_MODEL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// compute_vorticity
//
// Parameters:
//   vorticity_mf – output: 1-component MultiFab, ω_z [s⁻¹]
//   vel          – input:  velocity MultiFab (components 0=u, 1=v; must have
//                          ≥ 1 ghost cell)
//   geom         – grid geometry (provides cell sizes)
// ============================================================================
static void compute_vorticity(MultiFab&       vorticity_mf,
                               MultiFab&       vel,
                               const Geometry& geom)
{
    // Fill ghost cells so centred differences work at MPI subdomain boundaries
    vel.FillBoundary(geom.periodicity());

    const auto dx = geom.CellSizeArray();
    const Real inv2dx = Real(1.0) / (Real(2.0) * dx[0]);
    const Real inv2dy = Real(1.0) / (Real(2.0) * dx[1]);

    for (MFIter mfi(vorticity_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v   = vel.const_array(mfi);
        auto       omz = vorticity_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // ∂v/∂x: centred difference of north-south component (comp 1)
            Real dvdx = (v(i + 1, j, k, 1) - v(i - 1, j, k, 1)) * inv2dx;
            // ∂u/∂y: centred difference of east-west component (comp 0)
            Real dudy = (v(i, j + 1, k, 0) - v(i, j - 1, k, 0)) * inv2dy;
            omz(i, j, k) = dvdx - dudy;
        });
    }
}

#endif // VORTICITY_MODEL_H
