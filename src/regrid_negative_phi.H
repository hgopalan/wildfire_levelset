// ================================================================
// regrid_negative_phi.H
// Regrids the domain every N steps, keeping only boxes that
// contain at least one cell with phi < 0 (the interior region).
// Both phi and vel are re-created on the new BoxArray so they
// remain consistent.
// ================================================================
#ifndef REGRID_NEGATIVE_PHI_H
#define REGRID_NEGATIVE_PHI_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>

using namespace amrex;

// ----------------------------------------------------------------
// regrid_negative_phi
//   phi     – level-set field  (ng_phi ghost cells)
//   vel     – velocity field   (1 ghost cell)
//   geom    – geometry object
//   ng_phi  – ghost-cell count for phi (typically 3)
// ----------------------------------------------------------------
static void regrid_negative_phi(MultiFab& phi,
                                MultiFab& vel,
                                const Geometry& geom,
                                int ng_phi = 3)
{
    const BoxArray& old_ba = phi.boxArray();
    const int nboxes = old_ba.size();

    // ----- 1. Tag boxes that contain at least one phi < 0 cell -----
    Vector<Box> kept_boxes;
    kept_boxes.reserve(nboxes);

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& vbx = mfi.validbox();
        Real box_min = phi[mfi].min<RunOn::Host>(vbx, 0);
        if (box_min < 0.0) {
            kept_boxes.push_back(vbx);
        }
    }

    // ----- 2. Guard: if nothing is negative, skip -----
    if (kept_boxes.empty()) {
        amrex::Print() << "  regrid_negative_phi: no phi<0 boxes found, skipping.\n";
        return;
    }

    amrex::Print() << "  regrid_negative_phi: keeping "
                   << kept_boxes.size() << " / " << nboxes << " boxes.\n";

    // ----- 3. Build new BoxArray & DistributionMapping -----
    BoxArray new_ba(BoxList(std::move(kept_boxes)));
    DistributionMapping new_dm(new_ba);

    // ----- 4. Rebuild phi on new layout -----
    MultiFab new_phi(new_ba, new_dm, 1, ng_phi);
    new_phi.setVal(1.0e30);                        // sentinel for unfilled cells
    new_phi.ParallelCopy(phi, 0, 0, 1, 0, 0);      // copy valid data
    fill_boundary_extrap(new_phi, geom);            // refresh ghost cells

    // ----- 5. Rebuild vel on new layout -----
    MultiFab new_vel(new_ba, new_dm, 3, 1);
    new_vel.setVal(0.0);
    new_vel.ParallelCopy(vel, 0, 0, 3, 0, 0);
    new_vel.FillBoundary(geom.periodicity());

    // ----- 6. Swap into caller's MultiFabs -----
    std::swap(phi, new_phi);
    std::swap(vel, new_vel);
}

#endif // REGRID_NEGATIVE_PHI_H
