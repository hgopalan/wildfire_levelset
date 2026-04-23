#ifndef REGRID_NEGATIVE_PHI_H
#define REGRID_NEGATIVE_PHI_H

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

static bool regrid_negative_phi(MultiFab& phi,
                                MultiFab& vel,
                                BoxArray& ba,
                                DistributionMapping& dm,
                                Geometry& geom,
                                int ng_phi,
                                int refine_ratio,
                                amrex::Real tag_phi_threshold,
                                int max_refinements,
                                int& refinement_count)
{
    if (refinement_count >= max_refinements) { return false; }
    if (refine_ratio <= 1) { return false; }

    const amrex::Real phi_min = phi.min(0);
    if (phi_min >= tag_phi_threshold) { return false; }

    const Box old_domain = geom.Domain();
    const auto old_prob_lo = geom.ProbLoArray();
    const auto old_prob_hi = geom.ProbHiArray();
    const auto old_periodic = geom.isPeriodicArray();
    const int coord = geom.Coord();

    MultiFab old_phi(ba, dm, 1, ng_phi);
    MultiFab old_vel(ba, dm, 3, 1);
    MultiFab::Copy(old_phi, phi, 0, 0, 1, ng_phi);
    MultiFab::Copy(old_vel, vel, 0, 0, 3, 1);

    ba.refine(refine_ratio);
    // Keep processor ownership aligned with the pre-refinement layout so that
    // coarse-to-fine transfer can use matching local MFIter indices.
    dm = DistributionMapping(old_phi.DistributionMap().ProcessorMap());

    const Box new_domain = amrex::refine(old_domain, IntVect(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio)));
    RealBox rb({AMREX_D_DECL(old_prob_lo[0], old_prob_lo[1], old_prob_lo[2])},
               {AMREX_D_DECL(old_prob_hi[0], old_prob_hi[1], old_prob_hi[2])});
    Geometry new_geom(new_domain, &rb, coord, old_periodic.data());

    phi.define(ba, dm, 1, ng_phi);
    vel.define(ba, dm, 3, 1);

    const int rr = refine_ratio;
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p_new = phi.array(mfi);
        auto const v_new = vel.array(mfi);
        auto const p_old = old_phi.const_array(mfi);
        auto const v_old = old_vel.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int ic = i / rr;
            const int jc = j / rr;
            const int kc = k / rr;

            p_new(i,j,k) = p_old(ic,jc,kc);
            v_new(i,j,k,0) = v_old(ic,jc,kc,0);
            v_new(i,j,k,1) = v_old(ic,jc,kc,1);
            v_new(i,j,k,2) = v_old(ic,jc,kc,2);
        });
    }

    geom = new_geom;
    ++refinement_count;
    amrex::Print() << "Refined grid on negative phi tags (phi < " << tag_phi_threshold
                   << "), refinement_count=" << refinement_count << "\n";
    return true;
}

#endif
