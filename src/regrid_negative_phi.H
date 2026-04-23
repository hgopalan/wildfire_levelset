#ifndef REGRID_NEGATIVE_PHI_H
#define REGRID_NEGATIVE_PHI_H

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Print.H>
#include <algorithm>
#include <limits>
#include <memory>

static void initialize_fine_from_coarse(const MultiFab& coarse_phi,
                                        const MultiFab& coarse_vel,
                                        MultiFab& fine_phi,
                                        MultiFab& fine_vel,
                                        int refine_ratio)
{
    BoxArray coarse_from_fine_ba = fine_phi.boxArray();
    coarse_from_fine_ba.coarsen(refine_ratio);

    const DistributionMapping& fine_level_dm = fine_phi.DistributionMap();
    MultiFab coarse_phi_on_fine_dm(coarse_from_fine_ba, fine_level_dm, 1, 0);
    MultiFab coarse_vel_on_fine_dm(coarse_from_fine_ba, fine_level_dm, 3, 0);
    coarse_phi_on_fine_dm.ParallelCopy(coarse_phi, 0, 0, 1, 0, 0);
    coarse_vel_on_fine_dm.ParallelCopy(coarse_vel, 0, 0, 3, 0, 0);

    for (MFIter mfi(fine_phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p_fine = fine_phi.array(mfi);
        auto const v_fine = fine_vel.array(mfi);
        auto const p_coarse = coarse_phi_on_fine_dm.const_array(mfi);
        auto const v_coarse = coarse_vel_on_fine_dm.const_array(mfi);

        const int rr = refine_ratio;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int ic = i / rr;
            const int jc = j / rr;
            const int kc = k / rr;

            p_fine(i,j,k) = p_coarse(ic,jc,kc);
            v_fine(i,j,k,0) = v_coarse(ic,jc,kc,0);
            v_fine(i,j,k,1) = v_coarse(ic,jc,kc,1);
            v_fine(i,j,k,2) = v_coarse(ic,jc,kc,2);
        });
    }
}

static void synchronize_coarse_from_fine(MultiFab& coarse_phi,
                                         MultiFab& coarse_vel,
                                         const MultiFab& fine_phi,
                                         const MultiFab& fine_vel,
                                         int refine_ratio)
{
    const IntVect rr(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio));
    amrex::average_down(fine_phi, coarse_phi, 0, 1, rr);
    amrex::average_down(fine_vel, coarse_vel, 0, 3, rr);
}

static void synchronize_coarse_phi_from_fine(MultiFab& coarse_phi,
                                             const MultiFab& fine_phi,
                                             int refine_ratio)
{
    const IntVect rr(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio));
    amrex::average_down(fine_phi, coarse_phi, 0, 1, rr);
}

static int compute_fine_max_grid_size(int refine_ratio, int max_grid_size)
{
    const long long fine_max_grid_size_ll = static_cast<long long>(refine_ratio) *
                                            static_cast<long long>(max_grid_size);
    return static_cast<int>(std::max(1LL, std::min(fine_max_grid_size_ll,
                                                   static_cast<long long>(std::numeric_limits<int>::max()))));
}

static bool regrid_negative_phi(const MultiFab& coarse_phi,
                                const MultiFab& coarse_vel,
                                const Geometry& coarse_geom,
                                int ng_phi,
                                int refine_ratio,
                                amrex::Real tag_phi_threshold,
                                int max_refinements,
                                int max_grid_size,
                                bool& has_fine_level,
                                BoxArray& fine_ba,
                                DistributionMapping& fine_dm,
                                std::unique_ptr<Geometry>& fine_geom,
                                std::unique_ptr<MultiFab>& fine_phi,
                                std::unique_ptr<MultiFab>& fine_vel)
{
    if (refine_ratio <= 1 || max_refinements < 1) {
        if (has_fine_level) {
            has_fine_level = false;
            fine_geom.reset();
            fine_phi.reset();
            fine_vel.reset();
            amrex::Print() << "Removed fine level (local refinement disabled)\n";
            return true;
        }
        return false;
    }
    amrex::Print() << "Regridding for negative phi with tag threshold = " << tag_phi_threshold << "\n";
    BoxList tagged_cells;
    for (MFIter mfi(coarse_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (coarse_phi[mfi].min(bx, 0) < tag_phi_threshold) {
            tagged_cells.push_back(bx);
        }
    }
    amrex::Print() << "Number of tagged coarse boxes: " << tagged_cells.size() << "\n";
    if (tagged_cells.isEmpty()) {
        if (has_fine_level) {
            has_fine_level = false;
            fine_geom.reset();
            fine_phi.reset();
            fine_vel.reset();
            amrex::Print() << "Removed fine level (no phi < " << tag_phi_threshold << " tags)\n";
            return true;
        }
        return false;
    }
    amrex::Print() << "Refining " << tagged_cells.size() << " coarse boxes with refine ratio " << refine_ratio << "\n";
    BoxList fine_bl;
    for (const Box& bx : tagged_cells) {
        fine_bl.push_back(amrex::refine(bx, IntVect(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio))));
    }
    amrex::Print() << "Number of fine boxes before merging: " << fine_bl.size() << "\n";
    BoxArray new_fine_ba(fine_bl);
    new_fine_ba.maxSize(compute_fine_max_grid_size(refine_ratio, max_grid_size));
    amrex::Print() << "Number of fine boxes after merging: " << new_fine_ba.size() << "\n";
    const Box coarse_domain = coarse_geom.Domain();
    const Box fine_domain = amrex::refine(coarse_domain, IntVect(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio)));

    // Clip to fine_domain without relying on BoxArray::operator&=(Box),
    // which is not available in some AMReX versions.
    BoxList clipped_bl;
    for (int i = 0; i < new_fine_ba.size(); ++i) {
        Box clipped = new_fine_ba[i] & fine_domain;
        if (clipped.ok()) {
            clipped_bl.push_back(clipped);
        }
    }
    new_fine_ba = BoxArray(clipped_bl);
    amrex::Print() << "Number of fine boxes after clipping to fine domain: " << new_fine_ba.size() << "\n";
    if (new_fine_ba.size() == 0) {
        if (has_fine_level) {
            has_fine_level = false;
            fine_geom.reset();
            fine_phi.reset();
            fine_vel.reset();
            amrex::Print() << "Removed fine level (no valid tagged boxes)\n";
            return true;
        }
        return false;
    }
    amrex::Print() << "Creating fine level with domain " << fine_domain << " and " << new_fine_ba.size() << " boxes\n";
    const auto old_prob_lo = coarse_geom.ProbLoArray();
    const auto old_prob_hi = coarse_geom.ProbHiArray();
    const auto old_periodic = coarse_geom.isPeriodicArray();
    const int coord = coarse_geom.Coord();
    RealBox rb({AMREX_D_DECL(old_prob_lo[0], old_prob_lo[1], old_prob_lo[2])},
               {AMREX_D_DECL(old_prob_hi[0], old_prob_hi[1], old_prob_hi[2])});
    Geometry new_fine_geom(fine_domain, &rb, coord, old_periodic.data());
    DistributionMapping new_fine_dm(new_fine_ba);
    amrex::Print() << "New fine level geometry: domain = " << new_fine_geom.Domain()
                   << ", prob_lo = " << new_fine_geom.ProbLo()
                   << ", prob_hi = " << new_fine_geom.ProbHi() << "\n";
    std::unique_ptr<MultiFab> old_fine_phi = std::move(fine_phi);
    std::unique_ptr<MultiFab> old_fine_vel = std::move(fine_vel);

    fine_phi = std::make_unique<MultiFab>(new_fine_ba, new_fine_dm, 1, ng_phi);
    fine_vel = std::make_unique<MultiFab>(new_fine_ba, new_fine_dm, 3, 1);
    initialize_fine_from_coarse(coarse_phi, coarse_vel, *fine_phi, *fine_vel, refine_ratio);

    if (old_fine_phi && old_fine_vel) {
        // Preserve existing fine-level values in regions that remain refined.
        fine_phi->ParallelCopy(*old_fine_phi, 0, 0, 1, 0, 0);
        fine_vel->ParallelCopy(*old_fine_vel, 0, 0, 3, 0, 0);
    }

    fine_ba = new_fine_ba;
    fine_dm = new_fine_dm;
    fine_geom = std::make_unique<Geometry>(new_fine_geom);
    has_fine_level = true;

    amrex::Print() << "Updated local fine level for phi < " << tag_phi_threshold
                   << " with " << fine_ba.size() << " boxes\n";
    return true;
}

#endif