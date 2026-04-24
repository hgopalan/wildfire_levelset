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
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

static int compute_fine_max_grid_size(int refine_ratio, int max_grid_size);

namespace
{
    constexpr int box_int_count = 2 * AMREX_SPACEDIM;

    void append_box_bounds(std::vector<int>& box_data, const Box& bx)
    {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            box_data.push_back(bx.smallEnd(dir));
        }
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            box_data.push_back(bx.bigEnd(dir));
        }
    }

    Box box_from_bounds(const std::vector<int>& box_data, int offset)
    {
        IntVect lo = IntVect::TheZeroVector();
        IntVect hi = IntVect::TheZeroVector();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            lo[dir] = box_data[offset + dir];
            hi[dir] = box_data[offset + AMREX_SPACEDIM + dir];
        }
        return Box(lo, hi);
    }

    BoxArray build_global_fine_boxarray(const BoxList& local_tagged_cells,
                                        int refine_ratio,
                                        int max_grid_size,
                                        const Box& fine_domain)
    {
        if (amrex::ParallelDescriptor::NProcs() == 1) {
            BoxList fine_bl;
            for (const Box& bx : local_tagged_cells) {
                fine_bl.push_back(amrex::refine(bx,
                                                IntVect(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio))));
            }

            BoxArray fine_ba(fine_bl);
            fine_ba.maxSize(compute_fine_max_grid_size(refine_ratio, max_grid_size));

            BoxList clipped_bl;
            for (int i = 0; i < fine_ba.size(); ++i) {
                Box clipped = fine_ba[i] & fine_domain;
                if (clipped.ok()) {
                    clipped_bl.push_back(clipped);
                }
            }

            return BoxArray(clipped_bl);
        }

        std::vector<int> local_box_data;
        local_box_data.reserve(local_tagged_cells.size() * box_int_count);
        for (const Box& bx : local_tagged_cells) {
            append_box_bounds(local_box_data, bx);
        }

        const int io_rank = amrex::ParallelDescriptor::IOProcessorNumber();
        std::vector<int> recv_counts =
            amrex::ParallelDescriptor::Gather(static_cast<int>(local_box_data.size()), io_rank);
        std::vector<int> displacements;
        std::vector<int> gathered_box_data;

        if (amrex::ParallelDescriptor::IOProcessor()) {
            displacements.resize(recv_counts.size(), 0);
            int total_count = 0;
            for (int rank_idx = 0; rank_idx < static_cast<int>(recv_counts.size()); ++rank_idx) {
                displacements[rank_idx] = total_count;
                total_count += recv_counts[rank_idx];
            }
            gathered_box_data.resize(total_count);
        }

        amrex::ParallelDescriptor::Gatherv(local_box_data.data(),
                                           static_cast<int>(local_box_data.size()),
                                           gathered_box_data.data(),
                                           recv_counts,
                                           displacements,
                                           io_rank);

        std::vector<int> global_box_data;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            BoxList fine_bl;
            for (int offset = 0; offset < static_cast<int>(gathered_box_data.size()); offset += box_int_count) {
                fine_bl.push_back(amrex::refine(box_from_bounds(gathered_box_data, offset),
                                                IntVect(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio))));
            }

            BoxArray fine_ba(fine_bl);
            fine_ba.maxSize(compute_fine_max_grid_size(refine_ratio, max_grid_size));

            BoxList clipped_bl;
            for (int i = 0; i < fine_ba.size(); ++i) {
                Box clipped = fine_ba[i] & fine_domain;
                if (clipped.ok()) {
                    clipped_bl.push_back(clipped);
                }
            }

            BoxArray clipped_ba(clipped_bl);
            global_box_data.reserve(clipped_ba.size() * box_int_count);
            for (int i = 0; i < clipped_ba.size(); ++i) {
                append_box_bounds(global_box_data, clipped_ba[i]);
            }
        }

        int global_box_data_size = static_cast<int>(global_box_data.size());
        amrex::ParallelDescriptor::Bcast(&global_box_data_size, 1, io_rank);
        if (!amrex::ParallelDescriptor::IOProcessor()) {
            global_box_data.resize(global_box_data_size);
        }
        if (global_box_data_size > 0) {
            amrex::ParallelDescriptor::Bcast(global_box_data.data(), global_box_data_size, io_rank);
        }

        BoxList final_bl;
        for (int offset = 0; offset < global_box_data_size; offset += box_int_count) {
            final_bl.push_back(box_from_bounds(global_box_data, offset));
        }

        return BoxArray(final_bl);
    }
}

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
    int counter=0;
    for (MFIter mfi(coarse_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& phi_arr = coarse_phi.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (phi_arr(i, j, k, 0) < tag_phi_threshold) {
                // Tag only this individual cell, not the entire box
                tagged_cells.push_back(Box(IntVect(AMREX_D_DECL(i, j, k)),
                                        IntVect(AMREX_D_DECL(i, j, k))));
                counter=counter+1;
            }
        });
    }


    // Each rank tags only its own boxes (rank-local).
    /*BoxList tagged_cells;
    int counter=0;
    for (MFIter mfi(coarse_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (coarse_phi[mfi].min(bx, 0) < tag_phi_threshold) {
            counter=counter+1;
            tagged_cells.push_back(bx);
        }
    }*/
    amrex::Print() << "Number of tagged coarse boxes (local): " << counter << "\n";
    // Globally reduce the tag count so ALL ranks make the same branching decision.
    // Without this, ranks with no local tags would return early while others
    // proceed to collective MPI calls (e.g. DistributionMapping), causing a deadlock.
    int num_tagged = static_cast<int>(tagged_cells.size());
    amrex::ParallelDescriptor::ReduceIntSum(num_tagged);
    amrex::Print() << "Number of tagged coarse boxes (global): " << num_tagged << "\n";

    if (num_tagged == 0) {
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

    amrex::Print() << "Refining " << num_tagged << " coarse boxes (global) with refine ratio " << refine_ratio << "\n";
    const Box coarse_domain = coarse_geom.Domain();
    const Box fine_domain = amrex::refine(coarse_domain, IntVect(AMREX_D_DECL(refine_ratio, refine_ratio, refine_ratio)));
    BoxArray new_fine_ba = build_global_fine_boxarray(tagged_cells,
                                                      refine_ratio,
                                                      max_grid_size,
                                                      fine_domain);
    amrex::Print() << "Number of fine boxes after global merge/clip: " << new_fine_ba.size() << "\n";

    // Globally reduce to check if any valid boxes remain after clipping.
    int num_clipped = static_cast<int>(new_fine_ba.size());
    amrex::ParallelDescriptor::ReduceIntMax(num_clipped);

    if (num_clipped == 0) {
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
                << ", prob_lo = ("
                << new_fine_geom.ProbLo(0) << ", "
                << new_fine_geom.ProbLo(1) << ", "
                << new_fine_geom.ProbLo(2) << ")"
                << ", prob_hi = ("
                << new_fine_geom.ProbHi(0) << ", "
                << new_fine_geom.ProbHi(1) << ", "
                << new_fine_geom.ProbHi(2) << ")\n";

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
