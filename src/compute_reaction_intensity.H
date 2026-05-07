// ============================================================================
// Reaction Intensity Raster
// ============================================================================
// Computes a per-cell Rothermel (1972) reaction intensity [kW/m²] raster for
// output to plotfiles and GIS export.
//
// Reaction intensity I_R is the rate of heat release per unit area of the fuel
// bed [BTU/ft²/min in Rothermel's system; converted here to kW/m²].  It drives
// ROS through the propagating flux ratio ξ and is a key determinant of fireline
// intensity and fire behaviour hazard.
//
// Per-cell behaviour
// ------------------
//   When a FARSITE landscape file (.lcp) is active and a per-cell Rothermel
//   lookup table has been built, each cell is assigned the I_R of its fuel
//   model code.  When no landscape file is available (uniform fuel), all
//   cells receive the same scalar I_R computed from the global RothermelParams.
//
//   Crown-fire cells receive their surface I_R value (crown I_R is not
//   separately modelled here).
//
//   Cells with fuel model code 0 or outside the table range receive 0.
//
// Unit conversion
// ---------------
//   1 BTU/ft²/min = 0.18941 kW/m²
//
// References
// ----------
//   Rothermel, R.C. (1972). A mathematical model for predicting fire spread
//     in wildland fuels. USDA Forest Service Research Paper INT-115.
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator — Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef COMPUTE_REACTION_INTENSITY_H
#define COMPUTE_REACTION_INTENSITY_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

#include "parse_inputs.H"
#include "rothermel_model.H"

using namespace amrex;

// Unit conversion: BTU/ft²/min → kW/m²
static constexpr amrex::Real BTU_FT2_MIN_TO_KW_M2 = WildfireConst::BTU_FT2MIN_TO_KW_M2;

// ============================================================================
// compute_reaction_intensity
//
// Parameters:
//   reaction_intensity_mf – output: 1-component MultiFab [kW/m²]
//   fuel_model_mf         – input:  per-cell FBFM fuel model code
//                           (ignored when d_fuel_table_ptr is nullptr)
//   rp                    – global Rothermel parameters (used as fallback)
//   d_fuel_table_ptr      – per-fuel-model pre-computed Rothermel table on
//                           device (nullptr → use global rp)
//   fuel_table_size       – number of entries in d_fuel_table_ptr
// ============================================================================
inline void compute_reaction_intensity(
    MultiFab&                      reaction_intensity_mf,
    const MultiFab*                fuel_model_mf,
    const RothermelParams&         rp,
    const RothermelComputed*       d_fuel_table_ptr,
    int                            fuel_table_size)
{
    // Scalar reaction intensity from global parameters [kW/m²]
    const RothermelComputed rc_global = compute_rothermel_params(rp);
    const Real I_R_global_kw = rc_global.I_R * BTU_FT2_MIN_TO_KW_M2;

    const bool use_table = (d_fuel_table_ptr != nullptr && fuel_table_size > 0
                            && fuel_model_mf != nullptr);

    for (MFIter mfi(reaction_intensity_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ri = reaction_intensity_mf.array(mfi);

        if (!use_table) {
            // Uniform fuel: fill every cell with the global I_R
            const Real val = I_R_global_kw;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                ri(i, j, k) = val;
            });
        } else {
            auto const fm = fuel_model_mf->const_array(mfi);
            const RothermelComputed* tbl = d_fuel_table_ptr;
            const int tbl_size           = fuel_table_size;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                int code = static_cast<int>(fm(i, j, k));
                if (code >= 0 && code < tbl_size) {
                    ri(i, j, k) = tbl[code].I_R * BTU_FT2_MIN_TO_KW_M2;
                } else {
                    // Code out of range (non-burnable or unset): zero
                    ri(i, j, k) = Real(0.0);
                }
            });
        }
    }
}

#endif // COMPUTE_REACTION_INTENSITY_H
