// ============================================================================
// Heat Per Unit Area (HPA) Raster Output
// ============================================================================
// Computes the cumulative heat released per unit ground area [BTU/ft²] for
// every burned cell following the FARSITE convention:
//
//   HPA [BTU/ft²] = I_R [BTU/ft²·min] × τ [min]
//
// where
//   I_R  = Rothermel (1972) reaction intensity [BTU/ft²·min]
//   τ    = residence time [min]  (tau_residence [s] / 60)
//
// HPA is set to 0 for unburned cells (phi ≥ 0) and to the computed value for
// all burned cells (phi < 0).  This field is written to every plotfile as
// "heat_per_unit_area".
//
// Physical background
// -------------------
// Reaction intensity I_R [BTU/ft²·min] quantifies the rate of heat release
// per unit area of the fuel bed at the fire front.  Multiplied by the
// residence time (the duration a given location experiences active combustion),
// it gives the total energy deposited per unit area — equivalent to FARSITE's
// "heat per unit area" output and the numerator of Byram's (1959) fireline
// intensity formula when expressed per unit perimeter length.
//
// References:
//   Rothermel (1972). USDA Forest Service Research Paper INT-115.
//   Byram (1959). "Combustion of forest fuels." In Forest Fire: Control and
//     Use, pp. 61-89. McGraw-Hill.
//   Finney (1998). FARSITE: Fire Area Simulator — Model Development and
//     Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef HEAT_PER_UNIT_AREA_H
#define HEAT_PER_UNIT_AREA_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// compute_heat_per_unit_area
//
// Parameters:
//   hpa_mf        – output: 1-component MultiFab [BTU/ft²]
//   phi           – input:  level-set field (phi < 0 = burned)
//   I_R           – reaction intensity [BTU/ft²·min] (from Rothermel model)
//   tau_residence – fire residence time [s]
// ============================================================================
inline void compute_heat_per_unit_area(MultiFab&       hpa_mf,
                                       const MultiFab& phi,
                                       amrex::Real     I_R,
                                       amrex::Real     tau_residence)
{
    // Convert residence time from seconds to minutes and compute HPA
    const Real hpa_val = I_R * (tau_residence / Real(60.0));

    for (MFIter mfi(hpa_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto       hpa = hpa_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            hpa(i, j, k) = (p(i, j, k) < Real(0.0)) ? hpa_val : Real(0.0);
        });
    }
}

#endif // HEAT_PER_UNIT_AREA_H
