// ============================================================================
// Fire Intensity Class Raster (FARSITE-compatible)
// ============================================================================
// Classifies Byram (1959) fireline intensity into 6 discrete classes following
// the FARSITE / NFFL convention documented in Andrews et al. (2011) RMRS-GTR-249
// and used in the fire_report.H HTML output.
//
// Classes:
//   I:   < 10 kW/m        (< 0.5 m flame)    - Creeping / smouldering
//   II:  10-175 kW/m      (0.5-1.5 m flame)  - Low intensity surface fire
//   III: 175-500 kW/m     (1.5-2.5 m flame)  - Moderate intensity
//   IV:  500-2000 kW/m    (2.5-5 m flame)    - High intensity; spotting begins
//   V:   2000-10000 kW/m  (5-11 m flame)     - Very high; active crown fire
//   VI:  > 10000 kW/m     (> 11 m flame)     - Extreme; firestorm
//
// References:
//   Byram, G.M. (1959). Combustion of forest fuels. In: Davis, K.P. (ed.)
//     Forest Fire: Control and Use. McGraw-Hill.
//   Andrews, P.L. et al. (2011). BehavePlus fire modeling system, version 5.0:
//     Design and features. USDA Forest Service RMRS-GTR-249.
// ============================================================================

#ifndef FIRE_INTENSITY_CLASS_H
#define FIRE_INTENSITY_CLASS_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Classify Byram fireline intensity [kW/m] into 6 discrete classes (1-6)
// Returns integer class code:
//   1 = Class I    (< 10 kW/m)
//   2 = Class II   (10-175 kW/m)
//   3 = Class III  (175-500 kW/m)
//   4 = Class IV   (500-2000 kW/m)
//   5 = Class V    (2000-10000 kW/m)
//   6 = Class VI   (> 10000 kW/m)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int classify_fire_intensity(Real I_B_kW)
{
    if (I_B_kW < Real(10.0)) {
        return 1; // Class I: creeping / smouldering
    } else if (I_B_kW < Real(175.0)) {
        return 2; // Class II: low intensity surface fire
    } else if (I_B_kW < Real(500.0)) {
        return 3; // Class III: moderate intensity
    } else if (I_B_kW < Real(2000.0)) {
        return 4; // Class IV: high intensity; spotting begins
    } else if (I_B_kW < Real(10000.0)) {
        return 5; // Class V: very high; active crown fire
    } else {
        return 6; // Class VI: extreme; firestorm
    }
}

// ============================================================================
// Compute fire intensity class for entire MultiFab
//
// Parameters:
//   fire_intensity_class_mf – output: integer class code (1-6) per cell
//   fireline_intensity_mf   – input:  Byram fireline intensity [kW/m]
// ============================================================================
inline void compute_fire_intensity_class(
    MultiFab&       fire_intensity_class_mf,
    const MultiFab& fireline_intensity_mf)
{
    for (MFIter mfi(fire_intensity_class_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ic_arr = fire_intensity_class_mf.array(mfi);
        auto fi_arr = fireline_intensity_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real I_B = fi_arr(i, j, k);
            ic_arr(i, j, k) = static_cast<Real>(classify_fire_intensity(I_B));
        });
    }
}

#endif // FIRE_INTENSITY_CLASS_H
