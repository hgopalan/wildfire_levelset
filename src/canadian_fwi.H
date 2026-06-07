// ============================================================================
// Canadian Forest Fire Weather Index (FWI) System
// ============================================================================
// Implements the complete Canadian FWI System (Van Wagner 1987) including:
//   - Fine Fuel Moisture Code (FFMC) — 1-2 hour timelag fuel moisture
//   - Duff Moisture Code (DMC) — 10-20 day duff layer moisture
//   - Drought Code (DC) — seasonal deep organic layer moisture
//   - Buildup Index (BUI) — combined DMC and DC
//   - Fire Weather Index (FWI) — overall fire danger rating
//
// The FWI System is used operationally in Canada and many other countries
// for fire danger rating. It provides a standardized framework for assessing
// fire weather conditions based on temperature, relative humidity, wind speed,
// and 24-hour precipitation.
//
// References:
//   Van Wagner, C.E. (1987). Development and structure of the Canadian Forest
//     Fire Weather Index System. Canadian Forestry Service Technical Report 35.
//   Van Wagner, C.E. & Pickett, T.L. (1985). Equations and FORTRAN program for
//     the Canadian Forest Fire Weather Index System. Canadian Forestry Service
//     Technical Report 33.
//   Lawson, B.D. & Armitage, O.B. (2008). Weather Guide for the Canadian Forest
//     Fire Danger Rating System. Natural Resources Canada.
// ============================================================================

#ifndef CANADIAN_FWI_H
#define CANADIAN_FWI_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// FWI State Structure
// Stores the current state of all FWI System codes
// ============================================================================
struct FWIState {
    Real FFMC;  // Fine Fuel Moisture Code (dimensionless, typically 0-101)
    Real DMC;   // Duff Moisture Code (dimensionless, typically 0-300+)
    Real DC;    // Drought Code (dimensionless, typically 0-800+)
    Real ISI;   // Initial Spread Index (dimensionless)
    Real BUI;   // Buildup Index (dimensionless)
    Real FWI;   // Fire Weather Index (dimensionless)
};

// ============================================================================
// Default initial conditions for FWI codes
// These are typical spring startup values for temperate forests
// ============================================================================
namespace FWIConst {
    constexpr Real FFMC_DEFAULT = 85.0;  // Typical dry-day startup
    constexpr Real DMC_DEFAULT  = 6.0;   // Spring startup
    constexpr Real DC_DEFAULT   = 15.0;  // Spring startup
}

// ============================================================================
// compute_ffmc
//
// Fine Fuel Moisture Code (FFMC): represents moisture content of litter and
// fine fuels (1-2 hour timelag). Responds rapidly to weather changes.
//
// Van Wagner (1987) equations:
//   m_o = 147.2 × (101 - FFMC_prev) / (59.5 + FFMC_prev)  [%]
//   (apply wetting from rain)
//   (apply drying from EMC)
//   FFMC = 59.5 × (250 - m) / (147.2 + m)
//
// Parameters:
//   FFMC_prev - previous day's FFMC
//   T         - temperature [°C]
//   RH        - relative humidity [%] (0-100)
//   W         - wind speed [km/h]
//   rain      - 24-hour precipitation [mm]
//
// Returns: updated FFMC
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_ffmc(Real FFMC_prev, Real T, Real RH, Real W, Real rain)
{
    // Convert previous FFMC to moisture content [%]
    Real m_o = Real(147.2) * (Real(101.0) - FFMC_prev) / (Real(59.5) + FFMC_prev);
    
    // Apply wetting if rain > 0.5 mm
    Real m_r = m_o;
    if (rain > Real(0.5)) {
        if (m_o <= Real(150.0)) {
            m_r = m_o + Real(42.5) * rain * std::exp(Real(-100.0) / (Real(251.0) - m_o))
                      * (Real(1.0) - std::exp(Real(-6.93) / rain));
        } else {
            m_r = m_o + Real(42.5) * rain * std::exp(Real(-100.0) / (Real(251.0) - m_o))
                      * (Real(1.0) - std::exp(Real(-6.93) / rain))
                      + Real(0.0015) * (m_o - Real(150.0)) * (m_o - Real(150.0))
                        * std::sqrt(rain);
        }
        if (m_r > Real(250.0)) m_r = Real(250.0);
    }
    
    // Equilibrium Moisture Content (EMC) for drying
    Real E_d;
    if (RH < Real(100.0)) {
        const Real RH_safe = amrex::max(RH, Real(1.0));  // Prevent log(0)
        E_d = Real(0.942) * std::pow(RH_safe, Real(0.679))
            + Real(11.0) * std::exp((RH_safe - Real(100.0)) / Real(10.0))
            + Real(0.18) * (Real(21.1) - T) * (Real(1.0) - std::exp(Real(-0.115) * RH_safe));
    } else {
        E_d = m_r;  // No drying at 100% RH
    }
    
    // Logarithmic drying rate (function of wind and EMC)
    Real k_l, k_w;
    if (m_r > E_d) {
        k_l = Real(0.424) * (Real(1.0) - std::pow((RH / Real(100.0)), Real(1.7)))
            + Real(0.0694) * std::sqrt(W) * (Real(1.0) - std::pow((RH / Real(100.0)), Real(8.0)));
        k_w = k_l * (Real(0.581) * std::exp(Real(0.0365) * T));
        const Real m_1 = E_d + (m_r - E_d) * std::exp(Real(-2.303) * k_w);
        
        // Apply drying
        m_r = m_1;
    }
    
    // Equilibrium Moisture Content (EMC) for wetting
    Real E_w;
    if (RH < Real(100.0)) {
        const Real RH_safe = amrex::max(RH, Real(1.0));
        E_w = Real(0.618) * std::pow(RH_safe, Real(0.753))
            + Real(10.0) * std::exp((RH_safe - Real(100.0)) / Real(10.0))
            + Real(0.18) * (Real(21.1) - T) * (Real(1.0) - std::exp(Real(-0.115) * RH_safe));
    } else {
        E_w = m_r;  // No wetting beyond saturation
    }
    
    // Logarithmic wetting rate
    if (m_r < E_w) {
        k_l = Real(0.424) * (Real(1.0) - std::pow(((Real(100.0) - RH) / Real(100.0)), Real(1.7)))
            + Real(0.0694) * std::sqrt(W) * (Real(1.0) - std::pow(((Real(100.0) - RH) / Real(100.0)), Real(8.0)));
        k_w = k_l * (Real(0.581) * std::exp(Real(0.0365) * T));
        const Real m_2 = E_w - (E_w - m_r) * std::exp(Real(-2.303) * k_w);
        
        // Apply wetting
        m_r = m_2;
    }
    
    // Convert moisture content back to FFMC
    Real FFMC = Real(59.5) * (Real(250.0) - m_r) / (Real(147.2) + m_r);
    
    // Clamp to valid range
    if (FFMC > Real(101.0)) FFMC = Real(101.0);
    if (FFMC < Real(0.0))   FFMC = Real(0.0);
    
    return FFMC;
}

// ============================================================================
// compute_dmc
//
// Duff Moisture Code (DMC): represents moisture content of loosely compacted
// organic layers of moderate depth (10-20 day timelag). Responds more slowly
// than FFMC.
//
// Van Wagner (1987) equations:
//   P_r = (rain > 1.5 mm) ? (rain - 1.5) : 0
//   r_e = 0.92 × rain - 1.27
//   M_o = 20 + exp(5.6348 - DMC_prev / 43.43)
//   (apply wetting)
//   M = 244.72 - 43.43 × ln(M - 20)
//   K = 1.894 × (T + 1.1) × (100 - RH) × L_e × 10^(-6)
//   DMC = DMC_prev + 100 × K  (if no rain)
//
// Parameters:
//   DMC_prev - previous day's DMC
//   T        - temperature [°C]
//   RH       - relative humidity [%]
//   rain     - 24-hour precipitation [mm]
//   month    - month (1-12) for day-length factor
//
// Returns: updated DMC
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_dmc(Real DMC_prev, Real T, Real RH, Real rain, int month)
{
    // Day-length factor (average for mid-latitudes, 46°N)
    // Values from Van Wagner (1987) Table 3
    constexpr Real L_e_table[12] = {
        6.5,  // Jan
        7.5,  // Feb
        9.0,  // Mar
        12.8, // Apr
        13.9, // May
        13.9, // Jun
        12.4, // Jul
        10.9, // Aug
        9.4,  // Sep
        8.0,  // Oct
        7.0,  // Nov
        6.0   // Dec
    };
    
    const int m_idx = amrex::max(0, amrex::min(11, month - 1));
    const Real L_e = L_e_table[m_idx];
    
    Real DMC = DMC_prev;
    
    // Apply wetting if rain > 1.5 mm
    if (rain > Real(1.5)) {
        const Real r_e = Real(0.92) * rain - Real(1.27);
        const Real M_o = Real(20.0) + std::exp(Real(5.6348) - DMC_prev / Real(43.43));
        Real b;
        if (DMC_prev <= Real(33.0)) {
            b = Real(100.0) / (Real(0.5) + Real(0.3) * DMC_prev);
        } else if (DMC_prev <= Real(65.0)) {
            b = Real(14.0) - Real(1.3) * std::log(DMC_prev);
        } else {
            b = Real(6.2) * std::log(DMC_prev) - Real(17.2);
        }
        const Real M_r = M_o + Real(1000.0) * r_e / (Real(48.77) + b * r_e);
        const Real P_r = Real(244.72) - Real(43.43) * std::log(M_r - Real(20.0));
        if (P_r < Real(0.0)) {
            DMC = Real(0.0);
        } else {
            DMC = P_r;
        }
    }
    
    // Apply drying
    if (T > Real(-1.1)) {
        const Real T_safe = amrex::max(T, Real(-1.1));
        const Real RH_safe = amrex::min(RH, Real(100.0));
        const Real K = Real(1.894) * (T_safe + Real(1.1)) * (Real(100.0) - RH_safe)
                     * L_e * Real(1.0e-6);
        DMC = DMC + Real(100.0) * K;
    }
    
    // Clamp to valid range
    if (DMC < Real(0.0)) DMC = Real(0.0);
    
    return DMC;
}

// ============================================================================
// compute_dc
//
// Drought Code (DC): represents moisture content of deep, compact organic
// layers (seasonal timelag, 50+ days). Very slow response to weather.
//
// Van Wagner (1987) equations:
//   Q_o = 800 × exp(-DC_prev / 400)
//   (apply wetting from rain)
//   Q = Q_o - 3.937 × r_d
//   DC = 400 × ln(800 / Q)
//   V = 0.36 × (T + 2.8) + L_f
//   D = V × 10^(-6) / 1.6
//   DC = DC_prev + 100 × D
//
// Parameters:
//   DC_prev - previous day's DC
//   T       - temperature [°C]
//   rain    - 24-hour precipitation [mm]
//   month   - month (1-12) for day-length factor
//
// Returns: updated DC
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_dc(Real DC_prev, Real T, Real rain, int month)
{
    // Day-length factor (average for mid-latitudes, 46°N)
    // Values from Van Wagner (1987) Table 4
    constexpr Real L_f_table[12] = {
        -1.6, // Jan
        -1.6, // Feb
        -1.6, // Mar
        0.9,  // Apr
        3.8,  // May
        5.8,  // Jun
        6.4,  // Jul
        5.0,  // Aug
        2.4,  // Sep
        0.4,  // Oct
        -1.6, // Nov
        -1.6  // Dec
    };
    
    const int m_idx = amrex::max(0, amrex::min(11, month - 1));
    const Real L_f = L_f_table[m_idx];
    
    Real DC = DC_prev;
    
    // Apply wetting if rain > 2.8 mm
    if (rain > Real(2.8)) {
        const Real r_d = Real(0.83) * rain - Real(1.27);
        const Real Q_o = Real(800.0) * std::exp(-DC_prev / Real(400.0));
        const Real Q_r = Q_o + Real(3.937) * r_d;
        const Real D_r = Real(400.0) * std::log(Real(800.0) / Q_r);
        if (D_r < Real(0.0)) {
            DC = Real(0.0);
        } else {
            DC = D_r;
        }
    }
    
    // Apply drying
    if (T > Real(-2.8)) {
        const Real T_safe = amrex::max(T, Real(-2.8));
        const Real V = Real(0.36) * (T_safe + Real(2.8)) + L_f;
        if (V > Real(0.0)) {
            const Real D = V * Real(1.0e-6) / Real(1.6);
            DC = DC + Real(100.0) * D;
        }
    }
    
    // Clamp to valid range
    if (DC < Real(0.0)) DC = Real(0.0);
    
    return DC;
}

// ============================================================================
// compute_isi
//
// Initial Spread Index (ISI): combines FFMC and wind speed to represent
// rate of fire spread potential.
//
// Van Wagner (1987) equations:
//   m = 147.2 × (101 - FFMC) / (59.5 + FFMC)
//   f_F = 91.9 × exp(-0.1386 × m) × (1 + m^5.31 / 4.93×10^7)
//   f_W = exp(0.05039 × W)
//   ISI = 0.208 × f_F × f_W
//
// Parameters:
//   FFMC - Fine Fuel Moisture Code
//   W    - wind speed [km/h]
//
// Returns: Initial Spread Index
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_isi(Real FFMC, Real W)
{
    // Convert FFMC to moisture content [%]
    const Real m = Real(147.2) * (Real(101.0) - FFMC) / (Real(59.5) + FFMC);
    
    // Fine fuel moisture function
    const Real f_F = Real(91.9) * std::exp(Real(-0.1386) * m)
                   * (Real(1.0) + std::pow(m, Real(5.31)) / Real(4.93e7));
    
    // Wind function
    const Real f_W = std::exp(Real(0.05039) * W);
    
    // Initial Spread Index
    const Real ISI = Real(0.208) * f_F * f_W;
    
    return ISI;
}

// ============================================================================
// compute_bui
//
// Buildup Index (BUI): combines DMC and DC to represent total fuel available
// for combustion.
//
// Van Wagner (1987) equations:
//   if DMC <= 0.4 × DC:
//     BUI = 0.8 × DMC × DC / (DMC + 0.4 × DC)
//   else:
//     BUI = DMC - (1 - 0.8 × DC / (DMC + 0.4 × DC)) × (0.92 + (0.0114 × DMC)^1.7)
//
// Parameters:
//   DMC - Duff Moisture Code
//   DC  - Drought Code
//
// Returns: Buildup Index
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_bui(Real DMC, Real DC)
{
    Real BUI;
    if (DMC <= Real(0.4) * DC) {
        BUI = Real(0.8) * DMC * DC / (DMC + Real(0.4) * DC);
    } else {
        BUI = DMC - (Real(1.0) - Real(0.8) * DC / (DMC + Real(0.4) * DC))
            * (Real(0.92) + std::pow((Real(0.0114) * DMC), Real(1.7)));
    }
    
    // Clamp to valid range
    if (BUI < Real(0.0)) BUI = Real(0.0);
    
    return BUI;
}

// ============================================================================
// compute_fwi
//
// Fire Weather Index (FWI): combines ISI and BUI to represent overall fire
// intensity potential.
//
// Van Wagner (1987) equations:
//   if BUI <= 80:
//     f_D = 0.626 × BUI^0.809 + 2
//   else:
//     f_D = 1000 / (25 + 108.64 × exp(-0.023 × BUI))
//   B = 0.1 × ISI × f_D
//   if B > 1:
//     S = exp(2.72 × (0.434 × ln(B))^0.647)
//     FWI = S
//   else:
//     FWI = B
//
// Parameters:
//   ISI - Initial Spread Index
//   BUI - Buildup Index
//
// Returns: Fire Weather Index
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_fwi(Real ISI, Real BUI)
{
    // Fuel moisture function
    Real f_D;
    if (BUI <= Real(80.0)) {
        f_D = Real(0.626) * std::pow(BUI, Real(0.809)) + Real(2.0);
    } else {
        f_D = Real(1000.0) / (Real(25.0) + Real(108.64) * std::exp(Real(-0.023) * BUI));
    }
    
    // Intermediate value
    const Real B = Real(0.1) * ISI * f_D;
    
    // Fire Weather Index
    Real FWI;
    if (B > Real(1.0)) {
        const Real ln_B = std::log(B);
        const Real S = std::exp(Real(2.72) * std::pow((Real(0.434) * ln_B), Real(0.647)));
        FWI = S;
    } else {
        FWI = B;
    }
    
    return FWI;
}

// ============================================================================
// update_fwi_system
//
// Updates the complete FWI system for one timestep (one day).
//
// Parameters:
//   state_prev - previous day's FWI state
//   T          - temperature [°C]
//   RH         - relative humidity [%] (0-100)
//   W          - wind speed [km/h]
//   rain       - 24-hour precipitation [mm]
//   month      - month (1-12) for seasonal day-length adjustments
//
// Returns: updated FWI state
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
FWIState update_fwi_system(const FWIState& state_prev,
                           Real T, Real RH, Real W, Real rain, int month)
{
    FWIState state;
    
    // Update moisture codes
    state.FFMC = compute_ffmc(state_prev.FFMC, T, RH, W, rain);
    state.DMC  = compute_dmc(state_prev.DMC, T, RH, rain, month);
    state.DC   = compute_dc(state_prev.DC, T, rain, month);
    
    // Compute spread and buildup indices
    state.ISI = compute_isi(state.FFMC, W);
    state.BUI = compute_bui(state.DMC, state.DC);
    
    // Compute overall Fire Weather Index
    state.FWI = compute_fwi(state.ISI, state.BUI);
    
    return state;
}

// ============================================================================
// compute_fwi_moisture_codes
//
// Computes FWI system outputs and stores them in a MultiFab.
// Useful for fire danger mapping and fuel moisture scheduling.
//
// MultiFab layout:
//   Component 0 - FFMC (Fine Fuel Moisture Code)
//   Component 1 - DMC (Duff Moisture Code)
//   Component 2 - DC (Drought Code)
//   Component 3 - ISI (Initial Spread Index)
//   Component 4 - BUI (Buildup Index)
//   Component 5 - FWI (Fire Weather Index)
//
// Parameters:
//   fwi_mf     - output: 6-component MultiFab
//   T          - temperature [°C]
//   RH         - relative humidity [%]
//   W          - wind speed [km/h]
//   rain       - 24-hour precipitation [mm]
//   month      - month (1-12)
//   FFMC_init  - initial FFMC (default: 85.0)
//   DMC_init   - initial DMC (default: 6.0)
//   DC_init    - initial DC (default: 15.0)
// ============================================================================
inline void compute_fwi_moisture_codes(
    MultiFab& fwi_mf,
    Real T, Real RH, Real W, Real rain, int month,
    Real FFMC_init = FWIConst::FFMC_DEFAULT,
    Real DMC_init  = FWIConst::DMC_DEFAULT,
    Real DC_init   = FWIConst::DC_DEFAULT)
{
    // Initialize previous state
    FWIState state_prev;
    state_prev.FFMC = FFMC_init;
    state_prev.DMC  = DMC_init;
    state_prev.DC   = DC_init;
    state_prev.ISI  = Real(0.0);
    state_prev.BUI  = Real(0.0);
    state_prev.FWI  = Real(0.0);
    
    // Update FWI system
    FWIState state = update_fwi_system(state_prev, T, RH, W, rain, month);
    
    // Fill MultiFab
    for (MFIter mfi(fwi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto fwi_arr = fwi_mf.array(mfi);
        
        const Real ffmc = state.FFMC;
        const Real dmc  = state.DMC;
        const Real dc   = state.DC;
        const Real isi  = state.ISI;
        const Real bui  = state.BUI;
        const Real fwi  = state.FWI;
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            fwi_arr(i, j, k, 0) = ffmc;  // FFMC
            fwi_arr(i, j, k, 1) = dmc;   // DMC
            fwi_arr(i, j, k, 2) = dc;    // DC
            fwi_arr(i, j, k, 3) = isi;   // ISI
            fwi_arr(i, j, k, 4) = bui;   // BUI
            fwi_arr(i, j, k, 5) = fwi;   // FWI
        });
    }
}

// ============================================================================
// Spatially-varying FWI update dynamically at each timestep
//
// Parameters:
//   fwi_mf         - input/output: FWI System codes MultiFab (6 components)
//   temperature_mf - temperature field [°C]
//   humidity_mf    - relative humidity field [%]
//   vel            - wind velocity field [m/s]
//   rain_rate_mm_hr- current rainfall rate [mm/hr]
//   dt_seconds     - timestep [seconds]
//   month          - calendar month [1-12]
// ============================================================================
inline void update_fwi_field(
    MultiFab& fwi_mf,
    const MultiFab& temperature_mf,
    const MultiFab& humidity_mf,
    const MultiFab& vel,
    Real rain_rate_mm_hr,
    Real dt_seconds,
    int month)
{
    const Real dt_frac = dt_seconds / Real(86400.0);
    const Real rain = rain_rate_mm_hr * Real(24.0);
    const Real ms_to_kmh = Real(3.6);

    for (MFIter mfi(fwi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto fwi_arr = fwi_mf.array(mfi);
        auto const T_arr = temperature_mf.const_array(mfi);
        auto const RH_arr = humidity_mf.const_array(mfi);
        auto const vel_arr = vel.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Retrieve previous FWI state
            FWIState state_prev;
            state_prev.FFMC = fwi_arr(i, j, k, 0);
            state_prev.DMC  = fwi_arr(i, j, k, 1);
            state_prev.DC   = fwi_arr(i, j, k, 2);
            state_prev.ISI  = fwi_arr(i, j, k, 3);
            state_prev.BUI  = fwi_arr(i, j, k, 4);
            state_prev.FWI  = fwi_arr(i, j, k, 5);

            const Real T = T_arr(i, j, k);
            const Real RH = RH_arr(i, j, k);

            // Wind speed in km/h
            const Real ux = vel_arr(i, j, k, 0);
            const Real uy = vel_arr(i, j, k, 1);
            const Real wind_speed_ms = std::sqrt(ux*ux + uy*uy);
            const Real W = wind_speed_ms * ms_to_kmh;

            // Compute daily target state
            FWIState target = update_fwi_system(state_prev, T, RH, W, rain, month);

            // Fractionally update moisture codes
            const Real ffmc_new = state_prev.FFMC + (target.FFMC - state_prev.FFMC) * dt_frac;
            const Real dmc_new  = state_prev.DMC  + (target.DMC  - state_prev.DMC)  * dt_frac;
            const Real dc_new   = state_prev.DC   + (target.DC   - state_prev.DC)   * dt_frac;

            // Recompute diagnostic indices based on updated codes
            const Real isi_new = compute_isi(ffmc_new, W);
            const Real bui_new = compute_bui(dmc_new, dc_new);
            const Real fwi_new = compute_fwi(isi_new, bui_new);

            // Write back to MultiFab
            fwi_arr(i, j, k, 0) = ffmc_new;
            fwi_arr(i, j, k, 1) = dmc_new;
            fwi_arr(i, j, k, 2) = dc_new;
            fwi_arr(i, j, k, 3) = isi_new;
            fwi_arr(i, j, k, 4) = bui_new;
            fwi_arr(i, j, k, 5) = fwi_new;
        });
    }
}

#endif // CANADIAN_FWI_H
