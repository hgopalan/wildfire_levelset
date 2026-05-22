// ============================================================================
// Fuel Database for Rothermel Model
// ============================================================================
// Standard fuel models for continental United States based on Anderson (1982)
// "Aids to Determining Fuel Models For Estimating Fire Behavior"
// and Scott & Burgan (2005) "Standard Fire Behavior Fuel Models"
// ============================================================================

#ifndef FUEL_DATABASE_H
#define FUEL_DATABASE_H

#include "parse_inputs.H"
#include <map>
#include <string>
#include <algorithm>
#include <cctype>

// ============================================================================
// Anderson 13 Standard Fire Behavior Fuel Models
// ============================================================================
// These are the classical fuel models used across the continental US.
// Each model represents different fuel types and vegetation conditions.
// ============================================================================

struct FuelModel {
    std::string name;
    std::string description;
    amrex::Real w0;        // oven-dry total fuel load [lb/ft²]
    amrex::Real sigma;     // surface-area-to-volume ratio [ft⁻¹]
    amrex::Real delta;     // fuel bed depth [ft]
    amrex::Real M_x;       // moisture of extinction (fraction)
    amrex::Real h_heat;    // heat content [BTU/lb]
    amrex::Real S_T;       // total mineral content (fraction)
    amrex::Real S_e;       // effective mineral content (fraction)
    amrex::Real rho_p;     // particle density [lb/ft³]
    // Per-class fuel loads and SAV ratios (Rothermel 1972 multi-class path).
    // Default 0.0 – existing entries without per-class data retain the
    // single-class aggregate (w0/sigma) behaviour via the fallback path.
    // 10-hr dead SAV (109 ft⁻¹) and 100-hr dead SAV (30 ft⁻¹) are universal
    // constants and are not stored here; they are used directly in the model.
    amrex::Real w_d1    = 0.0;  // 1-hr  dead load [lb/ft²]
    amrex::Real sigma_d1 = 0.0;  // 1-hr  dead SAV  [ft⁻¹]
    amrex::Real w_d10   = 0.0;  // 10-hr dead load [lb/ft²]
    amrex::Real w_d100  = 0.0;  // 100-hr dead load [lb/ft²]
    amrex::Real w_lh    = 0.0;  // live herbaceous load [lb/ft²]
    amrex::Real sigma_lh = 0.0;  // live herbaceous SAV [ft⁻¹]
    amrex::Real w_lw    = 0.0;  // live woody load [lb/ft²]
    amrex::Real sigma_lw = 0.0;  // live woody SAV [ft⁻¹]
};

// ============================================================================
// Fuel Database - Anderson 13 Models
// ============================================================================
inline std::map<std::string, FuelModel> get_fuel_database() {
    std::map<std::string, FuelModel> database;
    
    // -------------------------------------------------------------------------
    // GRASS AND GRASS-DOMINATED MODELS (1-3)
    // -------------------------------------------------------------------------
    // Per-class loads from Rothermel (1972) Table A-1.
    // SAVs: 1-hr = model-specific; 10-hr = 109 ft⁻¹; 100-hr = 30 ft⁻¹;
    //       live herb/woody = 1500 ft⁻¹ (1550 ft⁻¹ for FM7).
    
    database["FM1"] = {
        "FM1", "Short Grass (1 ft)",
        0.034,    // w0 [lb/ft²]
        3500.0,   // sigma [ft⁻¹]
        1.0,      // delta [ft]
        0.12,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        // Per-class loads: 1-hr dead only
        0.034, 3500.0,   // w_d1, sigma_d1
        0.0,             // w_d10
        0.0,             // w_d100
        0.0, 1500.0,     // w_lh, sigma_lh
        0.0, 1500.0      // w_lw, sigma_lw
    };
    
    database["FM2"] = {
        "FM2", "Timber (Grass and Understory)",
        0.092,    // w0 [lb/ft²]
        3000.0,   // sigma [ft⁻¹]
        1.0,      // delta [ft]
        0.15,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.046, 3000.0,   // w_d1, sigma_d1
        0.023,           // w_d10
        0.023,           // w_d100
        0.023, 1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM3"] = {
        "FM3", "Tall Grass (2.5 ft)",
        0.138,    // w0 [lb/ft²]
        1500.0,   // sigma [ft⁻¹]
        2.5,      // delta [ft]
        0.25,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.138, 1500.0,   // w_d1, sigma_d1
        0.0,             // w_d10
        0.0,             // w_d100
        0.0,  1500.0,    // w_lh, sigma_lh
        0.0,  1500.0     // w_lw, sigma_lw
    };
    
    // -------------------------------------------------------------------------
    // CHAPARRAL AND SHRUB MODELS (4-7)
    // -------------------------------------------------------------------------
    
    database["FM4"] = {
        "FM4", "Chaparral (6 ft)",
        0.230,    // w0 [lb/ft²]
        1739.0,   // sigma [ft⁻¹]  (Rothermel 1972 Table A-1)
        6.0,      // delta [ft]
        0.20,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.230, 1739.0,   // w_d1, sigma_d1  (Rothermel 1972: σ₁ = 1739 ft⁻¹)
        0.184,           // w_d10
        0.092,           // w_d100
        0.230, 1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM5"] = {
        "FM5", "Brush (2 ft)",
        0.046,    // w0 [lb/ft²]
        1739.0,   // sigma [ft⁻¹]  (Rothermel 1972 Table A-1)
        2.0,      // delta [ft]
        0.20,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.046, 1739.0,   // w_d1, sigma_d1  (Rothermel 1972: σ₁ = 1739 ft⁻¹)
        0.023,           // w_d10
        0.0,             // w_d100
        0.092, 1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM6"] = {
        "FM6", "Dormant Brush, Hardwood Slash",
        0.069,    // w0 [lb/ft²]
        1750.0,   // sigma [ft⁻¹]
        2.5,      // delta [ft]
        0.25,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.069, 1750.0,   // w_d1, sigma_d1
        0.115,           // w_d10
        0.115,           // w_d100
        0.0,   1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM7"] = {
        "FM7", "Southern Rough",
        0.052,    // w0 [lb/ft²]
        1550.0,   // sigma [ft⁻¹]
        2.5,      // delta [ft]
        0.40,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.052, 1750.0,   // w_d1, sigma_d1
        0.086,           // w_d10
        0.069,           // w_d100
        0.017, 1550.0,   // w_lh, sigma_lh
        0.017, 1550.0    // w_lw, sigma_lw
    };
    
    // -------------------------------------------------------------------------
    // TIMBER LITTER MODELS (8-10)
    // -------------------------------------------------------------------------
    
    database["FM8"] = {
        "FM8", "Closed Timber Litter",
        0.069,    // w0 [lb/ft²]
        2000.0,   // sigma [ft⁻¹]
        0.2,      // delta [ft]
        0.30,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.069, 2000.0,   // w_d1, sigma_d1
        0.046,           // w_d10
        0.115,           // w_d100
        0.0,   1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM9"] = {
        "FM9", "Hardwood Litter",
        0.134,    // w0 [lb/ft²]
        2500.0,   // sigma [ft⁻¹]
        0.2,      // delta [ft]
        0.25,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.134, 2500.0,   // w_d1, sigma_d1
        0.019,           // w_d10
        0.007,           // w_d100
        0.0,   1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM10"] = {
        "FM10", "Timber (Litter and Understory)",
        0.138,    // w0 [lb/ft²]
        2000.0,   // sigma [ft⁻¹]
        1.0,      // delta [ft]
        0.25,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.138, 2000.0,   // w_d1, sigma_d1
        0.092,           // w_d10
        0.230,           // w_d100
        0.092, 1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    // -------------------------------------------------------------------------
    // SLASH MODELS (11-13)
    // -------------------------------------------------------------------------
    
    database["FM11"] = {
        "FM11", "Light Logging Slash",
        0.069,    // w0 [lb/ft²]
        1500.0,   // sigma [ft⁻¹]
        1.0,      // delta [ft]
        0.15,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.069, 1500.0,   // w_d1, sigma_d1
        0.207,           // w_d10
        0.253,           // w_d100
        0.0,   1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM12"] = {
        "FM12", "Medium Logging Slash",
        0.184,    // w0 [lb/ft²]
        1500.0,   // sigma [ft⁻¹]
        2.3,      // delta [ft]
        0.20,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.184, 1500.0,   // w_d1, sigma_d1
        0.644,           // w_d10
        0.759,           // w_d100
        0.0,   1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    database["FM13"] = {
        "FM13", "Heavy Logging Slash",
        0.322,    // w0 [lb/ft²]
        1500.0,   // sigma [ft⁻¹]
        3.0,      // delta [ft]
        0.25,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0,     // rho_p [lb/ft³]
        0.322, 1500.0,   // w_d1, sigma_d1
        1.058,           // w_d10
        1.288,           // w_d100
        0.0,   1500.0,   // w_lh, sigma_lh
        0.0,   1500.0    // w_lw, sigma_lw
    };
    
    // =========================================================================
    // Scott & Burgan (2005) 40 Standard Fire Behavior Fuel Models
    // =========================================================================
    // Reference: Scott, J.H. and Burgan, R.E. (2005). "Standard Fire Behavior
    // Fuel Models: A Comprehensive Set for Use with Rothermel's Surface Fire
    // Spread Model." Gen. Tech. Rep. RMRS-GTR-153. Fort Collins, CO: USDA
    // Forest Service, Rocky Mountain Research Station. 72 pp.
    //
    // Per-class fuel loads are from Scott & Burgan (2005) Table A-1
    // (ton/acre converted to lb/ft² via ×0.04591).
    // SAVs: 10-hr = 109 ft⁻¹, 100-hr = 30 ft⁻¹, live herb = 1800 ft⁻¹,
    //       live woody = 1500 ft⁻¹; 1-hr SAV is model-specific (σ_d1).
    // =========================================================================

    // -------------------------------------------------------------------------
    // NON-BURNABLE (NB) — no per-class data; single-class fallback (w0=0)
    // -------------------------------------------------------------------------
    database["NB1"] = {
        "NB1", "Urban/Developed",
        0.0, 0.0, 0.0, 0.30, 8000.0, 0.0555, 0.010, 32.0
    };
    database["NB2"] = {
        "NB2", "Snow/Ice",
        0.0, 0.0, 0.0, 0.30, 8000.0, 0.0555, 0.010, 32.0
    };
    database["NB3"] = {
        "NB3", "Agriculture",
        0.0, 0.0, 0.0, 0.30, 8000.0, 0.0555, 0.010, 32.0
    };
    database["NB8"] = {
        "NB8", "Open Water",
        0.0, 0.0, 0.0, 0.30, 8000.0, 0.0555, 0.010, 32.0
    };
    database["NB9"] = {
        "NB9", "Bare Ground",
        0.0, 0.0, 0.0, 0.30, 8000.0, 0.0555, 0.010, 32.0
    };

    // -------------------------------------------------------------------------
    // GRASS (GR) — σ_lh = 1800 ft⁻¹, σ_lw = 1500 ft⁻¹
    // -------------------------------------------------------------------------
    database["GR1"] = {
        "GR1", "Short, Sparse Dry Climate Grass",
        0.010, 2200.0, 0.4, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.00459, 2200.0,   // w_d1, sigma_d1
        0.0,               // w_d10
        0.0,               // w_d100
        0.01377, 1800.0,   // w_lh, sigma_lh
        0.0,     1500.0    // w_lw, sigma_lw
    };
    database["GR2"] = {
        "GR2", "Low Load, Dry Climate Grass",
        0.023, 2000.0, 1.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.00459, 2000.0,
        0.0,
        0.0,
        0.01836, 1800.0,
        0.0,     1500.0
    };
    database["GR3"] = {
        "GR3", "Low Load, Very Coarse, Humid Climate Grass",
        0.023, 1500.0, 2.0, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.00459, 1500.0,
        0.01836,
        0.0,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["GR4"] = {
        "GR4", "Moderate Load, Dry Climate Grass",
        0.087, 2000.0, 2.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.02296, 2000.0,
        0.0,
        0.0,
        0.02296, 1800.0,
        0.0,     1500.0
    };
    database["GR5"] = {
        "GR5", "Low Load, Humid Climate Grass",
        0.023, 1800.0, 1.5, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.00459, 1800.0,
        0.0,
        0.0,
        0.01836, 1800.0,
        0.0,     1500.0
    };
    database["GR6"] = {
        "GR6", "Moderate Load, Humid Climate Grass",
        0.165, 2200.0, 1.5, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.01836, 2200.0,
        0.0,
        0.0,
        0.05508, 1800.0,
        0.0,     1500.0
    };
    database["GR7"] = {
        "GR7", "High Load, Dry Climate Grass",
        0.230, 2000.0, 3.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 2000.0,
        0.0,
        0.0,
        0.05508, 1800.0,
        0.0,     1500.0
    };
    database["GR8"] = {
        "GR8", "High Load, Very Coarse, Humid Climate Grass",
        0.345, 1500.0, 4.0, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.02755, 1500.0,
        0.04591,
        0.0,
        0.05508, 1800.0,
        0.0,     1500.0
    };
    database["GR9"] = {
        "GR9", "Very High Load, Dry Climate Grass",
        0.414, 2200.0, 5.0, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 2200.0,
        0.0,
        0.0,
        0.13774, 1800.0,
        0.0,     1500.0
    };

    // -------------------------------------------------------------------------
    // GRASS-SHRUB (GS)
    // -------------------------------------------------------------------------
    database["GS1"] = {
        "GS1", "Low Load, Dry Climate Grass-Shrub",
        0.046, 2000.0, 0.9, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.00918, 2000.0,
        0.0,
        0.0,
        0.01836, 1800.0,
        0.00459, 1500.0
    };
    database["GS2"] = {
        "GS2", "Moderate Load, Dry Climate Grass-Shrub",
        0.046, 2000.0, 1.5, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.00918, 2000.0,
        0.00918,
        0.0,
        0.01836, 1800.0,
        0.01377, 1500.0
    };
    database["GS3"] = {
        "GS3", "Moderate Load, Humid Climate Grass-Shrub",
        0.115, 1800.0, 1.8, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.01377, 1800.0,
        0.01377,
        0.0,
        0.05508, 1800.0,
        0.02755, 1500.0
    };
    database["GS4"] = {
        "GS4", "High Load, Humid Climate Grass-Shrub",
        0.253, 1800.0, 2.1, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.07346, 2000.0,
        0.02755,
        0.0,
        0.06881, 1800.0,
        0.04591, 1500.0
    };

    // -------------------------------------------------------------------------
    // SHRUB (SH)
    // -------------------------------------------------------------------------
    database["SH1"] = {
        "SH1", "Low Load, Dry Climate Shrub",
        0.023, 2000.0, 1.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.00918, 2000.0,
        0.01836,
        0.0,
        0.0,     1800.0,
        0.01377, 1500.0
    };
    database["SH2"] = {
        "SH2", "Moderate Load, Dry Climate Shrub",
        0.138, 2000.0, 1.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 3000.0,
        0.04591,
        0.0,
        0.0,     1800.0,
        0.18364, 1500.0
    };
    database["SH3"] = {
        "SH3", "Moderate Load, Humid Climate Shrub",
        0.138, 1600.0, 2.4, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.02755, 1600.0,
        0.0,
        0.0,
        0.0,     1800.0,
        0.11478, 1500.0
    };
    database["SH4"] = {
        "SH4", "Low Load, Humid Climate Shrub",
        0.092, 2000.0, 3.0, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.02755, 2000.0,
        0.02755,
        0.0,
        0.00918, 1800.0,
        0.06427, 1500.0
    };
    database["SH5"] = {
        "SH5", "High Load, Dry Climate Shrub",
        0.345, 750.0,  6.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.12947, 750.0,
        0.09183,
        0.0,
        0.0,     1800.0,
        0.14327, 1500.0
    };
    database["SH6"] = {
        "SH6", "Low Load, Humid Climate Shrub",
        0.276, 750.0,  2.0, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.09183, 750.0,
        0.04591,
        0.0,
        0.0,     1800.0,
        0.12487, 1500.0
    };
    database["SH7"] = {
        "SH7", "Very High Load, Dry Climate Shrub",
        0.460, 750.0,  6.0, 0.15, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 750.0,
        0.18364,
        0.04591,
        0.0,     1800.0,
        0.18364, 1500.0
    };
    database["SH8"] = {
        "SH8", "High Load, Humid Climate Shrub",
        0.322, 1300.0, 3.0, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.09183, 750.0,
        0.18364,
        0.0,
        0.04591, 1800.0,
        0.04591, 1500.0
    };
    database["SH9"] = {
        "SH9", "Very High Load, Humid Climate Shrub",
        0.414, 1500.0, 4.4, 0.40, 8000.0, 0.0555, 0.010, 32.0,
        0.20660, 750.0,
        0.22957,
        0.0,
        0.04591, 1800.0,
        0.04591, 1500.0
    };

    // -------------------------------------------------------------------------
    // TIMBER-UNDERSTORY (TU)
    // -------------------------------------------------------------------------
    database["TU1"] = {
        "TU1", "Low Load, Dry Climate Timber-Grass-Shrub",
        0.023, 2000.0, 0.6, 0.20, 8000.0, 0.0555, 0.010, 32.0,
        0.00918, 2000.0,
        0.00459,
        0.00459,
        0.01377, 1800.0,
        0.00918, 1500.0
    };
    database["TU2"] = {
        "TU2", "Moderate Load, Humid Climate Timber-Shrub",
        0.069, 2000.0, 1.0, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 3000.0,
        0.02296,
        0.02296,
        0.00459, 1800.0,
        0.00918, 1500.0
    };
    database["TU3"] = {
        "TU3", "Moderate Load, Humid Climate Timber-Grass-Shrub",
        0.138, 1800.0, 1.3, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.01377, 1800.0,
        0.01836,
        0.01836,
        0.03673, 1800.0,
        0.03673, 1500.0
    };
    database["TU4"] = {
        "TU4", "Dwarf Conifer Understory",
        0.230, 2300.0, 0.5, 0.12, 8000.0, 0.0555, 0.010, 32.0,
        0.22957, 2300.0,
        0.0,
        0.0,
        0.0,     1800.0,
        0.00918, 1500.0
    };
    database["TU5"] = {
        "TU5", "Very High Load, Dry Climate Timber-Shrub",
        0.322, 1500.0, 1.0, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.09183, 1500.0,
        0.09183,
        0.09183,
        0.0,     1800.0,
        0.04591, 1500.0
    };

    // -------------------------------------------------------------------------
    // TIMBER LITTER (TL)
    // -------------------------------------------------------------------------
    database["TL1"] = {
        "TL1", "Low Load, Compact Conifer Litter",
        0.069, 2000.0, 0.2, 0.30, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 2000.0,
        0.02296,
        0.0,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL2"] = {
        "TL2", "Low Load, Broadleaf Litter",
        0.069, 1600.0, 0.2, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.02296, 2000.0,
        0.04591,
        0.01377,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL3"] = {
        "TL3", "Moderate Load Conifer Litter",
        0.115, 2000.0, 0.3, 0.20, 8000.0, 0.0555, 0.010, 32.0,
        0.01836, 2000.0,
        0.07346,
        0.01836,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL4"] = {
        "TL4", "Small Downed Logs",
        0.115, 2000.0, 0.4, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.02296, 2000.0,
        0.03673,
        0.03673,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL5"] = {
        "TL5", "High Load Conifer Litter",
        0.230, 2000.0, 0.6, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 2000.0,
        0.09183,
        0.09183,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL6"] = {
        "TL6", "High Load Broadleaf Litter",
        0.253, 2000.0, 0.3, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.05508, 2000.0,
        0.09183,
        0.07346,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL7"] = {
        "TL7", "Large Downed Logs",
        0.322, 2000.0, 0.4, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 2000.0,
        0.09183,
        0.09183,
        0.02296, 1800.0,
        0.02296, 1500.0
    };
    database["TL8"] = {
        "TL8", "Long, Needle Litter",
        0.184, 1800.0, 0.3, 0.35, 8000.0, 0.0555, 0.010, 32.0,
        0.09183, 1800.0,
        0.04591,
        0.02296,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["TL9"] = {
        "TL9", "Very High Load, Broadleaf Litter",
        0.322, 1600.0, 0.6, 0.35, 8000.0, 0.0555, 0.010, 32.0,
        0.09183, 1600.0,
        0.11479,
        0.06886,
        0.0,     1800.0,
        0.0,     1500.0
    };

    // -------------------------------------------------------------------------
    // SLASH-BLOWDOWN (SB)
    // -------------------------------------------------------------------------
    database["SB1"] = {
        "SB1", "Low Load Activity Fuel",
        0.069, 1500.0, 1.0, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.04591, 2000.0,
        0.18364,
        0.18364,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["SB2"] = {
        "SB2", "Moderate Load Activity Fuel or Low Load Blowdown",
        0.184, 1500.0, 1.0, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.09183, 2000.0,
        0.27546,
        0.18364,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["SB3"] = {
        "SB3", "High Load Activity Fuel or Moderate Load Blowdown",
        0.322, 1500.0, 1.2, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.22957, 2000.0,
        0.41323,
        0.18364,
        0.0,     1800.0,
        0.0,     1500.0
    };
    database["SB4"] = {
        "SB4", "High Load Blowdown",
        0.529, 1500.0, 2.7, 0.25, 8000.0, 0.0555, 0.010, 32.0,
        0.22957, 2000.0,
        0.41323,
        0.41323,
        0.22957, 1800.0,
        0.0,     1500.0
    };

    // -------------------------------------------------------------------------
    // COMMON ALIASES AND VARIATIONS
    // -------------------------------------------------------------------------
    
    // Numeric aliases (1-13) for Anderson 13 models
    for (int i = 1; i <= 13; ++i) {
        std::string fm_key = "FM" + std::to_string(i);
        if (database.count(fm_key) > 0) {
            database[std::to_string(i)] = database[fm_key];
        }
    }
    
    // Common name aliases
    database["SHORT_GRASS"] = database["FM1"];
    database["GRASS"] = database["FM1"];
    database["TIMBER_GRASS"] = database["FM2"];
    database["TALL_GRASS"] = database["FM3"];
    database["CHAPARRAL"] = database["FM4"];
    database["BRUSH"] = database["FM5"];
    database["DORMANT_BRUSH"] = database["FM6"];
    database["HARDWOOD_SLASH"] = database["FM6"];
    database["SOUTHERN_ROUGH"] = database["FM7"];
    database["TIMBER_LITTER"] = database["FM8"];
    database["CLOSED_TIMBER"] = database["FM8"];
    database["HARDWOOD_LITTER"] = database["FM9"];
    database["TIMBER_UNDERSTORY"] = database["FM10"];
    database["LIGHT_SLASH"] = database["FM11"];
    database["MEDIUM_SLASH"] = database["FM12"];
    database["HEAVY_SLASH"] = database["FM13"];
    
    return database;
}

// ============================================================================
// Lookup fuel model by name, numeric code, or alias (case-insensitive)
// ============================================================================
// Supports:
//   - Direct names: FM1, FM2, GR1, GR2, etc.
//   - Numeric codes: 1 → FM1, 2 → FM2, 101 → GR1, etc.
//   - Common aliases: SHORT_GRASS → FM1, CHAPARRAL → FM4, GRASS → FM3
// ============================================================================
inline bool lookup_fuel_model(const std::string& fuel_name, FuelModel& out_model) {
    static std::map<std::string, FuelModel> database = get_fuel_database();
    
    // Build alias map (static, computed once)
    static std::map<std::string, std::string> aliases = {
        // Anderson 13 common aliases
        {"SHORT_GRASS", "FM1"},
        {"GRASS", "FM3"},
        {"TALL_GRASS", "FM3"},
        {"CHAPARRAL", "FM4"},
        {"BRUSH", "FM5"},
        {"DORMANT_BRUSH", "FM6"},
        {"SOUTHERN_ROUGH", "FM7"},
        {"TIMBER_LITTER", "FM8"},
        {"CLOSED_TIMBER", "FM8"},
        {"HARDWOOD_LITTER", "FM9"},
        {"TIMBER_UNDERSTORY", "FM10"},
        {"LIGHT_SLASH", "FM11"},
        {"MEDIUM_SLASH", "FM12"},
        {"HEAVY_SLASH", "FM13"},
        // Scott & Burgan 40 common aliases
        {"GRASS_SHORT", "GR1"},
        {"GRASS_LOW", "GR2"},
        {"GRASS_MODERATE", "GR3"},
        {"GRASS_MODERATE_TALL", "GR4"},
        {"SHRUB_LOW", "SH1"},
        {"SHRUB_MODERATE", "SH2"}
    };
    
    // Convert to uppercase for case-insensitive lookup
    std::string upper_name = fuel_name;
    std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    
    // Try direct lookup first
    auto it = database.find(upper_name);
    if (it != database.end()) {
        out_model = it->second;
        return true;
    }
    
    // Try numeric code conversion (Anderson 13: 1-13 → FM1-FM13)
    try {
        int code = std::stoi(upper_name);
        
        // Anderson 13 codes (1-13)
        if (code >= 1 && code <= 13) {
            std::string key = "FM" + std::to_string(code);
            auto fm_it = database.find(key);
            if (fm_it != database.end()) {
                out_model = fm_it->second;
                return true;
            }
        }
        
        // Scott & Burgan 40 codes (101-204)
        if (code >= 101) {
            std::string key;
            if      (code >= 101 && code <= 109) key = "GR" + std::to_string(code - 100);
            else if (code >= 121 && code <= 124) key = "GS" + std::to_string(code - 120);
            else if (code >= 141 && code <= 149) key = "SH" + std::to_string(code - 140);
            else if (code >= 161 && code <= 165) key = "TU" + std::to_string(code - 160);
            else if (code >= 181 && code <= 189) key = "TL" + std::to_string(code - 180);
            else if (code >= 201 && code <= 204) key = "SB" + std::to_string(code - 200);
            
            if (!key.empty()) {
                auto fb_it = database.find(key);
                if (fb_it != database.end()) {
                    out_model = fb_it->second;
                    return true;
                }
            }
        }
    } catch (const std::exception&) {
        // Not a numeric code, continue to alias lookup
    }
    
    // Try alias lookup
    auto alias_it = aliases.find(upper_name);
    if (alias_it != aliases.end()) {
        auto db_it = database.find(alias_it->second);
        if (db_it != database.end()) {
            out_model = db_it->second;
            return true;
        }
    }
    
    return false;
}

// ============================================================================
// Apply fuel model to RothermelParams
// ============================================================================
// This function updates the fuel-related parameters in RothermelParams
// while preserving user-specified values for M_f (fuel moisture), slope, and
// unit conversions.  Per-class loads and SAVs are also copied so that the
// multi-class Rothermel path can be used when the model provides them.
// ============================================================================
inline void apply_fuel_model(const FuelModel& model, RothermelParams& params) {
    params.w0     = model.w0;
    params.sigma  = model.sigma;
    params.delta  = model.delta;
    params.M_x    = model.M_x;
    params.h_heat = model.h_heat;
    params.S_T    = model.S_T;
    params.S_e    = model.S_e;
    params.rho_p  = model.rho_p;
    // Per-class fuel loads and SAVs
    params.w_d1     = model.w_d1;
    params.sigma_d1 = model.sigma_d1;
    params.w_d10    = model.w_d10;
    params.w_d100   = model.w_d100;
    params.w_lh     = model.w_lh;
    params.sigma_lh = model.sigma_lh;
    params.w_lw     = model.w_lw;
    params.sigma_lw = model.sigma_lw;
    // Note: M_f (fuel moisture), per-class moistures, slope, and unit
    // conversions are NOT changed here as they are environmental/simulation-
    // specific parameters set by the user.
}

// ============================================================================
// Get list of all available fuel models
// ============================================================================
inline std::vector<std::string> get_available_fuel_models() {
    std::vector<std::string> models;
    // Anderson 13 models (FM1-FM13)
    for (int i = 1; i <= 13; ++i) {
        models.push_back("FM" + std::to_string(i));
    }
    // Scott & Burgan 40 models
    // Non-burnable
    for (const auto& m : {"NB1","NB2","NB3","NB8","NB9"})
        models.push_back(m);
    // Grass
    for (int i = 1; i <= 9; ++i)
        models.push_back("GR" + std::to_string(i));
    // Grass-Shrub
    for (int i = 1; i <= 4; ++i)
        models.push_back("GS" + std::to_string(i));
    // Shrub
    for (int i = 1; i <= 9; ++i)
        models.push_back("SH" + std::to_string(i));
    // Timber-Understory
    for (int i = 1; i <= 5; ++i)
        models.push_back("TU" + std::to_string(i));
    // Timber Litter
    for (int i = 1; i <= 9; ++i)
        models.push_back("TL" + std::to_string(i));
    // Slash-Blowdown
    for (int i = 1; i <= 4; ++i)
        models.push_back("SB" + std::to_string(i));
    return models;
}

// ============================================================================
// Print fuel model information
// ============================================================================
inline void print_fuel_model_info(const FuelModel& model) {
    amrex::Print() << "\n";
    amrex::Print() << "╔════════════════════════════════════════════════════════════════════╗\n";
    amrex::Print() << "║ Fuel Model: " << model.name << " - " << model.description << "\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║ Fuel load (w0):                " << model.w0 << " lb/ft²\n";
    amrex::Print() << "║ Surface-area-to-volume (σ):    " << model.sigma << " ft⁻¹\n";
    amrex::Print() << "║ Fuel bed depth (δ):            " << model.delta << " ft\n";
    amrex::Print() << "║ Moisture of extinction (Mx):   " << model.M_x << "\n";
    amrex::Print() << "║ Heat content (h):              " << model.h_heat << " BTU/lb\n";
    amrex::Print() << "║ Total mineral content (ST):    " << model.S_T << "\n";
    amrex::Print() << "║ Effective mineral content (Se): " << model.S_e << "\n";
    amrex::Print() << "║ Particle density (ρp):         " << model.rho_p << " lb/ft³\n";
    // Print per-class data if available
    const bool has_per_class = (model.w_d1 + model.w_d10 + model.w_d100
                               + model.w_lh + model.w_lw) > 0.0;
    if (has_per_class) {
        amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
        amrex::Print() << "║ Per-class fuel loads (multi-class Rothermel path active):          ║\n";
        if (model.w_d1   > 0.0) amrex::Print() << "║   1-hr  dead: " << model.w_d1   << " lb/ft²  σ=" << model.sigma_d1 << " ft⁻¹\n";
        if (model.w_d10  > 0.0) amrex::Print() << "║   10-hr dead: " << model.w_d10  << " lb/ft²  σ=109 ft⁻¹\n";
        if (model.w_d100 > 0.0) amrex::Print() << "║   100-hr dead:" << model.w_d100 << " lb/ft²  σ=30 ft⁻¹\n";
        if (model.w_lh   > 0.0) amrex::Print() << "║   Live herb:  " << model.w_lh   << " lb/ft²  σ=" << model.sigma_lh << " ft⁻¹\n";
        if (model.w_lw   > 0.0) amrex::Print() << "║   Live woody: " << model.w_lw   << " lb/ft²  σ=" << model.sigma_lw << " ft⁻¹\n";
    }
    amrex::Print() << "╚════════════════════════════════════════════════════════════════════╝\n";
    amrex::Print() << "\n";
}

// ============================================================================
// Print available fuel models
// ============================================================================
inline void print_available_fuel_models() {
    amrex::Print() << "\n";
    amrex::Print() << "╔════════════════════════════════════════════════════════════════════╗\n";
    amrex::Print() << "║         Anderson 13 Standard Fuel Models (Continental US)          ║\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║ GRASS AND GRASS-DOMINATED:                                         ║\n";
    amrex::Print() << "║   FM1  - Short Grass (1 ft)                                        ║\n";
    amrex::Print() << "║   FM2  - Timber (Grass and Understory)                             ║\n";
    amrex::Print() << "║   FM3  - Tall Grass (2.5 ft)                                       ║\n";
    amrex::Print() << "║ CHAPARRAL AND SHRUB:                                               ║\n";
    amrex::Print() << "║   FM4  - Chaparral (6 ft)                                          ║\n";
    amrex::Print() << "║   FM5  - Brush (2 ft)                                              ║\n";
    amrex::Print() << "║   FM6  - Dormant Brush, Hardwood Slash                             ║\n";
    amrex::Print() << "║   FM7  - Southern Rough                                            ║\n";
    amrex::Print() << "║ TIMBER LITTER:                                                     ║\n";
    amrex::Print() << "║   FM8  - Closed Timber Litter                                      ║\n";
    amrex::Print() << "║   FM9  - Hardwood Litter                                           ║\n";
    amrex::Print() << "║   FM10 - Timber (Litter and Understory)                            ║\n";
    amrex::Print() << "║ SLASH:                                                             ║\n";
    amrex::Print() << "║   FM11 - Light Logging Slash                                       ║\n";
    amrex::Print() << "║   FM12 - Medium Logging Slash                                      ║\n";
    amrex::Print() << "║   FM13 - Heavy Logging Slash                                       ║\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║       Scott & Burgan (2005) 40 Standard Fuel Models                ║\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║ NON-BURNABLE:  NB1, NB2, NB3, NB8, NB9                            ║\n";
    amrex::Print() << "║ GRASS:         GR1, GR2, GR3, GR4, GR5, GR6, GR7, GR8, GR9       ║\n";
    amrex::Print() << "║ GRASS-SHRUB:   GS1, GS2, GS3, GS4                                 ║\n";
    amrex::Print() << "║ SHRUB:         SH1, SH2, SH3, SH4, SH5, SH6, SH7, SH8, SH9       ║\n";
    amrex::Print() << "║ TIMBER-UNDERSTORY: TU1, TU2, TU3, TU4, TU5                            ║\n";
    amrex::Print() << "║ TIMBER LITTER:     TL1, TL2, TL3, TL4, TL5, TL6, TL7, TL8, TL9       ║\n";
    amrex::Print() << "║ SLASH-BLOWDOWN:    SB1, SB2, SB3, SB4                                 ║\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║ Usage: Add 'rothermel.fuel_model=FM1' to command line              ║\n";
    amrex::Print() << "║        Also accepts numeric codes: 1, 2, 3, ..., 101, 102, ...     ║\n";
    amrex::Print() << "║        Also accepts aliases: SHORT_GRASS, GRASS, CHAPARRAL, etc.   ║\n";
    amrex::Print() << "║        Examples: FM1, 1, SHORT_GRASS, GR2, 101, GRASS_SHORT        ║\n";
    amrex::Print() << "║ Landscape fuel type: set 'rothermel.landscape_fuel_type=13' or    ║\n";
    amrex::Print() << "║        '=40' to match the FBFM system in the landscape file.       ║\n";
    amrex::Print() << "╚════════════════════════════════════════════════════════════════════╝\n";
    amrex::Print() << "\n";
}

// ============================================================================
// Lookup fuel model by integer landscape code (FBFM13 or FBFM40)
// ============================================================================
// Maps the integer pixel codes stored in LANDFIRE rasters to FuelModel entries.
//
// FBFM13 (Anderson 13) integer codes:
//   1–13  → FM1–FM13 (burnable)
//   91–99 → non-burnable (returns false)
//
// FBFM40 (Scott & Burgan 40) integer codes:
//   101–109 → GR1–GR9   (Grass)
//   121–124 → GS1–GS4   (Grass-Shrub)
//   141–149 → SH1–SH9   (Shrub)
//   161–165 → TU1–TU5   (Timber-Understory)
//   181–189 → TL1–TL9   (Timber Litter)
//   201–204 → SB1–SB4   (Slash-Blowdown)
//   91–99   → non-burnable (returns false)
//
// Returns true when a burnable fuel model was found, false otherwise.
// ============================================================================
inline bool lookup_fuel_model_by_int_code(int code,
                                          const std::string& fuel_type,
                                          FuelModel& out_model)
{
    static std::map<std::string, FuelModel> database = get_fuel_database();

    std::string key;

    if (fuel_type == "40") {
        // Scott & Burgan 40 codes
        if      (code >= 101 && code <= 109) key = "GR" + std::to_string(code - 100);
        else if (code >= 121 && code <= 124) key = "GS" + std::to_string(code - 120);
        else if (code >= 141 && code <= 149) key = "SH" + std::to_string(code - 140);
        else if (code >= 161 && code <= 165) key = "TU" + std::to_string(code - 160);
        else if (code >= 181 && code <= 189) key = "TL" + std::to_string(code - 180);
        else if (code >= 201 && code <= 204) key = "SB" + std::to_string(code - 200);
        else return false; // non-burnable or unknown code
    } else {
        // Anderson 13 codes (default)
        if (code >= 1 && code <= 13) key = "FM" + std::to_string(code);
        else return false; // non-burnable or unknown code
    }

    auto it = database.find(key);
    if (it != database.end()) {
        out_model = it->second;
        return true;
    }
    return false;
}

// ============================================================================
// Check whether an integer landscape code is non-burnable
// (applies to both FBFM13 and FBFM40: codes 91–99)
// ============================================================================
inline bool is_nonburnable_fuel_code(int code)
{
    return (code == 91 || code == 92 || code == 93 ||
            code == 98 || code == 99);
}

#endif // FUEL_DATABASE_H
