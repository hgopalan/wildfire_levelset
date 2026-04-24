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
};

// ============================================================================
// Fuel Database - Anderson 13 Models
// ============================================================================
inline std::map<std::string, FuelModel> get_fuel_database() {
    std::map<std::string, FuelModel> database;
    
    // -------------------------------------------------------------------------
    // GRASS AND GRASS-DOMINATED MODELS (1-3)
    // -------------------------------------------------------------------------
    
    database["FM1"] = {
        "FM1", "Short Grass (1 ft)",
        0.034,    // w0 [lb/ft²]
        3500.0,   // sigma [ft⁻¹]
        1.0,      // delta [ft]
        0.12,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
    };
    
    // -------------------------------------------------------------------------
    // CHAPARRAL AND SHRUB MODELS (4-7)
    // -------------------------------------------------------------------------
    
    database["FM4"] = {
        "FM4", "Chaparral (6 ft)",
        0.230,    // w0 [lb/ft²]
        2000.0,   // sigma [ft⁻¹]
        6.0,      // delta [ft]
        0.20,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0      // rho_p [lb/ft³]
    };
    
    database["FM5"] = {
        "FM5", "Brush (2 ft)",
        0.046,    // w0 [lb/ft²]
        2000.0,   // sigma [ft⁻¹]
        2.0,      // delta [ft]
        0.20,     // M_x
        8000.0,   // h_heat [BTU/lb]
        0.0555,   // S_T
        0.010,    // S_e
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
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
        32.0      // rho_p [lb/ft³]
    };
    
    // -------------------------------------------------------------------------
    // COMMON ALIASES AND VARIATIONS
    // -------------------------------------------------------------------------
    
    // Numeric aliases (1-13)
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
// Lookup fuel model by name (case-insensitive)
// ============================================================================
inline bool lookup_fuel_model(const std::string& fuel_name, FuelModel& out_model) {
    static std::map<std::string, FuelModel> database = get_fuel_database();
    
    // Convert to uppercase for case-insensitive lookup
    std::string upper_name = fuel_name;
    std::transform(upper_name.begin(), upper_name.end(), upper_name.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    
    auto it = database.find(upper_name);
    if (it != database.end()) {
        out_model = it->second;
        return true;
    }
    
    return false;
}

// ============================================================================
// Apply fuel model to RothermelParams
// ============================================================================
// This function updates the fuel-related parameters in RothermelParams
// while preserving user-specified values for M_f (fuel moisture), slope, and
// unit conversions.
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
    // Note: M_f (fuel moisture), slope, and unit conversions are NOT changed
    // as they are typically environmental/simulation-specific parameters
}

// ============================================================================
// Get list of all available fuel models
// ============================================================================
inline std::vector<std::string> get_available_fuel_models() {
    std::vector<std::string> models;
    // Return only the primary FM1-FM13 models
    for (int i = 1; i <= 13; ++i) {
        models.push_back("FM" + std::to_string(i));
    }
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
    amrex::Print() << "╚════════════════════════════════════════════════════════════════════╝\n";
    amrex::Print() << "\n";
}

// ============================================================================
// Print available fuel models
// ============================================================================
inline void print_available_fuel_models() {
    auto models = get_available_fuel_models();
    auto database = get_fuel_database();
    
    amrex::Print() << "\n";
    amrex::Print() << "╔════════════════════════════════════════════════════════════════════╗\n";
    amrex::Print() << "║            Anderson 13 Standard Fuel Models (Continental US)       ║\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║                                                                    ║\n";
    amrex::Print() << "║ GRASS AND GRASS-DOMINATED:                                         ║\n";
    amrex::Print() << "║   FM1  - Short Grass (1 ft)                                        ║\n";
    amrex::Print() << "║   FM2  - Timber (Grass and Understory)                             ║\n";
    amrex::Print() << "║   FM3  - Tall Grass (2.5 ft)                                       ║\n";
    amrex::Print() << "║                                                                    ║\n";
    amrex::Print() << "║ CHAPARRAL AND SHRUB:                                               ║\n";
    amrex::Print() << "║   FM4  - Chaparral (6 ft)                                          ║\n";
    amrex::Print() << "║   FM5  - Brush (2 ft)                                              ║\n";
    amrex::Print() << "║   FM6  - Dormant Brush, Hardwood Slash                             ║\n";
    amrex::Print() << "║   FM7  - Southern Rough                                            ║\n";
    amrex::Print() << "║                                                                    ║\n";
    amrex::Print() << "║ TIMBER LITTER:                                                     ║\n";
    amrex::Print() << "║   FM8  - Closed Timber Litter                                      ║\n";
    amrex::Print() << "║   FM9  - Hardwood Litter                                           ║\n";
    amrex::Print() << "║   FM10 - Timber (Litter and Understory)                            ║\n";
    amrex::Print() << "║                                                                    ║\n";
    amrex::Print() << "║ SLASH:                                                             ║\n";
    amrex::Print() << "║   FM11 - Light Logging Slash                                       ║\n";
    amrex::Print() << "║   FM12 - Medium Logging Slash                                      ║\n";
    amrex::Print() << "║   FM13 - Heavy Logging Slash                                       ║\n";
    amrex::Print() << "║                                                                    ║\n";
    amrex::Print() << "╠════════════════════════════════════════════════════════════════════╣\n";
    amrex::Print() << "║ Usage: Add 'rothermel.fuel_model=FM1' to command line             ║\n";
    amrex::Print() << "║        Also accepts: 1, SHORT_GRASS, GRASS, CHAPARRAL, etc.       ║\n";
    amrex::Print() << "╚════════════════════════════════════════════════════════════════════╝\n";
    amrex::Print() << "\n";
}

#endif // FUEL_DATABASE_H
