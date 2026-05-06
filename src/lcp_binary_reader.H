// ============================================================================
// Binary FARSITE Landscape File (.lcp) Reader
// ============================================================================
// Reads the standard FARSITE binary LCP format as documented in:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator — Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
//
// Binary LCP structure:
//   Header  : 7316 bytes (little-endian)
//   Data    : NumRows × NumCols × NumLayers × sizeof(int16) bytes
//
// Always-present layers (5): elevation, slope, aspect, fuel model, canopy cover
// Crown fuel layers   (3): stand height, crown base height (CBH),
//                          crown bulk density (CBD)  — when CrownFuels != 0
// Ground fuel layers  (2): duff, coarse woody debris               — optional
//
// Layer data is stored top-to-bottom (northernmost row first), left-to-right
// (westernmost column first), as signed 16-bit integers.
//
// Usage:
//   LcpGridInfo info;
//   std::vector<float> elev, slope, aspect, cbh_m, cbd_kgm3, cc_pct, ht_m;
//   std::vector<int>   fuel;
//   bool ok = read_lcp_binary("fire.lcp", info, elev, slope, aspect, fuel,
//                              cc_pct, ht_m, cbh_m, cbd_kgm3);
//   if (!ok) { /* file not binary or unreadable – fall back to ASCII */ }
// ============================================================================

#ifndef LCP_BINARY_READER_H
#define LCP_BINARY_READER_H

#include <cstdint>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// ============================================================================
// Grid/geometry metadata extracted from the LCP header
// ============================================================================
struct LcpGridInfo {
    int32_t  numrows;       // number of rows    (north-south; row 0 = north)
    int32_t  numcols;       // number of columns (east-west;  col 0 = west)
    double   xllcorner;     // UTM easting  of lower-left corner [m or ft]
    double   yllcorner;     // UTM northing of lower-left corner [m or ft]
    double   cellsize;      // cell dimension [m or ft]
    int32_t  gridunits;     // 0 = metres, 1 = feet

    bool     has_crown;     // true when crown fuel layers are present
    bool     has_ground;    // true when ground fuel layers are present

    // Per-layer unit codes
    int16_t  EUnits;        // elevation:    0=metres, 1=feet
    int16_t  SUnits;        // slope:        0=degrees, 1=percent
    int16_t  AUnits;        // aspect:       0=degrees from north
    int16_t  CUnits;        // canopy cover: 0=percent
    int16_t  HUnits;        // stand height: 0=m*10, 1=ft, 2=m, 3=ft*10
    int16_t  BUnits;        // crown base h: 0=m*10, 1=ft, 2=m, 3=ft*10
    int16_t  PUnits;        // crown density:0=100*kg/m³, 1=10*t/ac, 2=kg/m³, 3=t/ac
};

// ============================================================================
// Internal helpers: unit conversion
// ============================================================================
namespace lcp_detail {

static inline double m_per_ft() { return 0.30480; }

// Elevation: raw int16 → metres
static inline float elev_to_m(int16_t raw, int16_t unit) {
    return (unit == 1) ? static_cast<float>(raw) * 0.30480f   // feet → m
                       : static_cast<float>(raw);              // metres
}

// Slope: raw int16 → degrees
static inline float slope_to_deg(int16_t raw, int16_t unit) {
    if (unit == 1) {
        // percent → degrees: atan(pct/100)
        return static_cast<float>(
            std::atan(static_cast<double>(raw) / 100.0) * (180.0 / M_PI));
    }
    return static_cast<float>(raw); // already degrees
}

// Aspect: raw int16 → degrees (compass from North, 0-360)
static inline float aspect_to_deg(int16_t raw) {
    // FARSITE stores aspect as degrees from North, 0-360.
    // No-data / flat is typically -1 or 0; clamp to [0,360].
    if (raw < 0) return 0.0f;
    return static_cast<float>(raw % 360);
}

// Height (stand height, CBH) → metres
static inline float height_to_m(int16_t raw, int16_t unit) {
    if (raw <= 0) return 0.0f;
    switch (unit) {
        case 0: return static_cast<float>(raw) * 0.1f;          // m*10 → m
        case 1: return static_cast<float>(raw) * 0.30480f;      // ft   → m
        case 2: return static_cast<float>(raw);                  // m
        case 3: return static_cast<float>(raw) * 0.030480f;     // ft*10→ m
        default: return static_cast<float>(raw);
    }
}

// Crown bulk density → kg/m³
static inline float cbd_to_kgm3(int16_t raw, int16_t unit) {
    if (raw <= 0) return 0.0f;
    switch (unit) {
        case 0: return static_cast<float>(raw) / 100.0f;        // 100*kg/m³ → kg/m³
        case 1: return static_cast<float>(raw) * 0.1f * 0.2242f;// 10*t/ac   → kg/m³ (~0.02242)
        case 2: return static_cast<float>(raw);                  // kg/m³
        case 3: return static_cast<float>(raw) * 0.2242f;       // t/ac      → kg/m³
        default: return static_cast<float>(raw);
    }
}

// Read a little-endian int32 from a byte buffer
static inline int32_t read_i32(const uint8_t* p) {
    int32_t v;
    std::memcpy(&v, p, 4);
    return v;
}

// Read a little-endian double from a byte buffer
static inline double read_f64(const uint8_t* p) {
    double v;
    std::memcpy(&v, p, 8);
    return v;
}

// Read a little-endian int16 from a byte buffer
static inline int16_t read_i16(const uint8_t* p) {
    int16_t v;
    std::memcpy(&v, p, 2);
    return v;
}

} // namespace lcp_detail

// ============================================================================
// is_binary_lcp
//
// Returns true if the file looks like a binary FARSITE LCP by checking:
//  1. The file is large enough to contain the 7316-byte header.
//  2. The first int32 (CrownFuels) is 0 or 20.
//  3. The second int32 (GroundFuels) is 0 or 20.
// ============================================================================
inline bool is_binary_lcp(const std::string& filename) {
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open()) return false;

    // Need at least the 7316-byte header
    f.seekg(0, std::ios::end);
    std::streamsize sz = f.tellg();
    if (sz < 7316) return false;

    f.seekg(0, std::ios::beg);
    uint8_t buf[8];
    f.read(reinterpret_cast<char*>(buf), 8);
    if (!f) return false;

    int32_t cf = lcp_detail::read_i32(buf + 0);
    int32_t gf = lcp_detail::read_i32(buf + 4);
    // Valid values: 0 (not present) or 20 (present)
    return (cf == 0 || cf == 20) && (gf == 0 || gf == 20);
}

// ============================================================================
// read_lcp_binary
//
// Reads a binary FARSITE LCP file.
//
// Output vectors are sized numrows × numcols (row-major, row 0 = north).
// All physical quantities are returned in SI:
//   elevation   [m], slope [deg], aspect [deg from north],
//   canopy_cover [%], stand_height [m], cbh [m], cbd [kg/m³].
//   fuel_model  is an integer code (FBFM13 or FBFM40).
//
// Crown layer vectors (canopy_cover_pct, stand_height_m, cbh_m, cbd_kgm3)
// are resized to numrows×numcols only when grid_info.has_crown == true.
// They are cleared otherwise.
//
// Returns true on success; false if the file is not a valid binary LCP.
// ============================================================================
inline bool read_lcp_binary(
    const std::string&  filename,
    LcpGridInfo&        grid_info,
    std::vector<float>& elevation_m,
    std::vector<float>& slope_deg,
    std::vector<float>& aspect_deg,
    std::vector<int>&   fuel_model,
    std::vector<float>& canopy_cover_pct,   // percent  [0,100]
    std::vector<float>& stand_height_m,     // metres
    std::vector<float>& cbh_m,              // metres
    std::vector<float>& cbd_kgm3)           // kg/m³
{
    using namespace lcp_detail;

    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open()) return false;

    // ---- Read the full 7316-byte header ----
    constexpr int HEADER_SIZE = 7316;
    std::vector<uint8_t> hdr(HEADER_SIZE);
    f.read(reinterpret_cast<char*>(hdr.data()), HEADER_SIZE);
    if (!f) return false;

    // Flags
    int32_t cr_flag = read_i32(hdr.data() + 0);
    int32_t gf_flag = read_i32(hdr.data() + 4);
    if ((cr_flag != 0 && cr_flag != 20) ||
        (gf_flag != 0 && gf_flag != 20)) return false;

    grid_info.has_crown  = (cr_flag == 20);
    grid_info.has_ground = (gf_flag == 20);

    // Grid dimensions and extent
    grid_info.numrows   = read_i32(hdr.data() + 4164);
    grid_info.numcols   = read_i32(hdr.data() + 4168);
    grid_info.xllcorner = read_f64(hdr.data() + 4172);
    grid_info.yllcorner = read_f64(hdr.data() + 4180);
    grid_info.cellsize  = read_f64(hdr.data() + 4188);
    grid_info.gridunits = read_i32(hdr.data() + 4196);

    if (grid_info.numrows <= 0 || grid_info.numcols <= 0) return false;

    // Unit codes (offset 4208+)
    grid_info.EUnits = read_i16(hdr.data() + 4208);
    grid_info.SUnits = read_i16(hdr.data() + 4210);
    grid_info.AUnits = read_i16(hdr.data() + 4212);
    // FOptions at 4214 (0=FBFM13, 1=custom) – informational only
    grid_info.CUnits = read_i16(hdr.data() + 4216);
    grid_info.HUnits = read_i16(hdr.data() + 4218);
    grid_info.BUnits = read_i16(hdr.data() + 4220);
    grid_info.PUnits = read_i16(hdr.data() + 4222);

    const int ncells = grid_info.numrows * grid_info.numcols;

    // Determine number of layers to read
    int nlayers = 5; // elevation, slope, aspect, fuel, canopy cover
    if (grid_info.has_crown)  nlayers += 3; // stand height, CBH, CBD
    if (grid_info.has_ground) nlayers += 2; // duff, coarse woody

    // Read all layer data
    std::vector<std::vector<int16_t>> layers(nlayers,
                                             std::vector<int16_t>(ncells));
    for (int l = 0; l < nlayers; ++l) {
        f.read(reinterpret_cast<char*>(layers[l].data()),
               static_cast<std::streamsize>(ncells) * sizeof(int16_t));
        if (!f) return false;
    }

    // ---- Convert layers to SI physical units ----
    elevation_m.resize(ncells);
    slope_deg.resize(ncells);
    aspect_deg.resize(ncells);
    fuel_model.resize(ncells);
    canopy_cover_pct.resize(ncells);

    const auto& raw_elev   = layers[0];
    const auto& raw_slope  = layers[1];
    const auto& raw_aspect = layers[2];
    const auto& raw_fuel   = layers[3];
    const auto& raw_cover  = layers[4];

    for (int idx = 0; idx < ncells; ++idx) {
        elevation_m[idx]      = elev_to_m(raw_elev[idx],   grid_info.EUnits);
        slope_deg[idx]        = slope_to_deg(raw_slope[idx], grid_info.SUnits);
        aspect_deg[idx]       = aspect_to_deg(raw_aspect[idx]);
        fuel_model[idx]       = static_cast<int>(raw_fuel[idx]);
        canopy_cover_pct[idx] = static_cast<float>(raw_cover[idx] > 0 ?
                                                    raw_cover[idx] : 0);
    }

    // ---- Crown layers (stand height, CBH, CBD) ----
    stand_height_m.clear();
    cbh_m.clear();
    cbd_kgm3.clear();
    if (grid_info.has_crown) {
        // Layer indices: 5 = stand height, 6 = CBH, 7 = CBD
        const auto& raw_ht   = layers[5];
        const auto& raw_cbh  = layers[6];
        const auto& raw_cbd  = layers[7];
        stand_height_m.resize(ncells);
        cbh_m.resize(ncells);
        cbd_kgm3.resize(ncells);
        for (int idx = 0; idx < ncells; ++idx) {
            stand_height_m[idx] = height_to_m(raw_ht[idx],  grid_info.HUnits);
            cbh_m[idx]          = height_to_m(raw_cbh[idx], grid_info.BUnits);
            cbd_kgm3[idx]       = cbd_to_kgm3(raw_cbd[idx], grid_info.PUnits);
        }
    }

    return true;
}

#endif // LCP_BINARY_READER_H
