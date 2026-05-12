// ============================================================================
// satellite_assimilation.H
// Real-time satellite fire detection assimilation (GOES / VIIRS / file).
// ============================================================================
//
// Supported data sources
// ----------------------
//  "goes"  – NOAA GOES-16/17/18 Active Fire products (ABI-L2-FDCF or
//             ABI-L2-FDCC) retrieved from the public AWS S3 bucket
//             s3://noaa-goes18 (or noaa-goes16 / noaa-goes17).
//             Detections are fetched via HTTP using `curl` (no library
//             dependency required).  The CSV listing from the S3 REST API
//             is parsed to find the most-recent NetCDF granule, which is
//             then converted to a lon/lat/confidence CSV locally using the
//             python3 helper `satellite_goes_to_csv.py` (optional; only
//             needed for GOES).
//
//  "viirs" – NASA FIRMS VIIRS-SNPP or VIIRS-NOAA20 active-fire data
//            fetched from the FIRMS REST API:
//              https://firms.modaps.eosdis.nasa.gov/api/area/csv/
//                  {MAP_KEY}/VIIRS_SNPP_NRT/{bbox}/{days}
//            Requires a free NASA FIRMS map key (api_key parameter).
//            The CSV response (latitude, longitude, confidence, ...) is
//            filtered by bounding box and confidence threshold, then
//            converted to simulation-domain coordinates.
//
//  "file"  – Pre-downloaded CSV (lon, lat, confidence [%]) for
//            offline / HPC environments without internet access.
//            Format:  longitude_deg,latitude_deg,confidence_pct
//            Comment lines starting with '#' are ignored.
//
// API keys
// --------
//  GOES  : No key required (public AWS S3 bucket).
//  VIIRS : Obtain a free map key at:
//            https://firms.modaps.eosdis.nasa.gov/api/map_key/
//          Then set  satellite.api_key = "YOUR_KEY"  in the inputs file.
//
// Coordinate alignment
// --------------------
// Detections arrive in WGS-84 lon/lat degrees.  A compact Karney (2011)
// / Bowring (1985) UTM forward projection converts them to UTM easting
// and northing [m] for the zone given by `satellite.utm_zone` and
// `satellite.utm_northern`.
//
// The simulation domain origin (prob_lo_x, prob_lo_y) must be aligned
// with the UTM origin so that:
//   sim_x = UTM_easting  - prob_lo_x
//   sim_y = UTM_northing - prob_lo_y
//
// For example, if your domain covers the area
//   UTM easting  [500 000, 560 000] m  → prob_lo_x = 500 000, prob_hi_x = 560 000
//   UTM northing [3 700 000, 3 760 000] m → prob_lo_y = 3 700 000, prob_hi_y = 3 760 000
// then satellite detections in that region will map directly to [0, 60 000] m
// in the simulation domain.
//
// Offline / HPC design
// --------------------
//  • satellite.source = "file" reads a local CSV with no network call.
//  • satellite.local_cache_file: when non-empty the downloaded data is
//    written to this path after each successful fetch; if a subsequent
//    network fetch fails the cache is read instead.
//  • Network errors never abort the simulation — a warning is printed and
//    the assimilation step is skipped for that interval.
//
// ParmParse keys
// --------------
//   satellite.enable                = 0       # master on/off switch
//   satellite.source                = "file"  # "goes", "viirs", or "file"
//   satellite.goes_product          = "ABI-L2-FDCF"   # GOES fire product
//   satellite.goes_bucket           = "noaa-goes18"   # AWS S3 bucket
//   satellite.viirs_url_base        = "https://firms.modaps.eosdis.nasa.gov/api/area/csv"
//   satellite.api_key               = ""      # NASA FIRMS map key
//   satellite.bbox_lon_min          = -120.0  # bounding box [deg]
//   satellite.bbox_lon_max          = -114.0
//   satellite.bbox_lat_min          =   33.0
//   satellite.bbox_lat_max          =   42.0
//   satellite.utm_zone              = 10      # UTM zone number 1-60
//   satellite.utm_northern          = 1       # 1=northern hemisphere
//   satellite.prob_lo_easting_m     = 0.0     # UTM easting  of prob_lo_x [m]
//   satellite.prob_lo_northing_m    = 0.0     # UTM northing of prob_lo_y [m]
//   satellite.fetch_interval_s      = 600.0   # re-fetch every N sim-seconds
//   satellite.use_as_ic             = 1       # apply as initial condition at t=0
//   satellite.use_mid_sim           = 1       # re-apply during simulation
//   satellite.confidence_threshold  = 50      # min detection confidence [%]
//   satellite.detection_radius_m    = 375.0   # fire disk radius per detection [m]
//   satellite.local_cache_file      = ""      # optional cache CSV path
//   satellite.suppress_if_burning   = 1       # skip already-burning cells mid-sim
//
// ============================================================================

#ifndef SATELLITE_ASSIMILATION_H
#define SATELLITE_ASSIMILATION_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

// ============================================================================
// SatelliteParams – all parmparse-exposed configuration
// ============================================================================
struct SatelliteParams {
    int         enable;                 // master on/off switch (default: 0)
    std::string source;                 // "goes", "viirs", or "file" (default: "file")
    std::string goes_product;           // GOES fire product (default: "ABI-L2-FDCF")
    std::string goes_bucket;            // AWS S3 bucket (default: "noaa-goes18")
    std::string viirs_url_base;         // FIRMS REST API base URL
    std::string api_key;                // NASA FIRMS map key
    double      bbox_lon_min;           // bounding-box min longitude [deg]
    double      bbox_lon_max;           // bounding-box max longitude [deg]
    double      bbox_lat_min;           // bounding-box min latitude  [deg]
    double      bbox_lat_max;           // bounding-box max latitude  [deg]
    int         utm_zone;               // UTM zone number 1-60 (default: 10)
    int         utm_northern;           // 1=northern hemisphere (default: 1)
    double      prob_lo_easting_m;      // UTM easting  of prob_lo_x [m] (default: 0)
    double      prob_lo_northing_m;     // UTM northing of prob_lo_y [m] (default: 0)
    amrex::Real fetch_interval_s;       // re-fetch interval [s] (default: 600)
    int         use_as_ic;              // apply as initial condition at t=0 (default: 1)
    int         use_mid_sim;            // re-apply during simulation (default: 1)
    int         confidence_threshold;   // minimum confidence [%] (default: 50)
    amrex::Real detection_radius_m;     // fire disk radius per detection [m] (default: 375)
    std::string local_cache_file;       // cache CSV path (default: "")
    std::string local_file;             // for source="file": CSV to read
    int         suppress_if_burning;    // skip already-burning cells (default: 1)
};

// ============================================================================
// SatelliteState – runtime state (tracks last fetch time)
// ============================================================================
struct SatelliteState {
    amrex::Real last_fetch_time = amrex::Real(-1.0);
};

// ============================================================================
// UTM forward projection (WGS-84 → UTM easting / northing)
// Implements the Transverse Mercator projection following the series
// approximation of Karney (2011) / IOGP (2019) EPSG-9807.
// Accuracy: < 1 mm within 3 degrees of the central meridian.
// ============================================================================
namespace satellite_detail {

inline void lonlat_to_utm(double lon_deg, double lat_deg,
                           int utm_zone, bool northern,
                           double& easting_m, double& northing_m)
{
    // WGS-84 ellipsoid constants
    constexpr double a  = 6378137.0;          // semi-major axis [m]
    constexpr double f  = 1.0 / 298.257223563;// flattening
    constexpr double k0 = 0.9996;             // scale factor at central meridian
    constexpr double E0 = 500000.0;           // false easting [m]
    constexpr double N0_S = 10000000.0;       // false northing for southern hemisphere [m]

    const double pi = std::acos(-1.0);
    const double deg2rad = pi / 180.0;

    const double lon0_deg = (utm_zone - 1) * 6.0 - 180.0 + 3.0;  // central meridian
    const double lon0 = lon0_deg * deg2rad;
    const double lon  = lon_deg  * deg2rad;
    const double lat  = lat_deg  * deg2rad;

    const double e2   = 2.0*f - f*f;                        // eccentricity squared
    const double ep2  = e2 / (1.0 - e2);                    // e'²
    const double n    = f / (2.0 - f);                      // third flattening
    const double A    = a / (1.0 + n) * (1.0 + n*n/4.0 + n*n*n*n/64.0);

    // Conformal latitude (isometric lat)
    const double e    = std::sqrt(e2);
    const double sin_lat = std::sin(lat);
    const double psi  = std::asinh(std::tan(lat)) -
                        e * std::atanh(e * sin_lat);
    const double chi  = std::atan(std::sinh(psi));

    const double dlon = lon - lon0;

    // Series coefficients (Karney 2011, up to n^3 accuracy)
    const double alpha1 = n/2.0 - 2.0*n*n/3.0 + 5.0*n*n*n/16.0;
    const double alpha2 = 13.0*n*n/48.0 - 3.0*n*n*n/5.0;
    const double alpha3 = 61.0*n*n*n/240.0;

    const double sin_chi = std::sin(chi);
    const double cos_chi = std::cos(chi);
    const double sin_dlon= std::sin(dlon);
    const double cos_dlon= std::cos(dlon);

    const double xi_p = std::atan2(sin_chi, cos_chi * cos_dlon);
    const double eta_p= std::atanh(cos_chi * sin_dlon / std::cosh(std::asinh(sin_chi)));

    // Apply Fourier-series corrections
    const double xi  = xi_p
        + alpha1 * std::sin(2.0*xi_p) * std::cosh(2.0*eta_p)
        + alpha2 * std::sin(4.0*xi_p) * std::cosh(4.0*eta_p)
        + alpha3 * std::sin(6.0*xi_p) * std::cosh(6.0*eta_p);
    const double eta = eta_p
        + alpha1 * std::cos(2.0*xi_p) * std::sinh(2.0*eta_p)
        + alpha2 * std::cos(4.0*xi_p) * std::sinh(4.0*eta_p)
        + alpha3 * std::cos(6.0*xi_p) * std::sinh(6.0*eta_p);

    easting_m  = k0 * A * eta + E0;
    northing_m = k0 * A * xi  + (northern ? 0.0 : N0_S);
}

// Run a shell command silently; return true on success (exit code 0).
inline bool run_shell_command(const std::string& cmd)
{
    int rc = std::system(cmd.c_str());
    return (rc == 0);
}

// Parse the satellite CSV from `filepath`.
// Expected columns (comma-separated, with optional header row):
//   longitude_deg, latitude_deg, confidence_pct
// Lines starting with '#' or 'l'/'L' (header like "latitude,...") are skipped.
// Returns a vector of (sim_x, sim_y) pairs in simulation-domain meters.
inline std::vector<std::pair<double,double>>
parse_detection_csv(const std::string& filepath,
                    const SatelliteParams& p)
{
    std::vector<std::pair<double,double>> pts;
    std::ifstream fin(filepath);
    if (!fin.is_open()) return pts;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        // Skip comment lines (starting with '#') and header rows.
        // A valid data line must start with a digit, '+', '-', or '.'.
        char first = line[0];
        if (first == '#') continue;
        if (!std::isdigit(static_cast<unsigned char>(first)) &&
            first != '+' && first != '-' && first != '.') continue;
        // Replace commas with spaces
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        double lon = 0.0, lat = 0.0;
        double conf = 100.0;
        if (!(iss >> lon >> lat)) continue;
        iss >> conf;   // optional; stays 100 if missing

        // Bounding-box filter
        if (lon < p.bbox_lon_min || lon > p.bbox_lon_max) continue;
        if (lat < p.bbox_lat_min || lat > p.bbox_lat_max) continue;

        // Confidence threshold
        if (conf < static_cast<double>(p.confidence_threshold)) continue;

        // Convert lon/lat → UTM
        double easting = 0.0, northing = 0.0;
        lonlat_to_utm(lon, lat, p.utm_zone, p.utm_northern == 1,
                      easting, northing);

        // Convert UTM → simulation domain coordinates
        double sim_x = easting  - p.prob_lo_easting_m;
        double sim_y = northing - p.prob_lo_northing_m;

        pts.emplace_back(sim_x, sim_y);
    }
    return pts;
}

} // namespace satellite_detail

// ============================================================================
// fetch_satellite_fire_points
//
// Main entry point: returns a vector of (sim_x, sim_y) detection points.
// Only the I/O rank does network/file I/O; results are broadcast.
// If fetch fails, returns an empty vector (simulation continues).
// ============================================================================
inline std::vector<std::pair<double,double>>
fetch_satellite_fire_points(const SatelliteParams& p,
                             amrex::Real sim_time_s)
{
    using satellite_detail::parse_detection_csv;
    using satellite_detail::run_shell_command;

    std::vector<std::pair<double,double>> pts;

    // All I/O is done on rank 0 and broadcast afterward.
    if (amrex::ParallelDescriptor::IOProcessor()) {

        std::string csv_path;

        if (p.source == "file") {
            // ---- Mode 1: pre-downloaded local file ----
            csv_path = p.local_file;
            if (csv_path.empty()) {
                amrex::Print() << "SatelliteAssimilation: source='file' but "
                                  "satellite.local_file is empty; skipping.\n";
            }
        }
        else if (p.source == "viirs") {
            // ---- Mode 2: NASA FIRMS VIIRS REST API ----
            // URL:  {base}/{key}/VIIRS_SNPP_NRT/world/1/  filtered by bbox
            // We fetch the last 1 day of detections and filter by bbox.
            if (p.api_key.empty()) {
                amrex::Print() << "SatelliteAssimilation: satellite.source='viirs' "
                                  "requires satellite.api_key; skipping.\n";
            } else {
                // Build bounding box string: lon_min,lat_min,lon_max,lat_max
                char bbox_str[128];
                std::snprintf(bbox_str, sizeof(bbox_str), "%.4f,%.4f,%.4f,%.4f",
                              p.bbox_lon_min, p.bbox_lat_min,
                              p.bbox_lon_max, p.bbox_lat_max);

                // Write to a temporary file
                csv_path = "/tmp/satellite_viirs_fetch.csv";
                char cmd[1024];
                std::snprintf(cmd, sizeof(cmd),
                    "curl -sS --max-time 30 --retry 2 "
                    "\"%s/%s/VIIRS_SNPP_NRT/%s/1\" "
                    "-o \"%s\" 2>/dev/null",
                    p.viirs_url_base.c_str(),
                    p.api_key.c_str(),
                    bbox_str,
                    csv_path.c_str());

                amrex::Print() << "SatelliteAssimilation: fetching VIIRS data "
                               << "(t=" << sim_time_s << " s)...\n";
                if (!run_shell_command(cmd)) {
                    amrex::Print() << "SatelliteAssimilation: WARNING - VIIRS fetch failed; "
                                      "using cache if available.\n";
                    csv_path = p.local_cache_file;
                } else {
                    // Cache the successful response
                    if (!p.local_cache_file.empty()) {
                        char cp_cmd[512];
                        std::snprintf(cp_cmd, sizeof(cp_cmd),
                                      "cp \"%s\" \"%s\" 2>/dev/null",
                                      csv_path.c_str(), p.local_cache_file.c_str());
                        run_shell_command(cp_cmd);
                    }
                }
            }
        }
        else if (p.source == "goes") {
            // ---- Mode 3: NOAA GOES S3 bucket ----
            // 1. List the latest GOES fire detection granule via the S3 REST API.
            // 2. Download the NetCDF file.
            // 3. Convert to lon/lat/confidence CSV using a helper script
            //    (satellite_goes_to_csv.py) if available, otherwise skip.
            //
            // Note: full GOES NetCDF ingestion requires netCDF4-python, which
            // may not be available in all environments.  The helper script path
            // can be customised; if the script is absent, we fall back to cache.

            // Build the S3 listing URL for today's GOES data
            char list_url[512];
            std::snprintf(list_url, sizeof(list_url),
                "https://s3.amazonaws.com/%s/%s/latest/",
                p.goes_bucket.c_str(), p.goes_product.c_str());

            csv_path = "/tmp/satellite_goes_fetch.csv";
            const std::string nc_path  = "/tmp/satellite_goes_latest.nc";
            const std::string list_out = "/tmp/satellite_goes_listing.xml";

            // Step 1: get S3 listing
            char cmd_list[1024];
            std::snprintf(cmd_list, sizeof(cmd_list),
                "curl -sS --max-time 30 \"%s\" -o \"%s\" 2>/dev/null",
                list_url, list_out.c_str());
            bool list_ok = run_shell_command(cmd_list);

            bool goes_ok = false;
            if (list_ok) {
                // Step 2: extract the most recent .nc key from S3 XML listing
                // (grep for <Key> tags containing the product name)
                char cmd_key[512];
                std::snprintf(cmd_key, sizeof(cmd_key),
                    "grep -oP '(?<=<Key>)[^<]+\\.nc' \"%s\" | tail -1 > "
                    "/tmp/satellite_goes_key.txt 2>/dev/null",
                    list_out.c_str());
                run_shell_command(cmd_key);

                std::ifstream fkey("/tmp/satellite_goes_key.txt");
                std::string nc_key;
                std::getline(fkey, nc_key);
                fkey.close();

                if (!nc_key.empty()) {
                    // Step 3: download the NetCDF granule
                    char cmd_dl[1024];
                    std::snprintf(cmd_dl, sizeof(cmd_dl),
                        "curl -sS --max-time 120 "
                        "\"https://s3.amazonaws.com/%s/%s\" "
                        "-o \"%s\" 2>/dev/null",
                        p.goes_bucket.c_str(), nc_key.c_str(), nc_path.c_str());
                    if (run_shell_command(cmd_dl)) {
                        // Step 4: convert NetCDF to lon/lat/confidence CSV
                        // using the companion Python helper (if present)
                        char cmd_conv[1024];
                        std::snprintf(cmd_conv, sizeof(cmd_conv),
                            "python3 satellite_goes_to_csv.py \"%s\" \"%s\" "
                            "%.4f %.4f %.4f %.4f 2>/dev/null",
                            nc_path.c_str(), csv_path.c_str(),
                            p.bbox_lon_min, p.bbox_lat_min,
                            p.bbox_lon_max, p.bbox_lat_max);
                        goes_ok = run_shell_command(cmd_conv);
                        if (!goes_ok) {
                            amrex::Print() << "SatelliteAssimilation: WARNING - GOES NetCDF "
                                             "conversion failed (satellite_goes_to_csv.py not "
                                             "found or error); falling back to cache.\n";
                        } else {
                            // Cache successful result
                            if (!p.local_cache_file.empty()) {
                                char cp_cmd[512];
                                std::snprintf(cp_cmd, sizeof(cp_cmd),
                                              "cp \"%s\" \"%s\" 2>/dev/null",
                                              csv_path.c_str(),
                                              p.local_cache_file.c_str());
                                run_shell_command(cp_cmd);
                            }
                        }
                    }
                }
            }
            if (!goes_ok) {
                // Fall back to cache
                amrex::Print() << "SatelliteAssimilation: WARNING - GOES fetch failed; "
                                  "using cache if available.\n";
                csv_path = p.local_cache_file;
            }
        }
        else {
            amrex::Print() << "SatelliteAssimilation: unknown source '"
                           << p.source << "'; skipping.\n";
            csv_path = "";
        }

        // Parse the CSV (lon, lat, confidence) → (sim_x, sim_y)
        if (!csv_path.empty()) {
            pts = parse_detection_csv(csv_path, p);
            amrex::Print() << "SatelliteAssimilation: loaded "
                           << pts.size() << " detection(s) from '"
                           << csv_path << "' (t=" << sim_time_s << " s)\n";
        }
    }

    // Broadcast point count and coordinates to all MPI ranks
    {
        int n = static_cast<int>(pts.size());
        amrex::ParallelDescriptor::Bcast(
            &n, 1, amrex::ParallelDescriptor::IOProcessorNumber());

        if (!amrex::ParallelDescriptor::IOProcessor()) {
            pts.resize(n);
        }
        if (n > 0) {
            // Pack into flat double arrays for Bcast
            std::vector<double> xs(n), ys(n);
            if (amrex::ParallelDescriptor::IOProcessor()) {
                for (int i = 0; i < n; ++i) {
                    xs[i] = pts[i].first;
                    ys[i] = pts[i].second;
                }
            } else {
                xs.resize(n); ys.resize(n);
            }
            amrex::ParallelDescriptor::Bcast(
                xs.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
            amrex::ParallelDescriptor::Bcast(
                ys.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
            if (!amrex::ParallelDescriptor::IOProcessor()) {
                for (int i = 0; i < n; ++i)
                    pts[i] = {xs[i], ys[i]};
            }
        }
    }

    return pts;
}

// ============================================================================
// apply_satellite_fire_points
//
// Marks cells as burning by merging the signed-distance function of each
// fire disk (radius = detection_radius_m) into phi:
//   phi(x) = min( phi(x), dist(x, centre) - radius )
//
// When suppress_if_burning = true, cells with phi < 0 are not touched
// (same semantics as IgnitionEvent::suppress_if_burning in ignition_schedule.H).
// ============================================================================
inline void apply_satellite_fire_points(
    MultiFab&                                       phi,
    const Geometry&                                 geom,
    const std::vector<std::pair<double,double>>&    pts,
    amrex::Real                                     detection_radius_m,
    bool                                            suppress_if_burning)
{
    if (pts.empty()) return;

    const int npts = static_cast<int>(pts.size());

    // Copy detection coordinates to device
    amrex::Gpu::DeviceVector<amrex::Real> d_px(npts), d_py(npts);
    {
        std::vector<amrex::Real> h_px(npts), h_py(npts);
        for (int i = 0; i < npts; ++i) {
            h_px[i] = static_cast<amrex::Real>(pts[i].first);
            h_py[i] = static_cast<amrex::Real>(pts[i].second);
        }
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         h_px.begin(), h_px.end(), d_px.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         h_py.begin(), h_py.end(), d_py.begin());
    }
    const amrex::Real* px  = d_px.data();
    const amrex::Real* py  = d_py.data();
    const amrex::Real  rad = detection_radius_m;
    const bool         sup = suppress_if_burning;

    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto       p  = phi.array(mfi);
        const int  np = npts;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real x = plo[0] + (i + amrex::Real(0.5)) * dx[0];
            const amrex::Real y = plo[1] + (j + amrex::Real(0.5)) * dx[1];

            // Skip already-burning cells when suppress_if_burning is set
            if (sup && p(i, j, k) < amrex::Real(0.0)) return;

            // Compute minimum distance from this cell to any detection point
            amrex::Real min_dist = amrex::Real(1.0e30);
            for (int n = 0; n < np; ++n) {
                amrex::Real dx2 = x - px[n];
                amrex::Real dy2 = y - py[n];
                amrex::Real d   = std::sqrt(dx2*dx2 + dy2*dy2);
                if (d < min_dist) min_dist = d;
            }

            // Signed distance to the fire-disk union: negative inside
            const amrex::Real sdf = min_dist - rad;
            if (sdf < p(i, j, k)) {
                p(i, j, k) = sdf;
            }
        });
    }
}

#endif // SATELLITE_ASSIMILATION_H
