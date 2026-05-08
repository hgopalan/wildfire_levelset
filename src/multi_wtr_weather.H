// ============================================================================
// Multiple Weather Stations with Spatial Interpolation
// ============================================================================
// Extends the single-station .wtr weather reader to support multiple weather
// stations, each with an (X, Y) position.  Per-cell T, RH, wind speed, and
// wind direction are produced via inverse-distance-weighting (IDW) of the
// time-interpolated station values at each timestep.
//
// Station file format (multi_wtr_stations.csv):
//   # station_id, x_m, y_m, wtr_file
//   1, 330000.0, 3775000.0, station1.wtr
//   2, 335000.0, 3775000.0, station2.wtr
//   3, 332500.0, 3780000.0, station3.wtr
//
// Each wtr_file is a standard FARSITE .wtr file as parsed by load_wtr_weather().
//
// Spatial interpolation:
//   For each grid cell, the station-averaged value is computed as:
//     V_cell = Σ_i (w_i × V_i) / Σ_i w_i
//   where
//     w_i = 1 / max(d_i^p, eps)
//     d_i = distance from cell centre to station i [m]
//     p   = IDW power (default: 2.0)
//     eps = floor distance [m] (default: 1.0) to avoid division-by-zero
//           when a cell coincides exactly with a station location.
//
// Wind vector interpolation:
//   Wind U and V components are interpolated separately (IDW in Cartesian
//   components), then converted back to speed and direction.
//
// Usage:
//   MultiWtrWeather mww;
//   load_multi_wtr_stations("stations.csv", mww,
//                           start_year, start_month, start_day, start_hour);
//
//   // At each timestep:
//   apply_multi_wtr_to_vel(vel_mf, geom, mww, sim_time_s);
//   auto [T_mean, RH_mean] = mww.get_domain_TRH_at_time(sim_time_s);
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   Shepard, D. (1968). A two-dimensional interpolation function for
//     irregularly-spaced data. Proc. 23rd ACM Conference, pp. 517–524.
// ============================================================================

#ifndef MULTI_WTR_WEATHER_H
#define MULTI_WTR_WEATHER_H

#include "wtr_weather.H"
#include "wind_dir_schedule.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Gpu.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

// ============================================================================
// Single weather-station descriptor
// ============================================================================
struct WtrStation {
    double x_m;             // station easting  [m]
    double y_m;             // station northing [m]
    int    id;              // station identifier (informational)
    WtrWeatherData data;    // parsed .wtr data for this station
};

// ============================================================================
// Multi-station container
// ============================================================================
struct MultiWtrWeather {
    std::vector<WtrStation> stations;  // all loaded stations
    double idw_power;                  // IDW exponent p (default: 2.0)
    double idw_eps;                    // distance floor [m] (default: 1.0)

    bool empty() const { return stations.empty(); }

    // Domain-mean (IDW centroid) T_C and RH_pct at time t_s.
    // Useful as a scalar fallback when spatial interpolation is not needed
    // (e.g., for the global diurnal moisture scalar in non-spatial mode).
    std::pair<double,double> get_domain_TRH_at_time(double t_s) const {
        if (stations.empty()) return {25.0, 30.0};
        double T_sum = 0.0, RH_sum = 0.0;
        for (const auto& st : stations) {
            auto [T, RH] = st.data.get_TRH_at_time(t_s);
            T_sum  += T;
            RH_sum += RH;
        }
        const double n = static_cast<double>(stations.size());
        return {T_sum / n, RH_sum / n};
    }

    // Domain-mean precipitation rate [mm/hr] at time t_s.
    double get_domain_precip_at_time(double t_s) const {
        if (stations.empty()) return 0.0;
        double sum = 0.0;
        for (const auto& st : stations) sum += st.data.get_precip_at_time(t_s);
        return sum / static_cast<double>(stations.size());
    }
};

// ============================================================================
// parse_station_list_file
//
// Reads a CSV file with columns:  id, x_m, y_m, wtr_file
// Returns a vector of (id, x_m, y_m, wtr_file_path) tuples for loading.
// ============================================================================
struct StationSpec { int id; double x_m, y_m; std::string wtr_path; };

inline std::vector<StationSpec>
parse_station_list_file(const std::string& list_file)
{
    std::vector<StationSpec> specs;
    std::ifstream ifs(list_file);
    if (!ifs.is_open()) {
        amrex::Abort("multi_wtr_weather: cannot open station list: " + list_file);
    }
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        // Replace commas with spaces for easy parsing
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        StationSpec s;
        if (!(iss >> s.id >> s.x_m >> s.y_m >> s.wtr_path)) continue;
        specs.push_back(s);
    }
    return specs;
}

// ============================================================================
// load_multi_wtr_stations
//
// Parses the station list file and loads each station's .wtr data.
// The .wtr paths are resolved relative to the directory of list_file.
//
// Parameters:
//   list_file    – path to station list CSV
//   mww          – output: populated MultiWtrWeather
//   start_year   – calendar year  of t=0 (0 = auto)
//   start_month  – calendar month of t=0 (0 = auto)
//   start_day    – calendar day   of t=0 (0 = auto)
//   start_hour   – HHMM hour      of t=0 (0 = auto)
//   idw_power    – IDW exponent (default: 2.0)
// ============================================================================
inline void load_multi_wtr_stations(const std::string& list_file,
                                     MultiWtrWeather&   mww,
                                     int  start_year  = 0,
                                     int  start_month = 0,
                                     int  start_day   = 0,
                                     int  start_hour  = 0,
                                     double idw_power = 2.0)
{
    mww.idw_power = idw_power;
    mww.idw_eps   = 1.0;  // 1 m floor

    // Resolve wtr paths relative to the directory containing list_file
    std::string dir;
    auto slash = list_file.find_last_of("/\\");
    if (slash != std::string::npos) dir = list_file.substr(0, slash + 1);

    const auto specs = parse_station_list_file(list_file);
    mww.stations.reserve(specs.size());

    for (const auto& sp : specs) {
        WtrStation st;
        st.id  = sp.id;
        st.x_m = sp.x_m;
        st.y_m = sp.y_m;
        const std::string wtr_path = dir.empty() ? sp.wtr_path : dir + sp.wtr_path;
        load_wtr_weather(wtr_path, st.data,
                         start_year, start_month, start_day, start_hour);
        mww.stations.push_back(std::move(st));
        amrex::Print() << "  Loaded weather station " << sp.id
                       << " at (" << sp.x_m << ", " << sp.y_m << ") m"
                       << " from " << wtr_path << "\n";
    }
    amrex::Print() << "Multi-station weather: " << mww.stations.size()
                   << " stations loaded (IDW power=" << mww.idw_power << ")\n";
}

// ============================================================================
// apply_multi_wtr_to_vel
//
// Updates the velocity (wind) MultiFab using IDW-interpolated station winds.
// Each station provides a (speed, dir) pair at time t_s from its WindDirSchedule.
// The U and V components are IDW-interpolated over the domain, then stored in
// the vel MultiFab.
//
// Parameters:
//   vel_mf       – velocity field (3 components; comps 0 and 1 updated)
//   geom         – AMReX geometry
//   mww          – multi-station container
//   t_s          – current simulation time [s]
// ============================================================================
inline void apply_multi_wtr_to_vel(amrex::MultiFab& vel_mf,
                                    const amrex::Geometry& geom,
                                    const MultiWtrWeather& mww,
                                    double t_s)
{
    if (mww.empty()) return;

    const int nst = static_cast<int>(mww.stations.size());
    const double idw_p = mww.idw_power;
    const double idw_e = mww.idw_eps;

    // Collect per-station (x, y, ux, uy) on host
    std::vector<double> sx(nst), sy(nst), su(nst), sv(nst);
    for (int k = 0; k < nst; ++k) {
        sx[k] = mww.stations[k].x_m;
        sy[k] = mww.stations[k].y_m;
        auto [ux_st, uy_st] = get_wind_at_time(mww.stations[k].data.wind_sched, t_s);
        su[k] = ux_st;
        sv[k] = uy_st;
    }

    // Copy station data to device
    amrex::Gpu::DeviceVector<amrex::Real> d_sx(nst), d_sy(nst), d_su(nst), d_sv(nst);
    for (int k = 0; k < nst; ++k) {
        d_sx[k] = static_cast<amrex::Real>(sx[k]);
        d_sy[k] = static_cast<amrex::Real>(sy[k]);
        d_su[k] = static_cast<amrex::Real>(su[k]);
        d_sv[k] = static_cast<amrex::Real>(sv[k]);
    }
    const amrex::Real* psx = d_sx.data();
    const amrex::Real* psy = d_sy.data();
    const amrex::Real* psu = d_su.data();
    const amrex::Real* psv = d_sv.data();
    const int nst_g = nst;
    const amrex::Real pidw = static_cast<amrex::Real>(idw_p);
    const amrex::Real peps = static_cast<amrex::Real>(idw_e);

    const auto dx   = geom.CellSizeArray();
    const auto plo  = geom.ProbLoArray();

    for (amrex::MFIter mfi(vel_mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto v = vel_mf.array(mfi);
        const auto gdx  = dx;
        const auto gplo = plo;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real cx = gplo[0] + (i + amrex::Real(0.5)) * gdx[0];
            const amrex::Real cy = gplo[1] + (j + amrex::Real(0.5)) * gdx[1];

            amrex::Real w_sum = amrex::Real(0.0);
            amrex::Real ux_sum = amrex::Real(0.0);
            amrex::Real uy_sum = amrex::Real(0.0);
            for (int s = 0; s < nst_g; ++s) {
                const amrex::Real dx_ = cx - psx[s];
                const amrex::Real dy_ = cy - psy[s];
                amrex::Real dist = std::sqrt(dx_*dx_ + dy_*dy_);
                if (dist < peps) dist = peps;  // avoid divide-by-zero
                const amrex::Real w = amrex::Real(1.0) / std::pow(dist, pidw);
                w_sum  += w;
                ux_sum += w * psu[s];
                uy_sum += w * psv[s];
            }
            const amrex::Real inv_w = (w_sum > amrex::Real(0.0))
                                      ? amrex::Real(1.0) / w_sum
                                      : amrex::Real(0.0);
            v(i, j, k, 0) = ux_sum * inv_w;
            v(i, j, k, 1) = uy_sum * inv_w;
        });
    }
}

// ============================================================================
// parse_inputs integration helper
//
// Add to InputParameters in parse_inputs.H:
//   std::string multi_wtr_file;  // path to station list CSV (default: "")
// ============================================================================

#endif // MULTI_WTR_WEATHER_H
