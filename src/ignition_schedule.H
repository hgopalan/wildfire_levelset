// ============================================================================
// Multiple Scheduled Ignitions
// ============================================================================
// Supports scheduling fire ignition events at arbitrary simulation times,
// enabling multi-start scenarios (e.g., lightning strikes at different times,
// spotfire confirmations, prescribed burn ignitions).
//
// File format (whitespace or comma separated; # or ! lines are comments):
//   time_s  cx  cy  radius  [source_type]
//   0       500 500 50
//   3600    1200 800 30
//   7200    300 1500 25
//
// Fields:
//   time_s      – simulation time [s] at which to ignite
//   cx          – ignition center X coordinate [m]
//   cy          – ignition center Y coordinate [m]
//   radius      – ignition sphere radius [m]
//   source_type – (optional) "sphere" (default) or "box"
//
// Usage:
//   IgnitionSchedule sched;
//   load_ignition_schedule("ignitions.csv", sched);
//   // at each timestep, call apply_scheduled_ignitions:
//   apply_scheduled_ignitions(phi, geom, sched, current_time, prev_time);
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
// ============================================================================

#ifndef IGNITION_SCHEDULE_H
#define IGNITION_SCHEDULE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

// One entry in the ignition schedule
struct IgnitionEvent {
    double      time_s;        // simulation time to ignite [s]
    double      cx;            // ignition center X [m]
    double      cy;            // ignition center Y [m]
    double      radius;        // sphere radius (or half-width for box) [m]
    std::string source_type;   // "sphere" or "box"
    bool        fired;         // true once this event has been applied
};

// The full schedule
struct IgnitionSchedule {
    std::vector<IgnitionEvent> events;
    bool empty() const { return events.empty(); }
};

// ============================================================================
// load_ignition_schedule
//
// Parse a CSV schedule file into an IgnitionSchedule.
// ============================================================================
inline void load_ignition_schedule(const std::string& filename, IgnitionSchedule& sched)
{
    sched.events.clear();
    std::ifstream f(filename);
    if (!f.is_open()) {
        amrex::Abort("IgnitionSchedule: cannot open file: " + filename);
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        IgnitionEvent ev;
        ev.source_type = "sphere";
        ev.fired       = false;
        if (iss >> ev.time_s >> ev.cx >> ev.cy >> ev.radius) {
            std::string st;
            if (iss >> st) ev.source_type = st;
            sched.events.push_back(ev);
        }
    }
    if (sched.events.empty()) {
        amrex::Abort("IgnitionSchedule: no valid records in: " + filename);
    }
    std::sort(sched.events.begin(), sched.events.end(),
              [](const IgnitionEvent& a, const IgnitionEvent& b){
                  return a.time_s < b.time_s; });
    amrex::Print() << "IgnitionSchedule: loaded " << sched.events.size()
                   << " ignition events from " << filename << "\n";
}

// ============================================================================
// apply_scheduled_ignitions
//
// For any ignition event with time in (prev_time, current_time], apply the
// ignition to phi (set phi = -radius in cells within radius of the center).
// Events are marked 'fired' to avoid re-application.
//
// Uses the same sphere initialization logic as init_phi_sphere in
// initial_conditions.H: cells whose centre lies within radius of the
// ignition point receive phi = min(phi, -(radius - dist)).
//
// Parameters:
//   phi          – level-set field (updated in place)
//   geom         – AMReX geometry
//   sched        – ignition schedule (events are marked 'fired' in place)
//   current_time – end of current timestep [s]
//   prev_time    – start of current timestep [s]
//   use_indicator – when true (FARSITE mode), write -1 to newly ignited cells
// ============================================================================
inline void apply_scheduled_ignitions(MultiFab&        phi,
                                      const Geometry&  geom,
                                      IgnitionSchedule& sched,
                                      amrex::Real       current_time,
                                      amrex::Real       prev_time,
                                      bool              use_indicator = false)
{
    if (sched.empty()) return;

    const auto dx   = geom.CellSizeArray();
    const auto plo  = geom.ProbLoArray();

    for (auto& ev : sched.events) {
        if (ev.fired) continue;
        // Fire if event time is in (prev_time, current_time]
        if (static_cast<amrex::Real>(ev.time_s) > prev_time &&
            static_cast<amrex::Real>(ev.time_s) <= current_time) {

            amrex::Print() << "IgnitionSchedule: firing event at t=" << ev.time_s
                           << " s  centre=(" << ev.cx << ", " << ev.cy
                           << ")  radius=" << ev.radius << " m"
                           << "  type=" << ev.source_type << "\n";

            const amrex::Real icx    = static_cast<amrex::Real>(ev.cx);
            const amrex::Real icy    = static_cast<amrex::Real>(ev.cy);
            const amrex::Real irad   = static_cast<amrex::Real>(ev.radius);
            const bool        use_ind = use_indicator;
            // Sphere ignition: update phi for cells within radius
            if (ev.source_type == "sphere" || ev.source_type == "box") {
                for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    auto p = phi.array(mfi);
                    const auto gdx  = dx;
                    const auto gplo = plo;
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        amrex::Real x = gplo[0] + (i + amrex::Real(0.5)) * gdx[0];
                        amrex::Real y = gplo[1] + (j + amrex::Real(0.5)) * gdx[1];
#if (AMREX_SPACEDIM == 3)
                        amrex::Real z = gplo[2] + (k + amrex::Real(0.5)) * gdx[2];
                        amrex::Real dist = std::sqrt((x-icx)*(x-icx) + (y-icy)*(y-icy) + z*z);
#else
                        amrex::Real dist = std::sqrt((x-icx)*(x-icx) + (y-icy)*(y-icy));
#endif
                        if (dist <= irad) {
                            const amrex::Real new_val = use_ind
                                ? amrex::Real(-1.0)
                                : amrex::min(p(i, j, k), -(irad - dist));
                            p(i, j, k) = amrex::min(p(i, j, k), new_val);
                        }
                    });
                }
            }
            ev.fired = true;
        }
    }
}

#endif // IGNITION_SCHEDULE_H
