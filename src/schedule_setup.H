// schedule_setup.H
// Loads and initialises all time-varying schedules and weather data used
// during the simulation: FMC, herb moisture, wind direction, WTR weather
// files, multi-station weather, retardant drops, ignition events,
// precipitation state, and the fuel-moisture conditioning spin-up.
//
// Returns a ScheduleData struct whose members replace the corresponding
// local variables in main.cpp.
#ifndef SCHEDULE_SETUP_H
#define SCHEDULE_SETUP_H

// All required headers provided by wildfire_includes.H (the umbrella header).
#include <utility>

using namespace amrex;

// All schedule / weather state loaded at startup.
struct ScheduleData {
    FMCSchedule          fmc_sched;
    HerbMoistureSchedule herb_sched;
    WindDirSchedule      wind_dir_sched;
    WtrWeatherData       wtr_data;
    bool                 wtr_active = false;
    MultiWtrWeather      multi_wtr;
    bool                 multi_wtr_active = false;
    RetardantDropList    retardant_drops;
    IgnitionSchedule     ignition_sched;
    PrecipState          precip_state;
    std::vector<std::pair<double,double>> precip_schedule; // (time_s, mm/hr)
};

// Load all schedules, initialise precipitation state, run conditioning spin-up,
// and apply t=0 values to vel / phi / inputs.
//
// Parameters:
//   inputs   – modified in-place (moisture, crown.FMC, ux/uy updated to t=0)
//   vel      – wind field updated to first schedule entry at t=0
//   phi      – phi updated when ignition_schedule has t=0 events
//   geom     – AMReX geometry
//   ba / dm  – needed to apply multi-wtr IDW wind
//   use_indicator – FARSITE indicator mode flag
inline ScheduleData setup_schedules(
    InputParameters& inputs,
    MultiFab& vel,
    MultiFab& phi,
    const Geometry& geom,
    const BoxArray& /*ba*/,
    const DistributionMapping& /*dm*/,
    bool use_indicator)
{
    ScheduleData sd;

    // ---- FMC seasonal schedule ----
    if (inputs.fmc_schedule.enable == 1) {
        if (!inputs.fmc_schedule.file.empty()) {
            load_fmc_schedule(inputs.fmc_schedule.file, sd.fmc_sched,
                              inputs.fmc_schedule.start_doy);
        } else if (inputs.fmc_schedule.use_farsite_curve == 1) {
            build_farsite_fmc_curve(sd.fmc_sched,
                                    inputs.fmc_schedule.start_doy,
                                    inputs.fmc_schedule.spring_start,
                                    inputs.fmc_schedule.summer_peak,
                                    inputs.fmc_schedule.fall_start,
                                    inputs.fmc_schedule.fall_end,
                                    static_cast<double>(inputs.fmc_schedule.fmc_min),
                                    static_cast<double>(inputs.fmc_schedule.fmc_max));
        }
        if (!sd.fmc_sched.empty()) {
            inputs.crown.FMC = static_cast<amrex::Real>(
                get_fmc_at_time(sd.fmc_sched, 0.0));
            amrex::Print() << "FMC at t=0: " << inputs.crown.FMC << " %\n";
        }
    }

    // ---- Live herbaceous moisture schedule ----
    if (inputs.herb_moisture_schedule.enable == 1) {
        if (!inputs.herb_moisture_schedule.file.empty()) {
            load_herb_moisture_schedule(inputs.herb_moisture_schedule.file,
                                        sd.herb_sched,
                                        inputs.herb_moisture_schedule.start_doy);
        } else if (inputs.herb_moisture_schedule.use_curing_curve == 1) {
            build_herb_curing_curve(sd.herb_sched,
                                    inputs.herb_moisture_schedule.start_doy,
                                    inputs.herb_moisture_schedule.spring_start,
                                    inputs.herb_moisture_schedule.summer_peak,
                                    inputs.herb_moisture_schedule.fall_start,
                                    inputs.herb_moisture_schedule.fall_end,
                                    static_cast<double>(inputs.herb_moisture_schedule.m_lh_min),
                                    static_cast<double>(inputs.herb_moisture_schedule.m_lh_max));
        }
        if (!sd.herb_sched.empty()) {
            double m_lh_pct0 = get_herb_moisture_at_time(sd.herb_sched, 0.0);
            inputs.rothermel.M_lh = static_cast<amrex::Real>(m_lh_pct0 / 100.0);
            amrex::Print() << "Herb moisture at t=0: " << m_lh_pct0
                           << " % (M_lh=" << inputs.rothermel.M_lh << " fraction)\n";
        }
    }

    // ---- Wind direction schedule ----
    if (!inputs.wind_dir_schedule_file.empty()) {
        load_wind_dir_schedule(inputs.wind_dir_schedule_file, sd.wind_dir_sched);
        auto [ux0, uy0] = get_wind_at_time(sd.wind_dir_sched, 0.0);
        inputs.ux = static_cast<amrex::Real>(ux0);
        inputs.uy = static_cast<amrex::Real>(uy0);
        init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
    }

    // ---- FARSITE .wtr hourly weather file ----
    sd.wtr_active = !inputs.wtr_file.empty();
    if (sd.wtr_active) {
        load_wtr_weather(inputs.wtr_file, sd.wtr_data,
                         inputs.wtr_start_year, inputs.wtr_start_month,
                         inputs.wtr_start_day,  inputs.wtr_start_hour);
        if (!sd.wind_dir_sched.empty())
            amrex::Print() << "WARNING: both wind_dir_schedule_file and wtr_file "
                              "set; wtr_file wind schedule takes precedence.\n";
        sd.wind_dir_sched = sd.wtr_data.wind_sched;
        auto [ux_w, uy_w] = get_wind_at_time(sd.wind_dir_sched, 0.0);
        inputs.ux = static_cast<amrex::Real>(ux_w);
        inputs.uy = static_cast<amrex::Real>(uy_w);
        init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
    }

    // ---- Aerial retardant drops ----
    if (!inputs.retardant_file.empty())
        load_retardant_drops(inputs.retardant_file, sd.retardant_drops);

    // ---- Multiple weather stations (IDW) ----
    sd.multi_wtr_active = !inputs.multi_wtr_file.empty();
    if (sd.multi_wtr_active) {
        load_multi_wtr_stations(inputs.multi_wtr_file, sd.multi_wtr,
                                inputs.wtr_start_year, inputs.wtr_start_month,
                                inputs.wtr_start_day,  inputs.wtr_start_hour,
                                static_cast<double>(inputs.multi_wtr_idw_power));
        apply_multi_wtr_to_vel(vel, geom, sd.multi_wtr, 0.0);
    }

    // ---- Multiple scheduled ignitions ----
    if (!inputs.ignition_schedule_file.empty()) {
        load_ignition_schedule(inputs.ignition_schedule_file, sd.ignition_sched);
        apply_scheduled_ignitions(phi, geom, sd.ignition_sched,
                                  Real(0.0), Real(0.0), use_indicator);
        fill_boundary_extrap(phi, geom);
    }

    // ---- Precipitation state ----
    sd.precip_state.M_d1    = static_cast<float>(inputs.rothermel.M_d1);
    sd.precip_state.M_d10   = static_cast<float>(inputs.rothermel.M_d10);
    sd.precip_state.M_d100  = static_cast<float>(inputs.rothermel.M_d100);
    sd.precip_state.M_d1000 = static_cast<float>(inputs.rothermel.M_d1000);
    sd.precip_state.initialized = true;

    // ---- Fuel moisture conditioning spin-up ----
    if (inputs.conditioning.n_days > 0 && inputs.diurnal_moisture.enable == 1) {
        amrex::Print() << "Fuel moisture conditioning: spinning up "
                       << inputs.conditioning.n_days << " day(s)...\n";
        const double cond_dt     = 3600.0;
        const int    n_steps_cond = inputs.conditioning.n_days * 24;
        WtrWeatherData cond_wtr;
        const bool cond_wtr_active =
            (!inputs.conditioning.wtr_file.empty()) ||
            (!inputs.wtr_file.empty() && sd.wtr_active);
        if (!inputs.conditioning.wtr_file.empty())
            load_wtr_weather(inputs.conditioning.wtr_file, cond_wtr,
                             inputs.wtr_start_year, inputs.wtr_start_month,
                             inputs.wtr_start_day,  inputs.wtr_start_hour);
        else if (sd.wtr_active)
            cond_wtr = sd.wtr_data;

        for (int ci = 0; ci < n_steps_cond; ++ci) {
            const double ct = static_cast<double>(ci) * cond_dt;
            float  rain_rate = static_cast<float>(inputs.precip_rain_rate_mm_hr);
            double T_cond    = 0.5*(inputs.diurnal_moisture.T_max + inputs.diurnal_moisture.T_min);
            double RH_cond   = 0.5*(inputs.diurnal_moisture.RH_max + inputs.diurnal_moisture.RH_min);
            if (cond_wtr_active && !cond_wtr.empty()) {
                auto [T_w, RH_w] = cond_wtr.get_TRH_at_time(ct);
                T_cond   = T_w;
                RH_cond  = RH_w;
                rain_rate = static_cast<float>(cond_wtr.get_precip_at_time(ct));
            }
            double rh_f = std::max(1.0, std::min(100.0, RH_cond)) / 100.0;
            double T_C  = std::max(-40.0, std::min(60.0, T_cond));
            double emc_pct;
            if      (RH_cond < 10.0) emc_pct = 0.03229 + 0.281073*rh_f - 0.000578*T_C*rh_f;
            else if (RH_cond < 50.0) emc_pct = 2.22749 + 0.160107*rh_f - 0.014784*T_C;
            else                     emc_pct = 21.0606 + 0.005565*rh_f*rh_f - 0.00035*T_C*rh_f - 0.483199*rh_f;
            emc_pct = std::max(0.5, emc_pct);
            RothermelMoistures emc_cond;
            emc_cond.M_d1   = static_cast<float>(emc_pct / 100.0);
            emc_cond.M_d10  = static_cast<float>(emc_pct * 1.10 / 100.0);
            emc_cond.M_d100 = static_cast<float>(emc_pct * 1.30 / 100.0);
            const double cond_frac = (n_steps_cond > 1)
                ? static_cast<double>(ci) / static_cast<double>(n_steps_cond - 1)
                : 1.0;
            const float M_lh_init   = static_cast<float>(inputs.rothermel.M_lh);
            const float M_lw_init   = static_cast<float>(inputs.rothermel.M_lw);
            const float M_lh_target = emc_cond.M_d100 * static_cast<float>(WildfireConst::COND_LH_EMC_MULT);
            const float M_lw_target = emc_cond.M_d100 * static_cast<float>(WildfireConst::COND_LW_EMC_MULT);
            emc_cond.M_lh = M_lh_init + static_cast<float>(cond_frac)*(M_lh_target - M_lh_init);
            emc_cond.M_lw = M_lw_init + static_cast<float>(cond_frac)*(M_lw_target - M_lw_init);
            apply_precipitation_moisture(sd.precip_state, emc_cond, rain_rate,
                                         static_cast<float>(cond_dt),
                                         static_cast<float>(inputs.precip_threshold_mm_hr),
                                         static_cast<float>(inputs.M_sat));
        }
        inputs.rothermel.M_d1    = static_cast<amrex::Real>(sd.precip_state.M_d1);
        inputs.rothermel.M_d10   = static_cast<amrex::Real>(sd.precip_state.M_d10);
        inputs.rothermel.M_d100  = static_cast<amrex::Real>(sd.precip_state.M_d100);
        inputs.rothermel.M_d1000 = static_cast<amrex::Real>(sd.precip_state.M_d1000);
        inputs.rothermel.M_f     = static_cast<amrex::Real>(sd.precip_state.M_d1);
        const float m_d100_f     = sd.precip_state.M_d100;
        inputs.rothermel.M_lh    = static_cast<amrex::Real>(m_d100_f * static_cast<float>(WildfireConst::COND_LH_EMC_MULT));
        inputs.rothermel.M_lw    = static_cast<amrex::Real>(m_d100_f * static_cast<float>(WildfireConst::COND_LW_EMC_MULT));
        amrex::Print() << "Conditioning complete: M_d1="   << sd.precip_state.M_d1
                       << " M_d10="  << sd.precip_state.M_d10
                       << " M_d100=" << sd.precip_state.M_d100
                       << " M_d1000="<< sd.precip_state.M_d1000
                       << " M_lh="   << inputs.rothermel.M_lh
                       << " M_lw="   << inputs.rothermel.M_lw << "\n";
    }

    // ---- Precipitation time series ----
    if (!inputs.precip_schedule_file.empty()) {
        std::ifstream pf(inputs.precip_schedule_file);
        if (!pf.is_open())
            amrex::Abort("Cannot open precip_schedule_file: " +
                         inputs.precip_schedule_file);
        std::string pline;
        while (std::getline(pf, pline)) {
            if (pline.empty() || pline[0] == '#' || pline[0] == '!') continue;
            for (char& c : pline) if (c == ',') c = ' ';
            std::istringstream piss(pline);
            double ts, rate;
            if (piss >> ts >> rate) sd.precip_schedule.push_back({ts, rate});
        }
        amrex::Print() << "Precipitation schedule: loaded "
                       << sd.precip_schedule.size() << " records from "
                       << inputs.precip_schedule_file << "\n";
    }

    return sd;
}

#endif // SCHEDULE_SETUP_H
