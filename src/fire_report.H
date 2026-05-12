// fire_report.H
// Self-contained HTML report generator for a wildfire simulation run.
//
// write_fire_report_html() reads the fire_stats CSV file written by
// append_fire_stats() / write_fire_stats_header() and produces a single
// HTML file that contains:
//   - Simulation parameter summary
//   - Final fire statistics
//   - Per-step time-series table
//   - SVG line charts (burned area, perimeter, head ROS, cumulative emissions)
//
// The generated file has no external dependencies (no CDN, no images) and
// can be opened in any modern browser without an internet connection.
//
// Call once at end-of-run, after all plotfile writes are complete:
//   write_fire_report_html(inputs, final_time_s, final_step);
//
// Requires:
//   inputs.fire_report_file – path for the HTML output (default: "" = disabled)
//   inputs.fire_stats_file  – path to the fire_stats CSV to embed
//
// Only the IO processor writes the file; the call is safe on all MPI ranks.

#ifndef FIRE_REPORT_H
#define FIRE_REPORT_H

#include "parse_inputs.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

namespace {

// ---- Internal CSV row parsed from fire_stats_file ----
struct FireStatsRow {
    int    step;
    double time_s;
    double burned_ha;
    double perimeter_km;
    int    active_front;
    double head_ros_ms;
    double co2_kg_m2;
    double co_kg_m2;
    double pm25_kg_m2;
};

// Read the fire_stats CSV; returns empty vector on failure.
inline std::vector<FireStatsRow> read_fire_stats_csv(const std::string& path)
{
    std::vector<FireStatsRow> rows;
    if (path.empty()) return rows;
    std::ifstream f(path);
    if (!f.is_open()) return rows;

    std::string line;
    // Skip header row
    if (!std::getline(f, line)) return rows;

    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        FireStatsRow r{};
        if (iss >> r.step >> r.time_s >> r.burned_ha >> r.perimeter_km
                >> r.active_front >> r.head_ros_ms
                >> r.co2_kg_m2 >> r.co_kg_m2 >> r.pm25_kg_m2) {
            rows.push_back(r);
        }
    }
    return rows;
}

// ---- Minimal SVG polyline chart ----
// Returns an <svg> element as a string.  x_vals / y_vals are parallel arrays.
// The chart has a fixed 600×200 px canvas with left/bottom margins for labels.
inline std::string make_svg_chart(const std::vector<double>& x_vals,
                                  const std::vector<double>& y_vals,
                                  const std::string& x_label,
                                  const std::string& y_label,
                                  const std::string& color)
{
    if (x_vals.empty()) return "<p>(no data)</p>";

    const int W = 600, H = 220;
    const int mx = 60, my = 15, bx = 40, by = 30; // margins: left, top, right, bottom
    const int pw = W - mx - bx; // plot width
    const int ph = H - my - by; // plot height

    double x_min = *std::min_element(x_vals.begin(), x_vals.end());
    double x_max = *std::max_element(x_vals.begin(), x_vals.end());
    double y_min = *std::min_element(y_vals.begin(), y_vals.end());
    double y_max = *std::max_element(y_vals.begin(), y_vals.end());
    if (x_max == x_min) x_max = x_min + 1.0;
    if (y_max == y_min) y_max = y_min + 1.0;

    // Map to SVG coords: origin at top-left, y increases downward.
    auto sx = [&](double x) { return mx + (x - x_min) / (x_max - x_min) * pw; };
    auto sy = [&](double y) { return H - by - (y - y_min) / (y_max - y_min) * ph; };

    std::ostringstream s;
    s << "<svg width=\"" << W << "\" height=\"" << H
      << "\" style=\"font:12px sans-serif;border:1px solid #ddd;border-radius:4px\">\n";

    // Axes
    s << "<line x1=\"" << mx << "\" y1=\"" << my << "\" x2=\"" << mx
      << "\" y2=\"" << H - by << "\" stroke=\"#888\" stroke-width=\"1\"/>\n";
    s << "<line x1=\"" << mx << "\" y1=\"" << H - by << "\" x2=\"" << W - bx
      << "\" y2=\"" << H - by << "\" stroke=\"#888\" stroke-width=\"1\"/>\n";

    // Y-axis label (rotated)
    s << "<text transform=\"rotate(-90," << (mx - 45) << "," << (H / 2)
      << ")\" x=\"" << (mx - 45) << "\" y=\"" << (H / 2)
      << "\" text-anchor=\"middle\" fill=\"#444\" font-size=\"11\">"
      << y_label << "</text>\n";

    // X-axis label
    s << "<text x=\"" << (mx + pw / 2) << "\" y=\"" << (H - 3)
      << "\" text-anchor=\"middle\" fill=\"#444\" font-size=\"11\">"
      << x_label << "</text>\n";

    // Y tick (max)
    s << "<text x=\"" << (mx - 4) << "\" y=\"" << (my + 4)
      << "\" text-anchor=\"end\" fill=\"#666\" font-size=\"10\">";
    s << std::fixed << std::setprecision(2) << y_max << "</text>\n";

    // Y tick (min)
    s << "<text x=\"" << (mx - 4) << "\" y=\"" << (H - by + 4)
      << "\" text-anchor=\"end\" fill=\"#666\" font-size=\"10\">";
    s << std::fixed << std::setprecision(2) << y_min << "</text>\n";

    // X ticks (first and last)
    s << "<text x=\"" << mx << "\" y=\"" << (H - by + 14)
      << "\" text-anchor=\"middle\" fill=\"#666\" font-size=\"10\">";
    s << std::fixed << std::setprecision(0) << x_min << "</text>\n";
    s << "<text x=\"" << (mx + pw) << "\" y=\"" << (H - by + 14)
      << "\" text-anchor=\"middle\" fill=\"#666\" font-size=\"10\">";
    s << std::fixed << std::setprecision(0) << x_max << "</text>\n";

    // Polyline
    s << "<polyline fill=\"none\" stroke=\"" << color
      << "\" stroke-width=\"2\" points=\"";
    for (std::size_t i = 0; i < x_vals.size(); ++i) {
        s << std::fixed << std::setprecision(1)
          << sx(x_vals[i]) << "," << sy(y_vals[i]) << " ";
    }
    s << "\"/>\n</svg>";
    return s.str();
}

} // anonymous namespace

// ============================================================================
// write_fire_report_html
//
// Generates a self-contained HTML fire report.
//
// Parameters:
//   inputs      – InputParameters (contains fire_report_file, fire_stats_file,
//                 rothermel, domain geometry, etc.)
//   final_time_s – total simulated time [s]
//   final_step   – number of timesteps completed
// ============================================================================
inline void write_fire_report_html(const InputParameters& inputs,
                                   double final_time_s,
                                   int    final_step)
{
    if (inputs.fire_report_file.empty()) return;
    if (!ParallelDescriptor::IOProcessor()) return;

    // Read fire stats CSV
    const std::vector<FireStatsRow> rows =
        read_fire_stats_csv(inputs.fire_stats_file);

    // Summary values from last row
    double burned_ha_final   = 0.0, perim_km_final = 0.0;
    double max_ros_ms        = 0.0, cum_co2        = 0.0;
    if (!rows.empty()) {
        burned_ha_final = rows.back().burned_ha;
        perim_km_final  = rows.back().perimeter_km;
        for (const auto& r : rows) {
            if (r.head_ros_ms > max_ros_ms) max_ros_ms = r.head_ros_ms;
            cum_co2 += r.co2_kg_m2;
        }
    }

    // Build time-series vectors for SVG charts
    std::vector<double> t_hrs, area_ha, perim_km_v, ros_v, co2_v;
    for (const auto& r : rows) {
        double t_h = r.time_s / 3600.0;
        t_hrs.push_back(t_h);
        area_ha.push_back(r.burned_ha);
        perim_km_v.push_back(r.perimeter_km);
        ros_v.push_back(r.head_ros_ms * 60.0); // m/s → m/min
        co2_v.push_back(r.co2_kg_m2);
    }

    // Generate SVG charts
    const std::string chart_area  = make_svg_chart(t_hrs, area_ha,    "Time (h)", "Area (ha)",     "#e63946");
    const std::string chart_perim = make_svg_chart(t_hrs, perim_km_v, "Time (h)", "Perimeter (km)","#457b9d");
    const std::string chart_ros   = make_svg_chart(t_hrs, ros_v,      "Time (h)", "Head ROS (m/min)","#f4a261");
    const std::string chart_co2   = make_svg_chart(t_hrs, co2_v,      "Time (h)", "CO₂ (kg/m²)",   "#2a9d8f");

    // Domain dimensions
    const double dom_x_km = (inputs.phi_x - inputs.plo_x) / 1000.0;
    const double dom_y_km = (inputs.phi_y - inputs.plo_y) / 1000.0;
    const double sim_hrs  = final_time_s / 3600.0;

    std::ofstream html(inputs.fire_report_file);
    if (!html.is_open()) {
        amrex::Print() << "WARNING: Cannot write fire report HTML: "
                       << inputs.fire_report_file << "\n";
        return;
    }

    // ---- HTML document ----
    html << R"(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Wildfire Simulation Report</title>
<style>
  body{font-family:'Segoe UI',Arial,sans-serif;margin:0;background:#f7f7f7;color:#222;}
  header{background:#b71c1c;color:#fff;padding:18px 32px;}
  header h1{margin:0;font-size:1.6em;letter-spacing:.5px;}
  header p{margin:4px 0 0;font-size:.95em;opacity:.85;}
  main{max-width:1000px;margin:24px auto;padding:0 16px;}
  section{background:#fff;border-radius:8px;box-shadow:0 1px 4px rgba(0,0,0,.1);
          margin-bottom:24px;padding:20px 24px;}
  h2{margin:0 0 14px;font-size:1.15em;color:#b71c1c;border-bottom:2px solid #b71c1c;padding-bottom:6px;}
  .kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;}
  .kpi{background:#fff5f5;border:1px solid #f5c6c6;border-radius:6px;padding:12px 16px;text-align:center;}
  .kpi .val{font-size:1.55em;font-weight:700;color:#c0392b;}
  .kpi .lbl{font-size:.82em;color:#666;margin-top:4px;}
  table{border-collapse:collapse;width:100%;font-size:.88em;}
  th{background:#b71c1c;color:#fff;padding:6px 10px;text-align:right;}
  th:first-child{text-align:left;}
  td{padding:5px 10px;border-bottom:1px solid #eee;text-align:right;}
  td:first-child{text-align:left;}
  tr:nth-child(even){background:#fafafa;}
  .param-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:8px;}
  .param-row{display:flex;justify-content:space-between;border-bottom:1px solid #f0f0f0;padding:4px 0;}
  .param-key{color:#555;font-size:.88em;}
  .param-val{font-weight:600;font-size:.88em;}
  .chart-grid{display:grid;grid-template-columns:1fr 1fr;gap:16px;}
  @media(max-width:700px){.chart-grid{grid-template-columns:1fr;}}
  .chart-box{background:#fafafa;border-radius:6px;padding:12px;text-align:center;}
  .chart-box h3{margin:0 0 8px;font-size:.95em;color:#555;}
  footer{text-align:center;font-size:.8em;color:#aaa;margin:32px 0 16px;}
</style>
</head>
<body>
<header>
  <h1>🔥 Wildfire Simulation Report</h1>
  <p>Generated by wildfire_levelset | Propagation: )";
    html << inputs.propagation_method << " | Model: " << inputs.fire_spread_model;
    html << R"(</p>
</header>
<main>

<!-- KPI cards -->
<section>
  <h2>Final Fire Statistics</h2>
  <div class="kpi-grid">
    <div class="kpi"><div class="val">)";
    html << std::fixed << std::setprecision(1) << burned_ha_final;
    html << R"(</div><div class="lbl">Burned Area (ha)</div></div>
    <div class="kpi"><div class="val">)";
    html << std::fixed << std::setprecision(2) << perim_km_final;
    html << R"(</div><div class="lbl">Perimeter (km)</div></div>
    <div class="kpi"><div class="val">)";
    html << std::fixed << std::setprecision(2) << sim_hrs;
    html << R"(</div><div class="lbl">Simulated Time (h)</div></div>
    <div class="kpi"><div class="val">)";
    html << final_step;
    html << R"(</div><div class="lbl">Time Steps</div></div>
    <div class="kpi"><div class="val">)";
    html << std::fixed << std::setprecision(3) << max_ros_ms;
    html << R"(</div><div class="lbl">Max Head ROS (m/s)</div></div>
    <div class="kpi"><div class="val">)";
    html << std::fixed << std::setprecision(4) << cum_co2;
    html << R"(</div><div class="lbl">Cumul. CO₂ (kg/m²)</div></div>
  </div>
</section>

<!-- Time-series charts -->
<section>
  <h2>Time Series</h2>
  <div class="chart-grid">
    <div class="chart-box"><h3>Burned Area</h3>)" << chart_area  << R"(</div>
    <div class="chart-box"><h3>Fire Perimeter</h3>)" << chart_perim << R"(</div>
    <div class="chart-box"><h3>Head ROS</h3>)" << chart_ros   << R"(</div>
    <div class="chart-box"><h3>CO₂ Emissions</h3>)" << chart_co2   << R"(</div>
  </div>
</section>

<!-- Simulation parameters -->
<section>
  <h2>Simulation Parameters</h2>
  <div class="param-grid">
    <div>
)";
    // Domain
    auto pr = [&](const std::string& k, const std::string& v) {
        html << "      <div class=\"param-row\"><span class=\"param-key\">" << k
             << "</span><span class=\"param-val\">" << v << "</span></div>\n";
    };
    auto prd = [&](const std::string& k, double v, int prec = 3) {
        std::ostringstream ss; ss << std::fixed << std::setprecision(prec) << v;
        pr(k, ss.str());
    };
    auto pri = [&](const std::string& k, int v) {
        pr(k, std::to_string(v));
    };

    pr("Propagation method", inputs.propagation_method);
    pr("Fire spread model",  inputs.fire_spread_model);
    prd("Domain X [km]",    dom_x_km, 2);
    prd("Domain Y [km]",    dom_y_km, 2);
    pri("Grid nx",          inputs.n_cell_x);
    pri("Grid ny",          inputs.n_cell_y);
    pri("Total steps",      final_step);
    prd("Final time [h]",   sim_hrs, 2);
    prd("Wind ux [m/s]",    static_cast<double>(inputs.ux), 2);
    prd("Wind uy [m/s]",    static_cast<double>(inputs.uy), 2);

    html << "    </div>\n    <div>\n";
    prd("w₀ [lb/ft²]",     static_cast<double>(inputs.rothermel.w0), 4);
    prd("sigma [ft⁻¹]",    static_cast<double>(inputs.rothermel.sigma), 1);
    prd("delta [ft]",       static_cast<double>(inputs.rothermel.delta), 2);
    prd("M_f (1-hr)",       static_cast<double>(inputs.rothermel.M_f), 3);
    prd("M_x (extinctn)",   static_cast<double>(inputs.rothermel.M_x), 3);
    prd("M_lh (live herb)", static_cast<double>(inputs.rothermel.M_lh), 3);
    prd("M_lw (live wood)", static_cast<double>(inputs.rothermel.M_lw), 3);
    prd("slope_x",          static_cast<double>(inputs.rothermel.slope_x), 4);
    prd("slope_y",          static_cast<double>(inputs.rothermel.slope_y), 4);
    pr("Slope×wind cross",  inputs.rothermel.use_slope_wind_cross ? "enabled" : "disabled");

    html << "    </div>\n    <div>\n";
    pr("Burn period",       inputs.burn_period.enable ? "enabled" : "disabled");
    if (inputs.burn_period.enable) {
        prd("  start hour",  static_cast<double>(inputs.burn_period.start_hour), 1);
        prd("  end hour",    static_cast<double>(inputs.burn_period.end_hour), 1);
    }
    pr("Spotting (basic)",  inputs.spotting.enable      ? "enabled" : "disabled");
    if (inputs.spotting.enable) {
        prd("  P_base",     static_cast<double>(inputs.spotting.P_base), 3);
        prd("  P_catch",    static_cast<double>(inputs.spotting.P_catch), 3);
    }
    pr("Spotting (Albini)", inputs.albini_spotting.enable ? "enabled" : "disabled");
    if (inputs.albini_spotting.enable) {
        prd("  P_base",     static_cast<double>(inputs.albini_spotting.P_base), 3);
        prd("  P_catch",    static_cast<double>(inputs.albini_spotting.P_catch), 3);
    }
    pr("Crown fire",        inputs.crown.enable ? "enabled" : "disabled");
    pr("Adj file",         !inputs.fuel_adj_file.empty() ? inputs.fuel_adj_file : "(none)");
    pr("Landscape file",   !inputs.rothermel.landscape_file.empty()
                             ? inputs.rothermel.landscape_file : "(none)");
    pr("Stats CSV",        !inputs.fire_stats_file.empty() ? inputs.fire_stats_file : "(none)");

    html << R"(    </div>
  </div>
</section>

<!-- Per-step table -->
<section>
  <h2>Per-Step Fire Statistics</h2>
  <div style="overflow-x:auto">
  <table>
    <thead><tr>
      <th>Step</th><th>Time (h)</th><th>Burned (ha)</th><th>Perim (km)</th>
      <th>Active cells</th><th>Head ROS (m/min)</th>
      <th>CO₂ (kg/m²)</th><th>CO (kg/m²)</th><th>PM₂.₅ (kg/m²)</th>
    </tr></thead>
    <tbody>
)";
    for (const auto& r : rows) {
        html << "    <tr>"
             << "<td>" << r.step << "</td>"
             << "<td>" << std::fixed << std::setprecision(2) << r.time_s / 3600.0 << "</td>"
             << "<td>" << std::fixed << std::setprecision(2) << r.burned_ha << "</td>"
             << "<td>" << std::fixed << std::setprecision(3) << r.perimeter_km << "</td>"
             << "<td>" << r.active_front << "</td>"
             << "<td>" << std::fixed << std::setprecision(4) << r.head_ros_ms * 60.0 << "</td>"
             << "<td>" << std::scientific << std::setprecision(3) << r.co2_kg_m2 << "</td>"
             << "<td>" << std::scientific << std::setprecision(3) << r.co_kg_m2 << "</td>"
             << "<td>" << std::scientific << std::setprecision(3) << r.pm25_kg_m2 << "</td>"
             << "</tr>\n";
    }
    html << R"(    </tbody>
  </table>
  </div>
</section>

</main>
<footer>wildfire_levelset | Report generated at end of simulation run</footer>
</body>
</html>
)";

    html.close();
    amrex::Print() << "Wrote HTML fire report to: " << inputs.fire_report_file << "\n";
}

#endif // FIRE_REPORT_H
