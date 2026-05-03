// velocity_field.H - Velocity field initialization functions
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <AMReX_Gpu.H>

#include "spatial_grid.H"

// ======================= Initialize velocity (constant) ==========================
static void init_velocity_constant(MultiFab& vel, const Geometry& geom, Real ux, Real uy, Real uz)
{
    (void)geom;  // Suppress unused parameter warning
#if (AMREX_SPACEDIM == 2)
    (void)uz;  // uz not used in 2D
#endif
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            v(i,j,k,0) = ux;
            v(i,j,k,1) = uy;
#if (AMREX_SPACEDIM == 3)
            v(i,j,k,2) = uz;
#else
            v(i,j,k,2) = Real(0.0);  // In 2D, set z-component to 0
#endif
        });
    }
}

// ======================= Initialize velocity from unstructured 2D data ===========
#if (AMREX_SPACEDIM == 2)
static void init_velocity_from_file(MultiFab& vel, const Geometry& geom, const std::string& filename)
{
    // Read unstructured data from file (X Y U V format)
    std::vector<Real> x_data, y_data, u_data, v_data;
    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open velocity file: " + filename);
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        Real x, y, u, v;
        if (iss >> x >> y >> u >> v) {
            x_data.push_back(x);
            y_data.push_back(y);
            u_data.push_back(u);
            v_data.push_back(v);
        }
    }
    infile.close();
    
    if (x_data.empty()) {
        amrex::Abort("No data read from velocity file: " + filename);
    }

    amrex::Print() << "Read " << x_data.size() << " velocity data points from " << filename << "\n";

    int n_points = static_cast<int>(x_data.size());

    // Characteristic data spacing
    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    // Build cell-list spatial index
    Gpu::DeviceVector<Real> d_sx, d_sy, d_su, d_sv;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> u_sorted = reorder_values(u_data, sidx);
    std::vector<Real> v_sorted = reorder_values(v_data, sidx);
    d_su.resize(n_points);
    d_sv.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, u_sorted.begin(), u_sorted.end(), d_su.begin());
    Gpu::copy(Gpu::hostToDevice, v_sorted.begin(), v_sorted.end(), d_sv.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* u_ptr = d_su.data();
    Real const* v_ptr = d_sv.data();

    // Get domain bounds
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // Interpolate onto MultiFab using Wendland C² RBF with cell-list index
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];

            v(i,j,k,0) = cell_list_interp(cl, u_ptr, x, y);
            v(i,j,k,1) = cell_list_interp(cl, v_ptr, x, y);
            v(i,j,k,2) = Real(0.0);
        });
    }

    amrex::Print() << "Velocity field interpolated from unstructured data\n";
}

// ======================= Update time-dependent velocity field ====================
static void load_wind_field_data(const std::string& filename,
                                  std::vector<Real>& x_data,
                                  std::vector<Real>& y_data,
                                  std::vector<Real>& u_data,
                                  std::vector<Real>& v_data)
{
    x_data.clear();
    y_data.clear();
    u_data.clear();
    v_data.clear();
    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open velocity file: " + filename);
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        Real x, y, u, v;
        if (iss >> x >> y >> u >> v) {
            x_data.push_back(x);
            y_data.push_back(y);
            u_data.push_back(u);
            v_data.push_back(v);
        }
    }
    infile.close();
    
    if (x_data.empty()) {
        amrex::Abort("No data read from velocity file: " + filename);
    }
    
    amrex::Print() << "Read " << x_data.size() << " velocity data points from " << filename << "\n";
}

// ======================= Interpolate velocity field using two time snapshots ====
static void interpolate_velocity_field(MultiFab& vel, const Geometry& geom,
                                        const std::vector<Real>& x_data1, const std::vector<Real>& y_data1,
                                        const std::vector<Real>& u_data1, const std::vector<Real>& v_data1,
                                        const std::vector<Real>& x_data2, const std::vector<Real>& y_data2,
                                        const std::vector<Real>& u_data2, const std::vector<Real>& v_data2,
                                        Real time_fraction)
{
    // Ensure both datasets have same size
    if (x_data1.size() != x_data2.size()) {
        amrex::Abort("Wind field files must have the same number of data points");
    }

    int n_points = static_cast<int>(x_data1.size());

    // Characteristic data spacing
    Real x_min_d = *std::min_element(x_data1.begin(), x_data1.end());
    Real x_max_d = *std::max_element(x_data1.begin(), x_data1.end());
    Real y_min_d = *std::min_element(y_data1.begin(), y_data1.end());
    Real y_max_d = *std::max_element(y_data1.begin(), y_data1.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    // Build cell-list spatial index (positions taken from snapshot 1)
    Gpu::DeviceVector<Real> d_sx, d_sy, d_su1, d_sv1, d_su2, d_sv2;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data1, y_data1, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> u1_sorted = reorder_values(u_data1, sidx);
    std::vector<Real> v1_sorted = reorder_values(v_data1, sidx);
    std::vector<Real> u2_sorted = reorder_values(u_data2, sidx);
    std::vector<Real> v2_sorted = reorder_values(v_data2, sidx);

    d_su1.resize(n_points); d_sv1.resize(n_points);
    d_su2.resize(n_points); d_sv2.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, u1_sorted.begin(), u1_sorted.end(), d_su1.begin());
    Gpu::copy(Gpu::hostToDevice, v1_sorted.begin(), v1_sorted.end(), d_sv1.begin());
    Gpu::copy(Gpu::hostToDevice, u2_sorted.begin(), u2_sorted.end(), d_su2.begin());
    Gpu::copy(Gpu::hostToDevice, v2_sorted.begin(), v2_sorted.end(), d_sv2.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* u1_ptr = d_su1.data();
    Real const* v1_ptr = d_sv1.data();
    Real const* u2_ptr = d_su2.data();
    Real const* v2_ptr = d_sv2.data();

    // Get domain bounds
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // Interpolate onto MultiFab using Wendland C² RBF with cell-list and temporal blend
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];

            Real u1_val = cell_list_interp(cl, u1_ptr, x, y);
            Real v1_val = cell_list_interp(cl, v1_ptr, x, y);
            Real u2_val = cell_list_interp(cl, u2_ptr, x, y);
            Real v2_val = cell_list_interp(cl, v2_ptr, x, y);

            // Linear interpolation in time
            v(i,j,k,0) = u1_val + time_fraction * (u2_val - u1_val);
            v(i,j,k,1) = v1_val + time_fraction * (v2_val - v1_val);
            v(i,j,k,2) = Real(0.0);
        });
    }
}

// ======================= Update time-dependent velocity field ====================
static void update_time_dependent_velocity(MultiFab& vel, const Geometry& geom,
                                            const std::string& base_filename,
                                            Real current_time,
                                            Real wind_time_spacing,
                                            std::vector<Real>& x_data1,
                                            std::vector<Real>& y_data1,
                                            std::vector<Real>& u_data1,
                                            std::vector<Real>& v_data1,
                                            std::vector<Real>& x_data2,
                                            std::vector<Real>& y_data2,
                                            std::vector<Real>& u_data2,
                                            std::vector<Real>& v_data2,
                                            int& current_field_index,
                                            int& next_field_index)
{
    // Determine which wind field indices to use based on current_time
    int new_current_index = static_cast<int>(current_time / wind_time_spacing);
    int new_next_index = new_current_index + 1;
    
    // Check if we need to load new wind field data
    bool need_reload = false;
    if (new_current_index != current_field_index || new_next_index != next_field_index) {
        need_reload = true;
        current_field_index = new_current_index;
        next_field_index = new_next_index;
    }
    
    if (need_reload) {
        // Construct filenames
        std::string file1, file2;
        
        // Extract base filename (remove .csv extension if present)
        std::string base = base_filename;
        size_t csv_pos = base.rfind(".csv");
        if (csv_pos != std::string::npos) {
            base = base.substr(0, csv_pos);
        }
        
        // First file: base filename if index is 0, otherwise base_N.csv
        if (current_field_index == 0) {
            file1 = base_filename;
        } else {
            file1 = base + "_" + std::to_string(current_field_index) + ".csv";
        }
        
        // Second file: base_N.csv
        file2 = base + "_" + std::to_string(next_field_index) + ".csv";
        
        amrex::Print() << "Loading wind fields: " << file1 << " and " << file2 
                       << " at time " << current_time << " s\n";
        
        // Load both wind field files
        load_wind_field_data(file1, x_data1, y_data1, u_data1, v_data1);
        load_wind_field_data(file2, x_data2, y_data2, u_data2, v_data2);
    }
    
    // Calculate time fraction within the current interval
    Real time_in_interval = current_time - (current_field_index * wind_time_spacing);
    Real time_fraction = time_in_interval / wind_time_spacing;
    
    // Clamp time_fraction to [0, 1]
    time_fraction = std::max(Real(0.0), std::min(Real(1.0), time_fraction));
    
    // Interpolate velocity field
    interpolate_velocity_field(vel, geom, 
                                x_data1, y_data1, u_data1, v_data1,
                                x_data2, y_data2, u_data2, v_data2,
                                time_fraction);
}
#endif  // AMREX_SPACEDIM == 2

// ======================= PLT Wind: Scan directory for file timestamps =============
// Probes for plt00000, plt00001, ... in plt_wind_dir and reads the simulation
// time stored in each file's Header.  Stops at the first missing index.
#include <AMReX_PlotFileUtil.H>

static void scan_plt_wind_files(const std::string& plt_wind_dir,
                                 std::vector<Real>& timestamps)
{
    timestamps.clear();
    for (int i = 0; ; ++i) {
        char hdr_path[512];
        std::snprintf(hdr_path, sizeof(hdr_path), "%s/plt%05d/Header",
                      plt_wind_dir.c_str(), i);
        std::ifstream hdr_test(hdr_path);
        if (!hdr_test.is_open()) break;
        hdr_test.close();

        char plt_path[512];
        std::snprintf(plt_path, sizeof(plt_path), "%s/plt%05d",
                      plt_wind_dir.c_str(), i);
        std::string plt_path_str(plt_path);
        PlotFileData pf(plt_path_str);
        timestamps.push_back(static_cast<Real>(pf.time()));
    }

    if (timestamps.empty()) {
        amrex::Abort("No plt wind files found in directory: " + plt_wind_dir);
    }

    amrex::Print() << "PLT wind: found " << timestamps.size()
                   << " files in " << plt_wind_dir << "\n";
    for (int i = 0; i < static_cast<int>(timestamps.size()); ++i) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "plt%05d", i);
        amrex::Print() << "  " << buf << "  t = " << timestamps[i] << " s\n";
    }
}

// ======================= PLT Wind: Extract 2D wind slice from a plt file ==========
// For AMREX_SPACEDIM==3 builds the data are sliced at the k-index that
// corresponds to fire_height in the plt file's coordinate system.
// For AMREX_SPACEDIM==2 builds the plt files must already be 2-D and
// fire_height is ignored.
static void extract_wind_from_plt(const std::string& plt_wind_dir, int file_idx,
                                    const std::string& u_var, const std::string& v_var,
                                    Real fire_height,
                                    std::vector<Real>& x_data,
                                    std::vector<Real>& y_data,
                                    std::vector<Real>& u_data,
                                    std::vector<Real>& v_data)
{
    char plt_path[512];
    std::snprintf(plt_path, sizeof(plt_path), "%s/plt%05d",
                  plt_wind_dir.c_str(), file_idx);
    std::string plt_path_str(plt_path);

    PlotFileData pf(plt_path_str);

    MultiFab mf_u = pf.get(0, u_var);
    MultiFab mf_v = pf.get(0, v_var);

    const auto problo = pf.probLo();
    const auto dx     = pf.cellSize(0);

    x_data.clear(); y_data.clear(); u_data.clear(); v_data.clear();

#if (AMREX_SPACEDIM == 3)
    // Find the k-index in the plt file that sits closest to fire_height
    Box domain0 = pf.probDomain(0);
    int target_k = static_cast<int>((fire_height - problo[2]) / dx[2]);
    int k_lo = domain0.smallEnd(2);
    int k_hi = domain0.bigEnd(2);
    target_k = std::max(k_lo, std::min(k_hi, target_k));

    amrex::Print() << "PLT wind: slicing " << plt_path_str
                   << " at k=" << target_k
                   << " (fire_height=" << fire_height << " m)\n";

    for (MFIter mfi(mf_u); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (target_k < bx.smallEnd(2) || target_k > bx.bigEnd(2)) continue;

        Array4<Real const> u_arr = mf_u.const_array(mfi);
        Array4<Real const> v_arr = mf_v.const_array(mfi);

        const int k = target_k;
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                Real x = problo[0] + (i + Real(0.5)) * dx[0];
                Real y = problo[1] + (j + Real(0.5)) * dx[1];
                x_data.push_back(x);
                y_data.push_back(y);
                u_data.push_back(u_arr(i, j, k, 0));
                v_data.push_back(v_arr(i, j, k, 0));
            }
        }
    }
#else
    // 2-D plt files: use all cells directly; fire_height is not used for slicing
    (void)fire_height;
    for (MFIter mfi(mf_u); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> u_arr = mf_u.const_array(mfi);
        Array4<Real const> v_arr = mf_v.const_array(mfi);

        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                Real x = problo[0] + (i + Real(0.5)) * dx[0];
                Real y = problo[1] + (j + Real(0.5)) * dx[1];
                x_data.push_back(x);
                y_data.push_back(y);
                u_data.push_back(u_arr(i, j, 0, 0));
                v_data.push_back(v_arr(i, j, 0, 0));
            }
        }
    }
#endif

    if (x_data.empty()) {
        amrex::Abort("No wind data extracted from plt file: " + plt_path_str);
    }

    amrex::Print() << "PLT wind: extracted " << x_data.size()
                   << " points from " << plt_path_str << "\n";
}

// ======================= PLT Wind: Interpolate scattered (x,y,u,v) to MultiFab ===
// Applies Wendland C² RBF spatial interpolation with temporal blending between
// two time snapshots.  Works for both 2-D and 3-D simulation domains: for 3-D
// builds the same (u,v) derived from the fire_height slice are applied at
// every k-level; the w-component is always set to zero.
static void interpolate_plt_wind_to_vel(MultiFab& vel, const Geometry& geom,
                                         const std::vector<Real>& x_data1,
                                         const std::vector<Real>& y_data1,
                                         const std::vector<Real>& u_data1,
                                         const std::vector<Real>& v_data1,
                                         const std::vector<Real>& x_data2,
                                         const std::vector<Real>& y_data2,
                                         const std::vector<Real>& u_data2,
                                         const std::vector<Real>& v_data2,
                                         Real time_fraction)
{
    int n_points = static_cast<int>(x_data1.size());

    Real x_min_d = *std::min_element(x_data1.begin(), x_data1.end());
    Real x_max_d = *std::max_element(x_data1.begin(), x_data1.end());
    Real y_min_d = *std::min_element(y_data1.begin(), y_data1.end());
    Real y_max_d = *std::max_element(y_data1.begin(), y_data1.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    Gpu::DeviceVector<Real> d_sx, d_sy, d_su1, d_sv1, d_su2, d_sv2;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data1, y_data1, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> u1_sorted = reorder_values(u_data1, sidx);
    std::vector<Real> v1_sorted = reorder_values(v_data1, sidx);
    std::vector<Real> u2_sorted = reorder_values(u_data2, sidx);
    std::vector<Real> v2_sorted = reorder_values(v_data2, sidx);

    d_su1.resize(n_points); d_sv1.resize(n_points);
    d_su2.resize(n_points); d_sv2.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, u1_sorted.begin(), u1_sorted.end(), d_su1.begin());
    Gpu::copy(Gpu::hostToDevice, v1_sorted.begin(), v1_sorted.end(), d_sv1.begin());
    Gpu::copy(Gpu::hostToDevice, u2_sorted.begin(), u2_sorted.end(), d_su2.begin());
    Gpu::copy(Gpu::hostToDevice, v2_sorted.begin(), v2_sorted.end(), d_sv2.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* u1_ptr = d_su1.data();
    Real const* v1_ptr = d_sv1.data();
    Real const* u2_ptr = d_su2.data();
    Real const* v2_ptr = d_sv2.data();

    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];

            Real u1_val = cell_list_interp(cl, u1_ptr, x, y);
            Real v1_val = cell_list_interp(cl, v1_ptr, x, y);
            Real u2_val = cell_list_interp(cl, u2_ptr, x, y);
            Real v2_val = cell_list_interp(cl, v2_ptr, x, y);

            v(i,j,k,0) = u1_val + time_fraction * (u2_val - u1_val);
            v(i,j,k,1) = v1_val + time_fraction * (v2_val - v1_val);
            v(i,j,k,2) = Real(0.0);
        });
    }
}

// ======================= PLT Wind: Update velocity field each time step ===========
// Manages the two-snapshot cache and drives temporal interpolation.  On the
// first call timestamps is populated by scanning plt_wind_dir.
static void update_plt_wind_velocity(MultiFab& vel, const Geometry& geom,
                                      const std::string& plt_wind_dir,
                                      const std::string& u_var,
                                      const std::string& v_var,
                                      Real fire_height,
                                      Real current_time,
                                      // state cached between calls:
                                      std::vector<Real>& timestamps,
                                      std::vector<Real>& x_data1,
                                      std::vector<Real>& y_data1,
                                      std::vector<Real>& u_data1,
                                      std::vector<Real>& v_data1,
                                      std::vector<Real>& x_data2,
                                      std::vector<Real>& y_data2,
                                      std::vector<Real>& u_data2,
                                      std::vector<Real>& v_data2,
                                      int& loaded_idx1,
                                      int& loaded_idx2)
{
    // Populate timestamps on first call
    if (timestamps.empty()) {
        scan_plt_wind_files(plt_wind_dir, timestamps);
    }

    const int n_files = static_cast<int>(timestamps.size());

    // Find the pair of file indices that bracket current_time
    int idx1 = 0, idx2 = 0;
    if (n_files == 1 || current_time <= timestamps[0]) {
        idx1 = idx2 = 0;
    } else if (current_time >= timestamps[n_files - 1]) {
        idx1 = idx2 = n_files - 1;
    } else {
        for (int i = 0; i < n_files - 1; ++i) {
            if (current_time >= timestamps[i] && current_time < timestamps[i + 1]) {
                idx1 = i;
                idx2 = i + 1;
                break;
            }
        }
    }

    // Reload wind data only when the bracketing pair changes
    if (idx1 != loaded_idx1 || idx2 != loaded_idx2) {
        char b1[32], b2[32];
        std::snprintf(b1, sizeof(b1), "plt%05d", idx1);
        std::snprintf(b2, sizeof(b2), "plt%05d", idx2);
        amrex::Print() << "PLT wind: loading " << b1 << " and " << b2
                       << " at t=" << current_time << " s\n";

        extract_wind_from_plt(plt_wind_dir, idx1, u_var, v_var, fire_height,
                               x_data1, y_data1, u_data1, v_data1);
        if (idx1 != idx2) {
            extract_wind_from_plt(plt_wind_dir, idx2, u_var, v_var, fire_height,
                                   x_data2, y_data2, u_data2, v_data2);
        } else {
            x_data2 = x_data1;
            y_data2 = y_data1;
            u_data2 = u_data1;
            v_data2 = v_data1;
        }

        loaded_idx1 = idx1;
        loaded_idx2 = idx2;
    }

    // Time fraction within the current interval
    Real time_fraction = Real(0.0);
    if (idx1 != idx2 && timestamps[idx2] > timestamps[idx1]) {
        time_fraction = (current_time - timestamps[idx1])
                      / (timestamps[idx2] - timestamps[idx1]);
        time_fraction = std::max(Real(0.0), std::min(Real(1.0), time_fraction));
    }

    interpolate_plt_wind_to_vel(vel, geom,
                                 x_data1, y_data1, u_data1, v_data1,
                                 x_data2, y_data2, u_data2, v_data2,
                                 time_fraction);
}