// velocity_field.H - Velocity field initialization functions
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <AMReX_Gpu.H>

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
    
    // Copy data to device-accessible memory
    Gpu::DeviceVector<Real> d_x(x_data.size());
    Gpu::DeviceVector<Real> d_y(y_data.size());
    Gpu::DeviceVector<Real> d_u(u_data.size());
    Gpu::DeviceVector<Real> d_v(v_data.size());
    
    Gpu::copy(Gpu::hostToDevice, x_data.begin(), x_data.end(), d_x.begin());
    Gpu::copy(Gpu::hostToDevice, y_data.begin(), y_data.end(), d_y.begin());
    Gpu::copy(Gpu::hostToDevice, u_data.begin(), u_data.end(), d_u.begin());
    Gpu::copy(Gpu::hostToDevice, v_data.begin(), v_data.end(), d_v.begin());
    
    Real const* x_ptr = d_x.data();
    Real const* y_ptr = d_y.data();
    Real const* u_ptr = d_u.data();
    Real const* v_ptr = d_v.data();
    int n_points = static_cast<int>(x_data.size());
    
    // Get domain bounds
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    // Interpolate onto MultiFab using inverse distance weighting
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute cell center coordinates
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            
            // Inverse distance weighting interpolation
            Real sum_weights = Real(0.0);
            Real sum_u = Real(0.0);
            Real sum_v = Real(0.0);
            Real epsilon = Real(1.0e-10);  // Small value to avoid division by zero
            bool exact_match = false;
            
            for (int n = 0; n < n_points; ++n) {
                Real dx_val = x - x_ptr[n];
                Real dy_val = y - y_ptr[n];
                Real dist_sq = dx_val * dx_val + dy_val * dy_val;
                
                // If very close to a data point, use its value directly
                if (dist_sq < epsilon * epsilon) {
                    v(i,j,k,0) = u_ptr[n];
                    v(i,j,k,1) = v_ptr[n];
                    v(i,j,k,2) = Real(0.0);
                    exact_match = true;
                    break;
                }
                
                // For power=2, weight = 1/dist^2 = 1/dist_sq (no sqrt needed)
                Real weight = Real(1.0) / dist_sq;
                sum_weights += weight;
                sum_u += weight * u_ptr[n];
                sum_v += weight * v_ptr[n];
            }
            
            // Apply weighted average if not an exact match
            if (!exact_match) {
                v(i,j,k,0) = sum_u / sum_weights;
                v(i,j,k,1) = sum_v / sum_weights;
                v(i,j,k,2) = Real(0.0);
            }
        });
    }
    
    amrex::Print() << "Velocity field interpolated from unstructured data\n";
}

// ======================= Load wind field data from file ==========================
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
                                        const std::vector<Real>& x_data2, const std::vector<Real>& /*y_data2*/,
                                        const std::vector<Real>& u_data2, const std::vector<Real>& v_data2,
                                        Real time_fraction)
{
    // Ensure both datasets have same size
    if (x_data1.size() != x_data2.size()) {
        amrex::Abort("Wind field files must have the same number of data points");
    }
    
    int n_points = static_cast<int>(x_data1.size());
    
    // Copy data to device-accessible memory
    Gpu::DeviceVector<Real> d_x(n_points);
    Gpu::DeviceVector<Real> d_y(n_points);
    Gpu::DeviceVector<Real> d_u1(n_points);
    Gpu::DeviceVector<Real> d_v1(n_points);
    Gpu::DeviceVector<Real> d_u2(n_points);
    Gpu::DeviceVector<Real> d_v2(n_points);
    
    Gpu::copy(Gpu::hostToDevice, x_data1.begin(), x_data1.end(), d_x.begin());
    Gpu::copy(Gpu::hostToDevice, y_data1.begin(), y_data1.end(), d_y.begin());
    Gpu::copy(Gpu::hostToDevice, u_data1.begin(), u_data1.end(), d_u1.begin());
    Gpu::copy(Gpu::hostToDevice, v_data1.begin(), v_data1.end(), d_v1.begin());
    Gpu::copy(Gpu::hostToDevice, u_data2.begin(), u_data2.end(), d_u2.begin());
    Gpu::copy(Gpu::hostToDevice, v_data2.begin(), v_data2.end(), d_v2.begin());
    
    Real const* x_ptr = d_x.data();
    Real const* y_ptr = d_y.data();
    Real const* u1_ptr = d_u1.data();
    Real const* v1_ptr = d_v1.data();
    Real const* u2_ptr = d_u2.data();
    Real const* v2_ptr = d_v2.data();
    
    // Get domain bounds
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    // Interpolate onto MultiFab using inverse distance weighting and temporal interpolation
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute cell center coordinates
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            
            // Inverse distance weighting interpolation for both time snapshots
            Real sum_weights = Real(0.0);
            Real sum_u1 = Real(0.0);
            Real sum_v1 = Real(0.0);
            Real sum_u2 = Real(0.0);
            Real sum_v2 = Real(0.0);
            Real epsilon = Real(1.0e-10);
            bool exact_match = false;
            Real u1_val = Real(0.0);
            Real v1_val = Real(0.0);
            Real u2_val = Real(0.0);
            Real v2_val = Real(0.0);
            
            for (int n = 0; n < n_points; ++n) {
                Real dx_val = x - x_ptr[n];
                Real dy_val = y - y_ptr[n];
                Real dist_sq = dx_val * dx_val + dy_val * dy_val;
                
                // If very close to a data point, use its value directly
                if (dist_sq < epsilon * epsilon) {
                    u1_val = u1_ptr[n];
                    v1_val = v1_ptr[n];
                    u2_val = u2_ptr[n];
                    v2_val = v2_ptr[n];
                    exact_match = true;
                    break;
                }
                
                // For power=2, weight = 1/dist^2
                Real weight = Real(1.0) / dist_sq;
                sum_weights += weight;
                sum_u1 += weight * u1_ptr[n];
                sum_v1 += weight * v1_ptr[n];
                sum_u2 += weight * u2_ptr[n];
                sum_v2 += weight * v2_ptr[n];
            }
            
            // Apply weighted average if not an exact match
            if (!exact_match) {
                u1_val = sum_u1 / sum_weights;
                v1_val = sum_v1 / sum_weights;
                u2_val = sum_u2 / sum_weights;
                v2_val = sum_v2 / sum_weights;
            }
            
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