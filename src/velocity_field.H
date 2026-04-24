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
#endif  // AMREX_SPACEDIM == 2