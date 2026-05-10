// ============================================================================
// plt_wind_reader.H
// Read 3-D wind data from an AMReX plotfile produced by massconsistent_amr
// (https://github.com/hgopalan/massconsistent_amr) and expose the result as
// one-dimensional GPU arrays of (x, y, z, u, v, w).
//
// The plotfile is a standard AMReX single-level plotfile written by
//   WriteSingleLevelPlotfile(plot_file, output, {"u","v","w",...}, geom, ...)
//
// This reader does NOT use AMReX MultiFab / Box / Geometry types so that it
// works inside the 2-D wildfire_levelset build (AMREX_SPACEDIM == 2) even
// though the source wind data live on a 3-D grid.
//
// For Albini spotting (which is a 2-D horizontal trajectory) the 3-D wind is
// projected to horizontal components by averaging u and v over all z levels
// at each (x, y) column.
// ============================================================================

#ifndef PLT_WIND_READER_H
#define PLT_WIND_READER_H

#include <AMReX_Gpu.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

// ============================================================================
// PltWindData
//
// Host-side flat 1-D arrays (one entry per 3-D grid point) and their device
// copies.  Populated by read_plt_wind_file().
// ============================================================================
struct PltWindData {
    // Coordinates and velocity components at every 3-D grid point
    std::vector<amrex::Real> x, y, z;
    std::vector<amrex::Real> u, v, w;

    // Device copies (usable in GPU kernels via raw pointers)
    amrex::Gpu::DeviceVector<amrex::Real> d_x, d_y, d_z;
    amrex::Gpu::DeviceVector<amrex::Real> d_u, d_v, d_w;

    int n_points = 0;   // total number of points ( = nx*ny*nz )

    // Column-averaged horizontal wind on the 2-D (x,y) grid
    // Size: nx_2d * ny_2d
    std::vector<amrex::Real> u2d, v2d;
    amrex::Gpu::DeviceVector<amrex::Real> d_u2d, d_v2d;
    int nx_2d = 0, ny_2d = 0;
    amrex::Real lo2d[2] = {0.0, 0.0};
    amrex::Real hi2d[2] = {1.0, 1.0};

    bool valid = false;
};

// ============================================================================
// Internal helper: trim trailing whitespace / carriage-returns from a string
// ============================================================================
static inline void trim_string(std::string& s)
{
    while (!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\n'))
        s.pop_back();
}

// ============================================================================
// Internal helper: parse the plt Header file.
//
// Returns false if parsing fails.
// ============================================================================
static bool parse_plt_header_3d(
    const std::string& plt_dir,
    std::vector<std::string>& var_names,
    double prob_lo[3], double prob_hi[3],
    int& nx, int& ny, int& nz,
    std::string& level0_cell_path)
{
    const std::string header_path = plt_dir + "/Header";
    std::ifstream hf(header_path);
    if (!hf.is_open()) {
        amrex::Print() << "[plt_wind_reader] Cannot open Header: "
                       << header_path << "\n";
        return false;
    }

    // Line 1: version string (e.g. "HyperCLaw-V1.1" or "AMReX Version...")
    std::string version;
    std::getline(hf, version);

    // Line 2: number of plotfile variables
    int nvar = 0;
    hf >> nvar;
    hf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Lines 3..2+nvar: variable names
    var_names.clear();
    var_names.reserve(nvar);
    for (int n = 0; n < nvar; ++n) {
        std::string vname;
        std::getline(hf, vname);
        trim_string(vname);
        var_names.push_back(vname);
    }

    // Space dimension
    int sdim = 3;
    hf >> sdim;

    // Time
    double ptime = 0.0;
    hf >> ptime;

    // Finest level
    int finest_level = 0;
    hf >> finest_level;

    // prob_lo (sdim values)
    for (int d = 0; d < 3; ++d) prob_lo[d] = 0.0;
    for (int d = 0; d < sdim; ++d) hf >> prob_lo[d];

    // prob_hi (sdim values)
    for (int d = 0; d < 3; ++d) prob_hi[d] = 1.0;
    for (int d = 0; d < sdim; ++d) hf >> prob_hi[d];

    // Refinement ratios (finest_level values; may be absent if finest_level==0)
    for (int lev = 0; lev < finest_level; ++lev) {
        int rr; hf >> rr;
    }

    // Bounding boxes per level (one line each)
    hf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    nx = 1; ny = 1; nz = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        std::string bline;
        std::getline(hf, bline);
        if (lev == 0) {
            // Format example: "((0,0,0) (63,63,31) (0,0,0))"
            int lo0=0,lo1=0,lo2=0, hi0=0,hi1=0,hi2=0;
            if (std::sscanf(bline.c_str(), "((%d,%d,%d) (%d,%d,%d)",
                            &lo0, &lo1, &lo2, &hi0, &hi1, &hi2) == 6) {
                nx = hi0 - lo0 + 1;
                ny = hi1 - lo1 + 1;
                nz = hi2 - lo2 + 1;
            } else {
                // 2-D fallback (should not happen for massconsistent_amr data)
                int lo0_2=0,lo1_2=0, hi0_2=0,hi1_2=0;
                if (std::sscanf(bline.c_str(), "((%d,%d) (%d,%d)",
                                &lo0_2, &lo1_2, &hi0_2, &hi1_2) == 4) {
                    nx = hi0_2 - lo0_2 + 1;
                    ny = hi1_2 - lo1_2 + 1;
                    nz = 1;
                }
            }
        }
    }

    // Steps per level
    for (int lev = 0; lev <= finest_level; ++lev) {
        int st; hf >> st;
    }

    // Cell sizes per level per direction
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int d = 0; d < sdim; ++d) {
            double cs; hf >> cs;
        }
    }

    // Coordinate system, then 0 (reserved)
    int coord_sys = 0, zero_val = 0;
    hf >> coord_sys >> zero_val;

    // Per-level data directory paths
    hf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    level0_cell_path = "";
    for (int lev = 0; lev <= finest_level; ++lev) {
        std::string lpath;
        std::getline(hf, lpath);
        trim_string(lpath);
        if (lev == 0) {
            level0_cell_path = plt_dir + "/" + lpath;
        }
    }

    return (nx > 0 && ny > 0 && nz > 0 && !level0_cell_path.empty());
}

// ============================================================================
// Internal helper: parse a single-line AMReX FAB ASCII header
//
// Example line produced by WriteSingleLevelPlotfile:
//   FAB ((8, (64 11 52 0 1 12 0 1023)),(8,(64 11 52 0 1 12 0 1023))) ((0,0,0) (63,63,31) (0,0,0))
//
// We need the box extents and whether data is single or double precision.
// Returns false if parsing fails.
// ============================================================================
static bool parse_fab_header(const std::string& line,
                             int lo[3], int hi[3],
                             int& bytes_per_val)
{
    // Extract the box substring between the second "(" and the last ")"
    // We look for the box descriptor: "((lo0,lo1,lo2) (hi0,hi1,hi2) (0,0,0))"
    // which appears after the format descriptor.
    std::size_t pos = line.rfind("((");
    if (pos == std::string::npos) return false;
    const char* s = line.c_str() + pos;
    int lo0=0,lo1=0,lo2=0, hi0=0,hi1=0,hi2=0;
    if (std::sscanf(s, "((%d,%d,%d) (%d,%d,%d)", &lo0,&lo1,&lo2, &hi0,&hi1,&hi2) == 6) {
        lo[0]=lo0; lo[1]=lo1; lo[2]=lo2;
        hi[0]=hi0; hi[1]=hi1; hi[2]=hi2;
    } else if (std::sscanf(s, "((%d,%d) (%d,%d)", &lo0,&lo1, &hi0,&hi1) == 4) {
        lo[0]=lo0; lo[1]=lo1; lo[2]=0;
        hi[0]=hi0; hi[1]=hi1; hi[2]=0;
    } else {
        return false;
    }

    // Extract bytes per value (8 = double, 4 = float) from the format descriptor
    // e.g. "((8, ..." → 8 bytes
    bytes_per_val = 8;  // default to double
    std::size_t paren = line.find('(');
    if (paren != std::string::npos) {
        int bpv = 8;
        if (std::sscanf(line.c_str() + paren + 1, "(%d,", &bpv) == 1) {
            bytes_per_val = bpv;
        }
    }
    return true;
}

// ============================================================================
// read_plt_wind_file
//
// Opens the AMReX plt directory at `plt_dir`, looks for variable components
// named "u", "v", and "w", reads all 3-D grid data, and fills `pwd`.
//
// On success, pwd.valid == true and:
//   pwd.x[i], pwd.y[i], pwd.z[i] – physical coordinates of point i
//   pwd.u[i], pwd.v[i], pwd.w[i] – velocity at point i
//   pwd.n_points                  – total number of points
//   pwd.u2d, pwd.v2d             – height-averaged horizontal wind on 2-D grid
//   pwd.nx_2d, pwd.ny_2d         – 2-D grid dimensions
//
// Device (GPU) copies are populated for all arrays.
// ============================================================================
inline bool read_plt_wind_file(const std::string& plt_dir, PltWindData& pwd)
{
    pwd.valid = false;

    // ---- 1. Parse the plt Header ----
    std::vector<std::string> var_names;
    double prob_lo[3] = {0.0, 0.0, 0.0};
    double prob_hi[3] = {1.0, 1.0, 1.0};
    int nx = 0, ny = 0, nz = 0;
    std::string level0_cell_path;

    if (!parse_plt_header_3d(plt_dir, var_names, prob_lo, prob_hi,
                              nx, ny, nz, level0_cell_path)) {
        amrex::Print() << "[plt_wind_reader] Failed to parse plt Header in "
                       << plt_dir << "\n";
        return false;
    }

    amrex::Print() << "[plt_wind_reader] Domain: " << nx << "×" << ny << "×" << nz
                   << "  prob_lo=(" << prob_lo[0] << "," << prob_lo[1] << "," << prob_lo[2]
                   << ")  prob_hi=(" << prob_hi[0] << "," << prob_hi[1] << "," << prob_hi[2]
                   << ")\n";
    amrex::Print() << "[plt_wind_reader] Variables:";
    for (const auto& v : var_names) amrex::Print() << " " << v;
    amrex::Print() << "\n";

    // ---- 2. Find component indices for u, v, w ----
    int idx_u = -1, idx_v = -1, idx_w = -1;
    for (int n = 0; n < static_cast<int>(var_names.size()); ++n) {
        if (var_names[n] == "u") idx_u = n;
        if (var_names[n] == "v") idx_v = n;
        if (var_names[n] == "w") idx_w = n;
    }
    if (idx_u < 0 || idx_v < 0) {
        amrex::Print() << "[plt_wind_reader] Could not find 'u' and 'v' in plt variables\n";
        return false;
    }
    if (idx_w < 0) {
        amrex::Print() << "[plt_wind_reader] 'w' component not found; setting w=0\n";
    }

    const int nvar = static_cast<int>(var_names.size());
    const int n_pts_total = nx * ny * nz;
    if (n_pts_total <= 0) {
        amrex::Print() << "[plt_wind_reader] Empty domain\n";
        return false;
    }

    // Allocate host arrays
    pwd.x.assign(n_pts_total, 0.0);
    pwd.y.assign(n_pts_total, 0.0);
    pwd.z.assign(n_pts_total, 0.0);
    pwd.u.assign(n_pts_total, 0.0);
    pwd.v.assign(n_pts_total, 0.0);
    pwd.w.assign(n_pts_total, 0.0);

    // Cell sizes
    const double dx3 = (prob_hi[0] - prob_lo[0]) / nx;
    const double dy3 = (prob_hi[1] - prob_lo[1]) / ny;
    const double dz3 = (nz > 1) ? (prob_hi[2] - prob_lo[2]) / nz : 1.0;

    // ---- 3. Read the VisMF binary data from Level_0/Cell_* ----
    // The VisMF header is at level0_cell_path + "_H"
    // (WriteSingleLevelPlotfile writes "Level_0/Cell_H" and "Level_0/Cell_D_*")
    const std::string vismf_header = level0_cell_path + "_H";
    std::ifstream vhf(vismf_header);
    if (!vhf.is_open()) {
        amrex::Print() << "[plt_wind_reader] Cannot open VisMF header: "
                       << vismf_header << "\n";
        return false;
    }

    // ---- 4. Parse VisMF header to find FABs ----
    // AMReX VisMF-V1 header format:
    //   VisMF-V1
    //   <nFAB>
    //   <nComp>
    //   <nGhost>
    //   FAB (lo)(hi) ncomp   file:offset
    //   ...
    // We collect (file, offset, box) tuples.

    struct FabInfo {
        std::string filename;
        long offset;
        int lo[3], hi[3];
        int bytes_per_val;
    };
    std::vector<FabInfo> fabs;

    std::string vhdr_line;
    std::getline(vhf, vhdr_line);  // version string (e.g. "VisMF-V1")
    trim_string(vhdr_line);

    // Newer AMReX VisMF format uses a different layout
    bool is_v1 = (vhdr_line.find("VisMF-V1") != std::string::npos);

    if (is_v1) {
        int nfab = 0, ncomp_vm = 0, nghost_vm = 0;
        vhf >> nfab >> ncomp_vm >> nghost_vm;
        vhf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Read FAB descriptors  (box and file:offset)
        for (int n = 0; n < nfab; ++n) {
            std::string fab_line;
            std::getline(vhf, fab_line);
            trim_string(fab_line);

            // Parse box
            FabInfo fi;
            fi.bytes_per_val = 8;
            fi.offset = 0;

            // box is in format like "((0,0,0) (63,63,31) (0,0,0))"
            // file:offset at end
            // Try to extract file:offset  (appears after box substring)
            // Look for "File_0000:offset" or "Cell_D_00000:offset" patterns
            std::size_t colon = fab_line.rfind(':');
            if (colon != std::string::npos) {
                // Everything before the last space before the colon is box info;
                // after colon is the offset
                std::size_t space = fab_line.rfind(' ', colon);
                std::string file_str = (space != std::string::npos)
                    ? fab_line.substr(space + 1, colon - space - 1)
                    : fab_line.substr(0, colon);
                fi.filename = level0_cell_path.substr(0, level0_cell_path.rfind('/') + 1)
                              + file_str;
                fi.offset = std::stol(fab_line.substr(colon + 1));
            } else {
                // No file:offset info – skip
                continue;
            }

            // Parse box from the line
            if (!parse_fab_header(fab_line, fi.lo, fi.hi, fi.bytes_per_val)) {
                // Try the box from the file line itself
                fi.lo[0]=fi.lo[1]=fi.lo[2]=0;
                fi.hi[0]=nx-1; fi.hi[1]=ny-1; fi.hi[2]=nz-1;
            }
            fabs.push_back(fi);
        }
    } else {
        // Modern VisMF format: the header lists FABs differently.
        // AMReX v20+ VisMF header starts with the version then fab count.
        // Fall through to direct binary read below.
        vhf.close();
    }

    if (fabs.empty()) {
        // Fall back: try to read data directly from Cell_D_00000 with the
        // full domain box (single-FAB layout, common for serial runs).
        amrex::Print() << "[plt_wind_reader] VisMF header parsing yielded no FABs; "
                       << "attempting direct binary read of Level_0 data.\n";

        // In the direct mode the binary file contains:
        //   ASCII header line (FAB descriptor)\n
        //   binary data (ncomp doubles per cell, Fortran order within box)
        // for each FAB in sequence.  For a single-rank run there is one FAB
        // covering the full domain in Cell_D_00000.
        const std::vector<std::string> data_candidates = {
            level0_cell_path + "_D_00000",
            level0_cell_path + "_00000",
            level0_cell_path + "00000",
        };
        std::string data_file;
        for (const auto& c : data_candidates) {
            std::ifstream test(c);
            if (test.is_open()) { data_file = c; break; }
        }

        if (data_file.empty()) {
            amrex::Print() << "[plt_wind_reader] Cannot find binary data files in "
                           << plt_dir << "/Level_0/\n";
            return false;
        }

        FabInfo fi;
        fi.filename = data_file;
        fi.offset = 0;
        fi.lo[0]=fi.lo[1]=fi.lo[2]=0;
        fi.hi[0]=nx-1; fi.hi[1]=ny-1; fi.hi[2]=nz-1;
        fi.bytes_per_val = 8;
        fabs.push_back(fi);
    }

    // ---- 5. Read binary FAB data and fill pwd arrays ----
    for (const auto& fi : fabs) {
        const int box_nx = fi.hi[0] - fi.lo[0] + 1;
        const int box_ny = fi.hi[1] - fi.lo[1] + 1;
        const int box_nz = fi.hi[2] - fi.lo[2] + 1;
        const int box_n  = box_nx * box_ny * box_nz;

        std::ifstream df(fi.filename, std::ios::binary);
        if (!df.is_open()) {
            amrex::Print() << "[plt_wind_reader] Cannot open data file: "
                           << fi.filename << "\n";
            return false;
        }

        // Seek to FAB offset (in the binary file the data may be preceded by
        // an ASCII FAB header line that starts at the stored offset).
        df.seekg(fi.offset);
        if (!df.good()) {
            amrex::Print() << "[plt_wind_reader] Seek failed in " << fi.filename << "\n";
            return false;
        }

        // Skip the ASCII FAB header line (starts with "FAB")
        {
            char first_char = '\0';
            df.get(first_char);
            if (first_char == 'F') {
                // It's an ASCII header line; consume until newline
                std::string skip_line;
                std::getline(df, skip_line);
            } else {
                // Already at data — seek back one byte
                df.seekg(-1, std::ios::cur);
            }
        }

        // Read all components for this box
        const long data_bytes = static_cast<long>(nvar) * box_n * fi.bytes_per_val;
        std::vector<char> raw(data_bytes);
        df.read(raw.data(), data_bytes);
        if (df.gcount() != data_bytes) {
            amrex::Print() << "[plt_wind_reader] Short read in " << fi.filename
                           << " (expected " << data_bytes
                           << " bytes, got " << df.gcount() << ")\n";
            return false;
        }

        // Data layout: component-major, then z, y, x (Fortran column-major order)
        // i.e. offset = comp * box_n + k * box_nx * box_ny + j * box_nx + i
        auto get_val = [&](int comp, int bi, int bj, int bk) -> double {
            const int flat_idx = comp * box_n
                               + bk * box_nx * box_ny
                               + bj * box_nx
                               + bi;
            if (fi.bytes_per_val == 4) {
                float fv = 0.0f;
                std::memcpy(&fv, raw.data() + flat_idx * 4, 4);
                return static_cast<double>(fv);
            } else {
                double dv = 0.0;
                std::memcpy(&dv, raw.data() + flat_idx * 8, 8);
                return dv;
            }
        };

        // Fill pwd arrays for this FAB's box
        for (int bk = 0; bk < box_nz; ++bk) {
            const int gk = fi.lo[2] + bk;
            for (int bj = 0; bj < box_ny; ++bj) {
                const int gj = fi.lo[1] + bj;
                for (int bi = 0; bi < box_nx; ++bi) {
                    const int gi = fi.lo[0] + bi;
                    // Global flat index (x-fastest, then y, then z)
                    const int gflat = gk * nx * ny + gj * nx + gi;

                    pwd.x[gflat] = prob_lo[0] + (gi + 0.5) * dx3;
                    pwd.y[gflat] = prob_lo[1] + (gj + 0.5) * dy3;
                    pwd.z[gflat] = prob_lo[2] + (gk + 0.5) * dz3;
                    pwd.u[gflat] = static_cast<amrex::Real>(get_val(idx_u, bi, bj, bk));
                    pwd.v[gflat] = static_cast<amrex::Real>(get_val(idx_v, bi, bj, bk));
                    pwd.w[gflat] = (idx_w >= 0)
                        ? static_cast<amrex::Real>(get_val(idx_w, bi, bj, bk))
                        : amrex::Real(0.0);
                }
            }
        }
    }

    pwd.n_points = n_pts_total;

    // ---- 6. Build height-averaged 2-D wind ----
    // For each (i, j) column average u and v over all k levels.
    pwd.nx_2d = nx;
    pwd.ny_2d = ny;
    pwd.lo2d[0] = static_cast<amrex::Real>(prob_lo[0]);
    pwd.lo2d[1] = static_cast<amrex::Real>(prob_lo[1]);
    pwd.hi2d[0] = static_cast<amrex::Real>(prob_hi[0]);
    pwd.hi2d[1] = static_cast<amrex::Real>(prob_hi[1]);

    pwd.u2d.assign(nx * ny, amrex::Real(0.0));
    pwd.v2d.assign(nx * ny, amrex::Real(0.0));

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double u_sum = 0.0, v_sum = 0.0;
            for (int k = 0; k < nz; ++k) {
                const int gflat = k * nx * ny + j * nx + i;
                u_sum += static_cast<double>(pwd.u[gflat]);
                v_sum += static_cast<double>(pwd.v[gflat]);
            }
            pwd.u2d[j * nx + i] = static_cast<amrex::Real>(u_sum / nz);
            pwd.v2d[j * nx + i] = static_cast<amrex::Real>(v_sum / nz);
        }
    }

    // ---- 7. Upload to GPU ----
    auto upload = [](const std::vector<amrex::Real>& h,
                     amrex::Gpu::DeviceVector<amrex::Real>& d) {
        d.resize(h.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h.begin(), h.end(), d.begin());
    };

    upload(pwd.x,   pwd.d_x);
    upload(pwd.y,   pwd.d_y);
    upload(pwd.z,   pwd.d_z);
    upload(pwd.u,   pwd.d_u);
    upload(pwd.v,   pwd.d_v);
    upload(pwd.w,   pwd.d_w);
    upload(pwd.u2d, pwd.d_u2d);
    upload(pwd.v2d, pwd.d_v2d);

    amrex::Print() << "[plt_wind_reader] Loaded " << n_pts_total
                   << " 3-D wind points from " << plt_dir << "\n";
    amrex::Print() << "[plt_wind_reader] 2-D column-averaged wind grid: "
                   << nx << "×" << ny << "\n";

    pwd.valid = true;
    return true;
}

#endif // PLT_WIND_READER_H
