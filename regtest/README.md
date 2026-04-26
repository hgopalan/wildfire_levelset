# Regression Test Suite

This directory contains regression tests for the wildfire level-set solver. Each test focuses on specific capabilities and features of the code.

## Test Directory Structure

```
regtest/
├── basic_levelset/       # Basic level-set advection with constant velocity
├── farsite_ellipse/      # FARSITE elliptical fire expansion (Richards 1990)
├── rothermel_fuel/       # Rothermel model with different fuel types
├── terrain_wind/         # External terrain (Gaussian hill) and wind field
├── anderson_lw/          # Anderson (1983) dynamic L/W ratio
├── reinitialization/     # Level-set reinitialization testing
├── ellipse_sdf/          # Elliptical SDF initial condition
├── eb_implicit/          # EB implicit function initial condition
└── 3d_sphere/           # Full 3D fire spread simulation
```

## Quick Start

### Running a Single Test

From the repository root:
```bash
./build/levelset regtest/<test_name>/inputs.i
```

For tests with data files (e.g., terrain_wind), run from the test directory:
```bash
cd regtest/terrain_wind
../../build/levelset inputs.i
```

### Running All Tests

You can run all tests with a simple script:
```bash
# From repository root
for test in regtest/*/; do
    echo "Running test: $test"
    if [ -f "$test/inputs.i" ]; then
        cd "$test"
        ../../build/levelset inputs.i > output.log 2>&1
        cd ../..
    fi
done
```

## Test Descriptions

### 1. basic_levelset
**Purpose**: Tests fundamental level-set advection capabilities.

**Features**:
- Spherical initial condition
- Constant velocity field
- Periodic reinitialization
- Basic advection verification

**Expected Runtime**: ~1 minute

**Key Parameters**:
- Domain: 1.0³ unit cube
- Grid: 64³ cells
- Time steps: 100

---

### 2. farsite_ellipse
**Purpose**: Tests FARSITE elliptical fire expansion model with fixed coefficients.

**Features**:
- Richards' (1990) ellipse model
- Fixed head/flank/backing fire coefficients
- Box (line) fire initial condition
- Wind-driven elliptical spread

**Expected Runtime**: ~1 minute

**Key Parameters**:
- FARSITE enabled with fixed L/W = 3.0
- Richards coefficients: a=1.0, b=0.4, c=0.2
- Constant wind field

---

### 3. rothermel_fuel
**Purpose**: Tests fire spread with different fuel models from NFFL database.

**Features**:
- Rothermel fire spread model
- NFFL fuel model database (FM1-FM13)
- Fuel moisture effects
- Wind-driven spread

**Expected Runtime**: ~1 minute

**Key Parameters**:
- Fuel: FM1 (Short Grass)
- Fuel moisture: 6%
- Moderate wind conditions

**Variants**: Modify `rothermel.fuel_model` to test FM1-FM13.

---

### 4. terrain_wind ⭐
**Purpose**: Tests external terrain and spatially-varying wind field (Gaussian hill).

**Features**:
- Gaussian hill terrain (100 m height)
- Spatially-varying wind field with speedup over crest
- Terrain slope calculation from elevation data
- Unstructured data interpolation (IDW)
- Anderson L/W ratio with terrain effects

**Expected Runtime**: ~2-3 minutes

**Key Parameters**:
- Domain: 1000 m × 1000 m
- Grid: 100 × 100 cells
- Terrain: Gaussian hill with σ=150m
- Wind: Base 5 m/s with up to 50% speedup

**Build Requirements**: Requires 2D build (`-DLEVELSET_DIM_2D=ON`)

**Data Files**:
- `gaussian_hill_terrain.csv`: 2500 elevation points
- `gaussian_hill_wind.csv`: 2500 wind velocity points

---

### 5. anderson_lw
**Purpose**: Tests dynamic L/W ratio calculation based on wind speed (Anderson 1983).

**Features**:
- Anderson (1983) L/W formula
- Dynamic ellipse elongation
- Wind speed dependent behavior
- Automatic coefficient calculation

**Expected Runtime**: ~1 minute

**Key Parameters**:
- Wind: ~10 mph (expected L/W ≈ 2.5)
- Dynamic coefficient computation

**Formula**: L/W = 0.936·exp(0.2566·U) + 0.461·exp(-0.1548·U) - 0.397

---

### 6. reinitialization
**Purpose**: Tests level-set reinitialization to maintain signed distance property.

**Features**:
- Aggressive reinitialization (every 5 steps)
- Signed distance preservation
- Interface sharpness
- Comparison with less frequent reinitialization

**Expected Runtime**: ~2 minutes

**Key Parameters**:
- Reinit frequency: Every 5 steps
- Reinit iterations: 30
- Shape preservation analysis

---

### 7. 3d_sphere
**Purpose**: Tests full 3D fire spread simulation.

**Features**:
- 3D level-set advection
- 3D FARSITE elliptical expansion
- 3D Rothermel model
- 3D terrain effects
- 3D reinitialization

**Expected Runtime**: ~5-10 minutes

**Key Parameters**:
- Domain: 1.0³ unit cube
- Grid: 48³ cells
- 3D wind vector
- Upslope terrain

**Build Requirements**: Requires 3D build (default)

**Visualization**: Use ParaView/VisIt for 3D isosurfaces

---

### 8. ellipse_sdf
**Purpose**: Tests elliptical initial condition with signed distance function.

**Features**:
- Elliptical SDF initial condition
- Approximate signed distance for ellipse
- Advection of elliptical shape
- Multi-axis radii specification

**Expected Runtime**: ~1 minute

**Key Parameters**:
- Domain: 1.0³ unit cube
- Grid: 64³ cells
- Ellipse center: (0.4, 0.5, 0.5)
- Semi-axes: rx=0.25, ry=0.15, rz=0.10
- Constant velocity field

**Formula**: Approximate SDF for ellipsoid in normalized coordinates

---

### 9. eb_implicit
**Purpose**: Tests embedded boundary capabilities using implicit function representations.

**Features**:
- EB implicit function for initial conditions
- Support for plane, cylinder, sphere, ellipsoid geometries
- Exact signed distance for simple geometries
- Extensible to complex implicit functions

**Expected Runtime**: ~1 minute

**Key Parameters**:
- Domain: 1.0³ unit cube
- Grid: 64³ cells
- EB type: ellipsoid
- Parameterized geometry specification

**Supported EB Types**:
- `plane`: Defined by normal vector and offset
- `cylinder`: Circular cylinder along z-axis
- `sphere`: Spherical geometry
- `ellipsoid`: Ellipsoidal geometry

---

## Testing Checklist

Use this checklist to verify all capabilities:

- [ ] **Basic Advection**: `basic_levelset` runs successfully
- [ ] **FARSITE**: `farsite_ellipse` produces elliptical spread
- [ ] **Fuel Models**: `rothermel_fuel` correctly applies fuel properties
- [ ] **Terrain**: `terrain_wind` handles external elevation data
- [ ] **Wind Fields**: `terrain_wind` interpolates spatially-varying wind
- [ ] **Anderson L/W**: `anderson_lw` computes dynamic coefficients
- [ ] **Reinitialization**: `reinitialization` maintains |∇φ|=1
- [ ] **Elliptical SDF**: `ellipse_sdf` creates elliptical initial conditions
- [ ] **EB Capabilities**: `eb_implicit` uses implicit function geometries
- [ ] **3D Capability**: `3d_sphere` runs in 3D mode

## Build Configurations

### Default 3D Build
```bash
cmake -S . -B build
cmake --build build -j
```

**Compatible tests**: `basic_levelset`, `farsite_ellipse`, `rothermel_fuel`, `anderson_lw`, `reinitialization`, `ellipse_sdf`, `eb_implicit`, `3d_sphere`

### 2D Build
```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
```

**Required for**: `terrain_wind` (uses 2D wind/terrain data files)

**Compatible tests**: All except `3d_sphere`

### Running Tests with CMake/CTest

After building, you can run regression tests using CMake's testing framework:

```bash
# Run all regression tests
cd build
ctest -L regtest --output-on-failure

# Or use the custom target
make regtest

# Run a specific test
ctest -R ellipse_sdf --output-on-failure
```

## Output

Each test produces:
- **Plotfiles**: `plt####` directories containing field data
- **Console output**: Simulation progress and diagnostics

### Plotfile Contents
- `phi`: Level-set function (signed distance or indicator)
- `velx`, `vely`, `velz`: Velocity field components
- `farsite_dx`, `farsite_dy`, `farsite_dz`: FARSITE displacements (if enabled)

### Visualization
Use ParaView, VisIt, or AMReX-native tools:
```bash
# With ParaView
paraview plt*
```

## Test Validation

### Automated Checks (Future)
- Verify plotfiles are generated
- Check conservation properties
- Compare against reference solutions
- Verify ellipse aspect ratios

### Manual Verification
1. **basic_levelset**: Sphere should advect without distortion
2. **farsite_ellipse**: Fire should form ellipse with L/W ≈ 3.0
3. **rothermel_fuel**: Different fuels show different spread rates
4. **terrain_wind**: Fire accelerates upslope and over crest
5. **anderson_lw**: L/W ratio matches Anderson formula
6. **reinitialization**: |∇φ| remains close to 1.0
7. **3d_sphere**: 3D ellipsoid formation

## Adding New Tests

To add a new regression test:

1. Create a new directory: `regtest/<test_name>/`
2. Add `inputs.i` with test parameters
3. Add `README.md` documenting the test
4. Add any required data files
5. Update this main README.md

### Test Template

```
# Test Name
# Tests: Brief description

# Grid & domain
n_cell = 64
prob_lo_x = 0.0
prob_hi_x = 1.0

# Time & output
nsteps = 100
plot_int = 10

# Add test-specific parameters...
```

## References

1. **Richards, G.D. (1990)**. "An elliptical growth model of forest fire fronts and its numerical solution." International Journal of Numerical Methods in Engineering.

2. **Anderson, H.E. (1983)**. "Predicting wind-driven wild land fire size and shape." Research Paper INT-305, USDA Forest Service.

3. **Rothermel, R.C. (1972)**. "A mathematical model for predicting fire spread in wildland fuels." Research Paper INT-115, USDA Forest Service.

4. **Osher, S. & Fedkiw, R. (2003)**. "Level Set Methods and Dynamic Implicit Surfaces." Springer.

## Troubleshooting

### Test fails to run
- Check build configuration (2D vs 3D)
- Verify input file syntax
- Check for required data files

### No plotfiles generated
- Verify `plot_int` is set
- Check `nsteps` > 0
- Ensure sufficient disk space

### Unexpected results
- Review README for expected behavior
- Check parameter values
- Verify initial conditions
- Compare with similar tests

## Contact

For questions or issues with tests, please open an issue on the GitHub repository.
