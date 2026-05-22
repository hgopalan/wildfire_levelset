# Fire Whirl/Vorticity-Driven Spotting - Implementation Status

## Summary

**Status:** Fire whirl and vorticity models are **implemented but NOT integrated with spotting**.

The codebase contains fully functional fire whirl and vorticity calculation models, but these are currently independent diagnostic outputs and do not affect spotting distance calculations.

---

## What's Currently Implemented

### 1. Vorticity Model (`src/vorticity_model.H`)
- ✅ Computes vertical vorticity component: ω_z = ∂v/∂x − ∂u/∂y [s⁻¹]
- ✅ Cell-centered finite difference calculation
- ✅ GPU-ready kernel (CUDA/HIP/SYCL compatible)
- ✅ Written to plotfiles as `vorticity_z` field
- **Purpose:** Diagnostic field for identifying fire whirls, horizontal roll vortices, and wind shear

### 2. Weise & Biging (1996) Fire Whirl Model (`src/weise_biging_whirl.H`)
- ✅ Full fire whirl characterization model
- ✅ Computes 6 output fields:
  - `weise_flame_height` - vertical flame height [m]
  - `weise_flame_tilt` - flame tilt angle from vertical [rad]
  - `weise_whirl_height` - fire whirl column height [m]
  - `weise_whirl_radius` - fire whirl core radius [m]
  - `weise_angular_velocity` - angular velocity [rad/s]
  - `weise_max_tang_vel` - maximum tangential velocity [m/s]
- ✅ Based on modified Froude number and Rankine vortex scaling
- ✅ Enabled via `weise_biging.enable = 1` in input file
- ✅ All outputs written to plotfiles
- **Purpose:** Fire whirl geometry and kinematics prediction

### 3. Spotting Models (`src/firebrand_spotting.H`, `src/albini_spotting.H`)
- ✅ Albini (1983) trajectory-based spotting
- ✅ Firebrand trajectory calculation with 3-D wind
- ✅ Stochastic ignition probability
- ✅ Torching-tree spotting (Albini 1979)
- ❌ **NO connection to vorticity or fire whirl fields**

---

## What's Missing

### Integration Gap

**The fire whirl and vorticity fields are computed but NOT used to modify spotting behavior.**

Current spotting models use:
- Wind speed at multiple heights (3-D wind field)
- Fireline intensity
- Fuel moisture
- Tree height (for torching)

They **do NOT use**:
- Vertical vorticity (ω_z)
- Fire whirl characteristics (radius, angular velocity, tangential velocity)
- Flame tilt angle

---

## Why This Matters (FARSITE Context)

In FARSITE and other advanced fire models, fire whirls can significantly enhance spotting by:

1. **Lofting firebrands higher** - Whirl updrafts carry embers to greater heights
2. **Increasing horizontal transport** - Whirl tangential velocities add to ambient wind
3. **Extending spotting distance** - Combined effect can increase spotting distance by 2-5×

### Physics-Based Connection

Fire whirls create:
- Strong updrafts (5-20 m/s vertical velocity)
- Concentrated vorticity (angular velocity 1-10 rad/s)
- Enhanced lofting of firebrands

The spotting distance enhancement could be modeled as:
```
d_spot_whirl = d_spot_base × (1 + α × (v_θ/U) × (Ω × r_w / U))
```

Where:
- `d_spot_base` = baseline Albini spotting distance
- `α` = empirical enhancement factor (~ 0.5-2.0)
- `v_θ` = whirl tangential velocity
- `Ω` = angular velocity
- `r_w` = whirl radius
- `U` = ambient wind speed

---

## Implementation Recommendation

### Difficulty: **MODERATE to HARD**

To integrate vorticity/whirl with spotting:

#### Approach 1: Simple Enhancement Factor (MODERATE)
1. In `firebrand_spotting.H`, check if `weise_biging.enable == 1`
2. If enabled, read `weise_max_tang_vel` and `weise_whirl_radius` at source cell
3. Compute whirl enhancement factor
4. Multiply base spotting distance by enhancement factor

**Changes needed:**
- Modify `compute_firebrand_spotting()` to accept optional `weise_data` MultiFab
- Add whirl enhancement calculation in spotting kernel
- Update documentation

**Pros:** Simple, physically motivated
**Cons:** Lacks detailed trajectory physics

#### Approach 2: Full Trajectory Integration (HARD)
1. Modify Albini trajectory solver to include whirl vertical velocity
2. Add rotational component to firebrand motion
3. Track firebrands through whirl core

**Changes needed:**
- Significant rewrite of trajectory integration
- 3-D vorticity field (not just vertical component)
- Empirical validation against whirl spotting data

**Pros:** Physically accurate
**Cons:** Complex, requires validation data

---

## Recommended Action

### Short-term (EASY - for documentation only)
✅ Add note to documentation that vorticity/whirl are diagnostic only
✅ Recommend users analyze whirl fields post-simulation to identify high-risk zones

### Medium-term (MODERATE - ~1-2 days work)
- Implement Approach 1: Simple enhancement factor
- Add parameter `weise_biging.enhance_spotting` (default: 0)
- Test with high-wind scenarios where whirls are likely

### Long-term (HARD - research project)
- Full trajectory coupling with whirl dynamics
- Validate against field observations (if available)
- Publish methodology

---

## References

- Weise, D.R. & Biging, G.S. (1996). "Effects of wind velocity and slope on flame properties." *Canadian Journal of Forest Research*, 26(10), 1849–1858.
- Albini, F.A. (1983). "Potential spotting distance from wind-driven surface fires." USDA Forest Service Research Paper INT-309.
- Cunningham, P. & Linn, R.R. (2007). "Numerical simulations of grass fires using a coupled atmosphere-fire model." *Journal of Geophysical Research*, 112, D05108.

---

## Conclusion

**Fire whirl/vorticity models are implemented and functional** but serve as **diagnostic outputs only**. They do **not affect spotting calculations**. This is consistent with most operational fire models (including FARSITE), which treat vorticity as a post-analysis diagnostic rather than a predictive component.

For enhanced realism in extreme fire scenarios, integrating vorticity effects into spotting would be valuable but requires moderate to significant development effort.
