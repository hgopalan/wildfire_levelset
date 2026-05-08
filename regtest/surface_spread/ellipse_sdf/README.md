# Elliptical SDF Test

## Purpose
Tests elliptical initial condition with signed distance function (SDF).

## Features
- Elliptical initial condition
- Approximate SDF for ellipse
- Constant velocity field
- Periodic reinitialization
- Basic advection of elliptical shape

## Expected Runtime
~1 minute

## Key Parameters
- Domain: 1.0³ unit cube
- Grid: 64³ cells
- Ellipse center: (0.4, 0.5, 0.5)
- Semi-axes: rx=0.25, ry=0.15, rz=0.10
- Velocity: (0.2, 0.15, 0.0) m/s
- Time steps: 100
- Reinit every 20 steps

## Expected Behavior
The elliptical fire front should advect in the direction of the velocity field while maintaining its elliptical shape. The signed distance property should be preserved through periodic reinitialization.

## Visualization
Use ParaView or VisIt to visualize the phi field. The zero level set should show an ellipse that moves and potentially deforms slightly due to advection.
