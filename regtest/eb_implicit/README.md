# EB Implicit Function Test

## Purpose
Tests embedded boundary (EB) capabilities using implicit function representations for initial conditions.

## Features
- EB implicit function for initial condition
- Support for multiple geometry types (plane, cylinder, sphere, ellipsoid)
- Signed distance function from implicit representation
- Compatible with level-set advection
- This example uses an ellipsoid geometry

## Expected Runtime
~1 minute

## Key Parameters
- Domain: 1.0³ unit cube
- Grid: 64³ cells
- EB type: ellipsoid
- Ellipsoid center: (0.3, 0.5, 0.5)
- Semi-axes: rx=0.20, ry=0.15, rz=0.10
- Velocity: (0.25, 0.1, 0.0) m/s
- Time steps: 100

## Supported EB Types

### 1. Plane
```
source_type = eb
eb_type = plane
eb_param1 = 1.0    # normal_x
eb_param2 = 0.0    # normal_y
eb_param3 = 0.0    # normal_z
eb_param4 = -0.5   # offset d (plane: nx*x + ny*y + nz*z + d = 0)
```

### 2. Cylinder
```
source_type = eb
eb_type = cylinder
eb_param1 = 0.5    # center_x
eb_param2 = 0.5    # center_y
eb_param3 = 0.0    # (unused for cylinder, axis along z)
eb_param4 = 0.2    # radius
```

### 3. Sphere
```
source_type = eb
eb_type = sphere
eb_param1 = 0.5    # center_x
eb_param2 = 0.5    # center_y
eb_param3 = 0.5    # center_z
eb_param4 = 0.25   # radius
```

### 4. Ellipsoid
```
source_type = eb
eb_type = ellipsoid
eb_param1 = 0.5    # center_x
eb_param2 = 0.5    # center_y
eb_param3 = 0.5    # center_z
eb_param4 = 0.25   # radius_x
eb_param5 = 0.15   # radius_y
eb_param6 = 0.10   # radius_z
```

## Expected Behavior
The EB implicit function creates a signed distance field that can be used as the initial condition for the level-set method. The zero level set represents the geometry boundary (e.g., ellipsoid surface in this example).

The fire front should advect from the ellipsoidal initial condition in the direction of the velocity field.

## Notes
- The EB implicit function provides an exact signed distance for simple geometries (plane, sphere, cylinder)
- For ellipsoids, an approximate signed distance is used (same as the ellipse source type)
- This capability can be extended to support more complex implicit function representations
- Unlike AMReX's full EB support, this is a simplified approach using analytical implicit functions

## Visualization
Use ParaView or VisIt to visualize the phi field. The zero level set should show the specified geometry that advects with the flow.
