# Fuel Model Database Tests

This directory contains test cases for the Rothermel fuel model database.

## Running Tests

### Test 1: Short Grass (FM1)
```bash
../../build/levelset nsteps=10 plot_int=10 rothermel.fuel_model=FM1
```

### Test 2: Chaparral (FM4) - Default
```bash
../../build/levelset nsteps=10 plot_int=10 rothermel.fuel_model=FM4
```

### Test 3: Using Alias (CHAPARRAL)
```bash
../../build/levelset nsteps=10 plot_int=10 rothermel.fuel_model=CHAPARRAL
```

### Test 4: Heavy Slash (FM13)
```bash
../../build/levelset nsteps=10 plot_int=10 rothermel.fuel_model=FM13
```

### Test 5: Fuel Model with Custom Moisture
```bash
../../build/levelset nsteps=10 plot_int=10 rothermel.fuel_model=FM1 rothermel.M_f=0.15
```

### Test 6: Numeric Alias
```bash
../../build/levelset nsteps=10 plot_int=10 rothermel.fuel_model=5
```

## Expected Behavior

- When a valid fuel model is specified, the program should display the fuel model information
- When an invalid fuel model is specified, a warning should be displayed with available options
- Individual parameters can override fuel model values
- Aliases (numeric or descriptive) should work the same as FM codes

## Available Fuel Models

See the fuel_database.H file for the complete Anderson 13 fuel model database:

- **FM1**: Short Grass (1 ft)
- **FM2**: Timber (Grass and Understory)  
- **FM3**: Tall Grass (2.5 ft)
- **FM4**: Chaparral (6 ft) - Default
- **FM5**: Brush (2 ft)
- **FM6**: Dormant Brush, Hardwood Slash
- **FM7**: Southern Rough
- **FM8**: Closed Timber Litter
- **FM9**: Hardwood Litter
- **FM10**: Timber (Litter and Understory)
- **FM11**: Light Logging Slash
- **FM12**: Medium Logging Slash
- **FM13**: Heavy Logging Slash

## Aliases

Common aliases that work for fuel models:
- Numeric: `1`, `2`, `3`, ..., `13`
- Descriptive: `SHORT_GRASS`, `GRASS`, `TALL_GRASS`, `CHAPARRAL`, `BRUSH`, etc.
