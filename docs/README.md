# Documentation

This directory contains the Sphinx documentation for the Wildfire Level-Set Solver.

## Building Documentation Locally

### Prerequisites

Install the required Python packages:

```bash
pip install -r requirements.txt
```

### Build HTML

```bash
cd docs
make html
```

The generated HTML files will be in `_build/html/`. Open `_build/html/index.html` in a web browser.

### Build PDF

```bash
cd docs
make latexpdf
```

The generated PDF will be in `_build/latex/`.

### Clean Build

```bash
cd docs
make clean
```

## Building with CMake

You can also build the documentation using CMake:

```bash
cmake -S . -B build -DLEVELSET_BUILD_DOCS=ON
cmake --build build --target docs
```

The generated HTML will be in `build/docs/index.html`.

## Online Documentation

The documentation is automatically built and deployed to GitHub Pages when changes are pushed to the main branch.

Visit: https://hgopalan.github.io/wildfire_levelset/

## Documentation Structure

* `index.rst` - Main documentation entry point
* `overview.rst` - Project overview and summary
* `mathematical_models.rst` - Mathematical equations and models
* `code_structure.rst` - Code organization and architecture
* `building.rst` - Build instructions
* `usage.rst` - Usage guide and input parameters
* `api_reference.rst` - API reference for functions and data structures
* `conf.py` - Sphinx configuration
* `requirements.txt` - Python dependencies for building docs
