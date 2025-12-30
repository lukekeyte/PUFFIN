# PUFFIN: Python Utility For FUV Irradiated disk deNsities

`PUFFIN` is a parametric model for computing gas density structures in externally FUV-irradiated protoplanetary disks. It provides both 1D radial and 2D cylindrical models of photoevaporating disks, designed to be fast and accessible for applications in disk chemistry.

> [!IMPORTANT]
> **Note on this release:** This is the initial working version (v1.0) of `PUFFIN`. The physics implementation and model outputs are complete and validated. A refactored version with improved code structure (class-based design, enhanced documentation, and additional output options) is planned for release in January 2026. The current version is fully functional for research applications.

## Features

- **1D radial density profiles**: Fast parametric models with disk, wind, and transition components
- **2D cylindrical density profiles**: Full vertical structure with iterative hydrostatic equilibrium
- **Photoevaporative wind modeling**: Spherically diverging winds launched from disk surface
- **Smooth transitions**: Disk and wind are smoothly connected based on parameters calibrated against hundreds of hydrodynamical models
- **FRIED grid integration**: Automatic mass-loss rate interpolation from hydrodynamic simulations
- **Flexible parameterization**: Adjustable grid resolution, mass loss rates, and transition parameters

## Installation

`PUFFIN` is a standalone Python script. Simply clone or download the repository:
```bash
git clone https://github.com/lukekeyte/puffin.git
```

Ensure you have the required dependencies installed:
```bash
pip install numpy scipy matplotlib
```

## Repository Structure
```
puffin/
├── FRIEDV2_ALL_fPAH1p0_growth.dat   # FRIED grid lookup table
├── helpers.py                       # Helper functions
└── puffin.py                        # Main code (1D and 2D models)
```

## Dependencies

- numpy
- scipy
- matplotlib

## Quick Start

### 1D Model

For fast radial density profiles:
```python
import puffin

# Basic usage with default parameters
r_array, rho = puffin.DiskModel1D(
    m_star=0.5,      # stellar mass (M_sun)
    r_d=100,         # characteristic radius (AU)
    sigma_1au=100,   # surface density at 1 AU (g/cm^2)
    FFUV_G0=1000     # external FUV field (G0)
)

# Custom grid and parameters
r_array, rho = puffin.DiskModel1D(
    m_star=0.5,
    r_d=100,
    sigma_1au=100,
    FFUV_G0=1000,
    n_points=500,    # grid resolution
    gridsize=800,    # outer radius (AU)
    gamma=3.5,       # transition parameter
    p=0.3,           # transition parameter
    q=0.5            # transition parameter
)
```

### 2D Model

For full cylindrical density structures with vertical hydrostatic equilibrium:
```python
# Basic usage
r_array, z_array, rho_total = puffin.DiskModel2D(
    m_star=0.5,
    FFUV_G0=1000,
    r_d=100,
    sigma_1au=100
)

# Custom parameters
r_array, z_array, rho_total = puffin.DiskModel2D(
    m_star=0.5,
    FFUV_G0=1000,
    r_d=100,
    sigma_1au=100,
    n_points=1000,   
    N_ITER=20, 
    gamma=3.5,    
    p=0.2,
    q=0.4
)
```


## Citation

If you use `PUFFIN` in your research, please cite:

Keyte & Haworth (submitted) - 'A parametric model for externally irradiated protoplanetary disks with photoevaporative winds'

## Author

**Luke Keyte**  
Postdoctoral Researcher  
Queen Mary University of London  
l.keyte@qmul.ac.uk

## Contributing

Contributions, bug reports, and feature requests are welcome! Please feel free to open an issue or contact the author directly.

## License

This project is licensed under the MIT License.