# PUFFIN: Python Utility For FUV Irradiated disk deNsities

`PUFFIN` is a parametric model for computing gas density structures in externally FUV-irradiated protoplanetary disks. It provides both 1D radial and 2D cylindrical models of photoevaporating disks, designed to be fast and accessible for applications in disk chemistry.

[![Documentation](https://readthedocs.org/projects/puffin/badge/?version=latest)](https://puffin.readthedocs.io)
[![PyPI version](https://badge.fury.io/py/puffin-disk.svg)](https://badge.fury.io/py/puffin-disk)

## Features

- **1D radial density profiles**: Fast parametric models with disk, wind, and transition components
- **2D cylindrical density profiles**: Full vertical structure with iterative hydrostatic equilibrium
- **FRIED grid integration**: Automatic mass-loss rate interpolation from hydrodynamic simulations
- **Flexible parameterization**: Adjustable grid resolution, mass loss rates, and model parameters

## Installation

Install `puffin_disk` directly from PyPI using pip:

```bash
pip install puffin_disk
```

This will automatically install all required dependencies (numpy, scipy, matplotlib).

### Development Installation

To install the latest development version from GitHub:

```bash
pip install git+https://github.com/lukekeyte/PUFFIN.git
```

### Requirements

- Python 3.8 or higher
- numpy >= 1.20.0
- scipy >= 1.7.0
- matplotlib >= 3.3.0

## Quick Start

### 1D Model

For fast radial density profiles:

```python
import puffin_disk

# Basic usage with default parameters
r_array, rho = puffin_disk.DiskModel1D(
    m_star=0.5,      # stellar mass (M_sun)
    r_d=100,         # characteristic radius (AU)
    sigma_1au=100,   # surface density at 1 AU (g/cm^2)
    FFUV_G0=1000     # external FUV field (G0)
)

# Custom grid and parameters
r_array, rho = puffin_disk.DiskModel1D(
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
import puffin_disk

# Basic usage
r_array, z_array, rho_total = puffin_disk.DiskModel2D(
    m_star=0.5,
    FFUV_G0=1000,
    r_d=100,
    sigma_1au=100
)

# Custom parameters
r_array, z_array, rho_total = puffin_disk.DiskModel2D(
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

# Visualize the result
puffin_disk.plot_density(rho_total, r_array, z_array)
```


## Documentation

Full documentation including detailed API reference, tutorials, and examples is available at [puffin.readthedocs.io](https://puffin.readthedocs.io).

## Physical Model

The density structure consists of three components:

1. **Disk**: Hydrostatic disk with power-law surface density (Σ ∝ r⁻¹) and exponential outer truncation at the gravitational radius r_d
2. **Wind**: Spherical photoevaporative outflow (ρ ∝ r⁻²) launched from the τ=1 FUV surface, with density set by mass loss rate
3. **Transition**: Smooth exponential taper and plateau region blending disk to wind

The 2D model iteratively solves for vertical hydrostatic equilibrium with temperature-dependent scale heights, accounting for FUV heating in the photodissociation region (PDR).

## Citation

If you use `PUFFIN` in your research, please cite:

**Keyte & Haworth (2026)** - *A parametric model for externally irradiated protoplanetary disks with photoevaporative winds*

The FRIED grid mass loss rates are from:

**Haworth et al. (2018)** - *The FRIED grid of mass-loss rates for externally irradiated
protoplanetary discs* - MNRAS, 481, 452  
**Haworth et al. (2023)** - *FRIED v2: a new grid of mass-loss rates for externally irradiated protoplanetary discs* - MNRAS, 526, 4315

## Contributing

Contributions, bug reports, and feature requests are welcome! Please feel free to:

- Open an issue on [GitHub](https://github.com/lukekeyte/PUFFIN/issues)
- Submit a pull request
- Contact the author directly at l.keyte@qmul.ac.uk

## Support

If you encounter any problems or have questions:

- Check the [documentation](https://puffin.readthedocs.io)
- Search existing [GitHub issues](https://github.com/lukekeyte/PUFFIN/issues)
- Open a new issue with a minimal reproducible example

## Author

**Luke Keyte**  
Postdoctoral Researcher  
Queen Mary University of London  
l.keyte@qmul.ac.uk

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.