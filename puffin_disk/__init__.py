"""
PUFFIN: Python Utility For FUV Irradiated disk deNsities
=========================================================

A parametric model for calculating the density structure of externally 
FUV-irradiated protoplanetary disks undergoing photoevaporation.

This package provides tools to compute 1D and 2D gas density profiles for 
protoplanetary disks subject to external far-ultraviolet (FUV) radiation fields,
with applications to externally irradiated environments like the Orion Nebula Cluster.

Main Functions
--------------
DiskModel1D : Compute 1D midplane density structure
DiskModel2D : Compute 2D (r,z) density structure with vertical hydrostatic equilibrium

References
----------
If using this code, please cite:
  - Keyte & Haworth (2026) - PUFFIN overview paper
  - Haworth et al. (2018, 2023) - FRIED grid mass loss rates
"""

__version__ = "0.1.0"
__author__  = "Luke Keyte"

from .puffin import DiskModel1D, DiskModel2D

from .helpers import (
    plot_density,           # Visualize 2D density structures
    plot_all,               # Comprehensive plotting of model outputs
    calculate_disk_mass,    # Calculate total disk mass from model
    interpolate_mdot,       # Interpolate mass loss rates from FRIED grid
    get_stellar_properties  # Get stellar parameters from mass (Eker et al. 2018)
)

__all__ = [
    # Core modeling functions
    'DiskModel1D',
    'DiskModel2D',
    # Visualization tools
    'plot_density',
    'plot_all',
    # Analysis utilities
    'calculate_disk_mass',
    'interpolate_mdot',
    'get_stellar_properties',
]
