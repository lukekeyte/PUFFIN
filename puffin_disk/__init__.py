"""
PUFFIN: Python Utility For FUV Irradiated disk deNsities
=========================================================
A parametric model for externally FUV-irradiated protoplanetary disks.

Author: Luke Keyte
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Luke Keyte"

from .puffin import DiskModel1D, DiskModel2D

from .helpers import (
    smooth_temperature_profile,
    eker_mlr,
    get_stellar_properties,
    planck_wavelength,
    fuv_fraction,
    calculate_disk_mass,
    interpolate_mdot,
    smooth_cosine,
)

__all__ = [
    "DiskModel1D",
    "DiskModel2D",
    "smooth_temperature_profile",
    "eker_mlr",
    "get_stellar_properties",
    "planck_wavelength",
    "fuv_fraction",
    "calculate_disk_mass",
    "interpolate_mdot",
    "smooth_cosine",
]
