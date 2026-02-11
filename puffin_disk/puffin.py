"""
PUFFIN: Python Utility For FUV Irradiated disk deNsities
========================================================
A parametric model for calculating the density structure of externally 
FUV-irradiated protoplanetary disks undergoing photoevaporation.

This module provides functions to compute 1D and 2D gas density profiles for 
protoplanetary disks subject to external far-ultraviolet (FUV) radiation fields.

Main Functions
--------------
DiskModel1D : 
    Compute 1D midplane density structure
DiskModel2D : 
    Compute 2D (r,z) density structure with vertical hydrostatic equilibrium

Physical Model
--------------
The density structure consists of three components:
  1. Disk: Hydrostatic disk with power-law surface density Sigma proportional to r^(-1) and 
     exponential outer truncation at the gravitational radius r_d
  2. Wind: Spherical photoevaporative outflow (rho proportional to r^(-2)) launched from 
     the tau=1 FUV surface, with density set by mass loss rate
  3. Transition: Smooth exponential taper and 'plateau' region blending disk to wind

The 2D model iteratively solves for vertical hydrostatic equilibrium with 
temperature-dependent scale heights, accounting for FUV heating in the 
photodissociation region (PDR).

Dependencies
------------
numpy : Array operations and mathematical functions
scipy.spatial.cKDTree : Nearest-neighbor searches for wind geometry
scipy.interpolate.griddata : Spatial interpolation
helpers : Custom utilities for I/O, stellar properties, and interpolation

Mass Loss Rates
---------------
If not explicitly provided, mass loss rates are interpolated from the FRIED 
grid (Haworth et al. 2018, 2023) using the helpers.interpolate_mdot function, which 
queries a pre-computed lookup table spanning:
  - Stellar mass: 0.3 - 3.0 Msun
  - Disk radius: 10 - 150 AU
  - Surface density: 10 - 10^4 g/cm^2
  - FUV field: 100 - 10^5 G0

References
----------
If using this code, please cite:
  - Keyte & Haworth (2026) - PUFFIN overview paper
  - Haworth et al. (2018, 2023) - FRIED grid for mass-loss rates

Author: Luke Keyte
Institution: Queen Mary University of London, Astronomy Unit
Contact: l.keyte@qmul.ac.uk
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import platform
from . import helpers


class DiskModel1D:
    """
    Calculate 1D radial density structure for externally FUV-irradiated protoplanetary disks.
    
    This function computes the gas density profile for a disk undergoing external 
    photoevaporation. The model includes three components: the power-law disk structure,
    an outflowing wind, and a transitional region that smoothly connects them.
    
    Parameters
    ----------
    m_star : float
        Stellar mass in solar masses
    r_d : float
        Disk outer radius in AU
    sigma_1au : float
        Surface density at 1 AU in g/cm^2
    FFUV_G0 : float
        External FUV field strength in Habing units (G0)
    n_points : int, optional
        Number of radial grid points (default: 200)
    gridsize : float, optional
        Outer edge of radial grid in AU (default: 8 * r_d)
    m_dot : float, optional
        Mass loss rate in M_sun/yr (default: computed from interpolation of FRIED grid)
    gamma : float, optional
        Exponential taper parameter (default: computed from prescription)
    p : float, optional
        Plateau parameter (default: 0.2)
    q : float, optional
        Plateau parameter (default: 0.4)
    
    Attributes
    ----------
    radius : ndarray
        Radial grid in AU
    density : ndarray
        Total midplane density in g/cm^3
    temperature : ndarray
        Disk temperature in K
    rho_disk : ndarray
        Disk component density in g/cm^3
    rho_wind : ndarray
        Wind component density in g/cm^3
    rho_plateau : ndarray
        Plateau component density in g/cm^3
    """
    
    # Physical constants
    MU_DISK = 2.3
    MU_PDR  = 1.3
    BIG_G   = 6.67259e-8      # cm^3 g^-1 s^-2
    K_B     = 1.380626e-16    # erg K^-1
    M_SUN   = 1.9891e33       # g
    AU_CM   = 1.495979e13     # cm
    M_P     = 1.6726219e-24    # g
    YR_SEC  = 3.15576e7       # s
    
    
    def __init__(self, m_star, r_d, sigma_1au, FFUV_G0, n_points=200, gridsize=None, 
                 m_dot=None, gamma=None, p=0.2, q=0.4, verbose=True):

        # Initialise model
        self.verbose = verbose
        self._validate_inputs(m_star, r_d, sigma_1au, FFUV_G0, n_points, m_dot, gamma, p, q, gridsize)

        # Store input parameters
        self.m_star = m_star
        self.r_d = r_d
        self.sigma_1au = sigma_1au
        self.FFUV_G0 = FFUV_G0
        self.n_points = n_points
        self.gridsize = gridsize if gridsize is not None else r_d * 8
        self.m_dot = m_dot
        self.gamma = gamma
        self.p = p
        self.q = q  
        
        # Initialize result arrays (computed when run() is called)
        self.radius = None
        self.density = None
        self.temperature = None
        self.rho_disk = None
        self.rho_wind = None
        self.rho_plateau = None
        
        
        
    def _validate_inputs(self, m_star, r_d, sigma_1au, FFUV_G0, n_points, 
                    m_dot, gamma, p, q, gridsize):
        """Validate input parameters and raise errors for invalid values."""
        
        # Physical parameter checks
        if m_star < 0.3 or m_star > 3.0:
            raise ValueError(f"Stellar mass must be between 0.3 to 3.0 M_sun")
        
        if r_d < 20 or r_d > 150:
            raise ValueError(f"Disk radius must be between 10 to 150 au")
        
        if sigma_1au < 10 or sigma_1au > 1e4:
            raise ValueError(f"Surface density must be between 10 to 10,000 g/cm2")
        
        if FFUV_G0 < 1e2 or FFUV_G0 > 1e5:
            raise ValueError(f"FUV field strength must be between 100 to 100,000 G0")
        
        # Grid checks
        if n_points < 200:
            raise ValueError(f"Need at least 200 grid points for reliable computation. Regrid afterwards if desired.")
        
        if gridsize is not None and gridsize <= r_d:
            raise ValueError(f"Grid size ({gridsize} AU) must be larger than disk radius ({r_d} AU)")
        
        # Optional parameter checks
        if m_dot is not None:
            if m_dot <= 0:
                raise ValueError(f"Mass loss rate must be positive (got m_dot={m_dot})")
            if m_dot > 1e-5:
                print(f"WARNING: Unusually high mass loss rate (M_dot = {m_dot:.2e} M_sun/yr). "
                    f"Typical values are < 1e-5 M_sun/yr")
        
        if gamma is not None and gamma <= 0:
            raise ValueError(f"Taper parameter must be positive (got gamma={gamma})")
            
        
    def _setup_grid(self):
        """Create logarithmic radial grid."""
        self.radius = np.logspace(np.log10(0.1), np.log10(self.gridsize), self.n_points)
        
    def _interpolate_mass_loss_rate(self):
        """Interpolate mass loss rate from input parameters."""
        if self.m_dot is None:
            m_dot_log = helpers.interpolate_mdot(self.m_star, self.r_d, 
                                                  self.sigma_1au, self.FFUV_G0)
            self.m_dot = 10.0**(m_dot_log)
    
    def _compute_gamma(self):
        """Compute exponential taper parameter from Keyte & Haworth (2026) prescription."""
        if self.gamma is None:
            self.gamma = (0.77 * self.m_star**0.10 * self.r_d**0.78 * 
                         (self.FFUV_G0 / 100)**0.07)
    
    def _compute_physical_properties(self):
        """Compute PDR temperature, sound speed, and geometric factor."""
        if self.verbose:
            print(' > Computing physical properties...')

        t0_PDR      = 200.0
        self.t_PDR  = min(max(t0_PDR * (self.FFUV_G0 / 1000.0)**0.2, 10), 3000.0)
        self.cs_PDR = np.sqrt((self.K_B * self.t_PDR) / (self.MU_PDR * self.M_P))
        
        # Calculate geometric factor at disk edge
        t_1au = 150.0 * self.m_star**0.25
        t_d   = max(t_1au * self.r_d**(-0.5), 10)
        cs_d  = np.sqrt((self.K_B * t_d) / (self.MU_DISK * self.M_P))
        omega = np.sqrt((self.BIG_G * self.m_star * self.M_SUN) / (self.r_d * self.AU_CM)**3)
        h_d   = cs_d / omega
        self.mathcal_F = h_d / np.sqrt((self.r_d * self.AU_CM)**2 + h_d**2)
        
        
    def _compute_disk_component(self):
        """Compute standard disk density and temperature profile."""
        if self.verbose:
            print(' > Computing disk density...')
        t_1au = 150.0 * self.m_star**0.25
        
        self.rho_disk    = np.zeros(self.n_points)
        self.temperature = np.zeros(self.n_points)
        
        for i in range(self.n_points):
            # Surface density with exponential taper
            sigma = (self.sigma_1au * self.radius[i]**(-1.0) * 
                    np.exp(-((self.radius[i]) / (self.r_d * 1.1))**self.gamma))
            
            # Temperature
            self.temperature[i] = max(t_1au * self.radius[i]**(-0.5), 10)
            
            # Scale height and midplane density
            cs    = np.sqrt((self.K_B * self.temperature[i]) / (self.MU_DISK * self.M_P))
            omega = np.sqrt((self.BIG_G * self.m_star * self.M_SUN) / (self.radius[i] * self.AU_CM)**3)
            H     = cs / omega
            self.rho_disk[i] = sigma / H
    
    def _compute_wind_component(self):
        """Compute spherically diverging wind density profile."""
        if self.verbose:
            print(' > Computing spherically diverging wind...')
        self.rho_wind = np.zeros(self.n_points)
        
        for i in range(self.n_points):
            if self.radius[i] >= self.r_d:
                self.rho_wind[i] = (self.m_dot * self.M_SUN / 
                                   (4.0 * np.pi * self.radius[i]**2 * self.AU_CM**2 * 
                                    self.YR_SEC * self.mathcal_F * self.cs_PDR))
    
    def _compute_plateau_component(self):
        """Compute plateau interpolation between disk and wind."""
        if self.verbose:
            print(' > Computing transition from disk to wind...')

        self.rho_plateau = np.zeros(self.n_points)
        idx_rd = (np.abs(self.radius - self.r_d)).argmin()
        
        for i in range(self.n_points):
            if self.radius[i] > self.r_d:
                x_norm = np.log10(self.radius[i] / self.r_d)
                lamda  = self.r_d**self.p + (self.FFUV_G0 / 100)**self.q
                f      = (1 - np.exp(-lamda * x_norm)) / (1 - np.exp(-lamda))
                self.rho_plateau[i] = (self.rho_disk[idx_rd] * (self.rho_wind[i] / self.rho_disk[idx_rd])**f)
    
    def _assemble_density(self):
        """Combine density components to create final density profile."""
        if self.verbose:
            print(' > Assembling total density structure...')

        self.density = np.zeros(self.n_points)
        
        for i in range(self.n_points):
            if self.radius[i] >= self.r_d:
                self.density[i] = np.maximum.reduce([self.rho_disk[i], 
                                                     self.rho_wind[i], 
                                                     self.rho_plateau[i]])
            else:
                self.density[i] = self.rho_disk[i]
    
    
    def compute(self):
        """
        Run the disk model calculation.
        
        Returns
        -------
        radius : ndarray
            Radial grid in AU
        density : ndarray
            Total midplane density in g/cm^3
            
        Notes
        -----
        Results are also stored as instance attributes:
        - self.radius, self.density
        - self.temperature
        - self.rho_disk, self.rho_wind, self.rho_plateau
        
        Examples
        --------
        import puffin_disk
        disk = puffin_disk.DiskModel1D(m_star=1.0, r_d=100, sigma_1au=100, FFUV_G0=1000)
        r, rho = disk.compute()
        plt.loglog(r, rho)
        """
        
        self._setup_grid()
        self._interpolate_mass_loss_rate()
        self._compute_gamma()
        
        if self.verbose:
            helpers.print_initialisation_1D(self.m_star, self.r_d, self.sigma_1au, self.FFUV_G0, self.n_points, self.gridsize, self.m_dot, self.gamma, self.p, self.q)
            helpers.log_section('COMPUTING 1D MODEL')
        
        self._compute_physical_properties()
        self._compute_disk_component()
        self._compute_wind_component()
        self._compute_plateau_component()
        self._assemble_density()
        
        if self.verbose:
            print(' > MODEL COMPLETE')

        return self.radius, self.density


    
    def __repr__(self):
        return (f"DiskModel1D(m_star={self.m_star}, r_d={self.r_d}, "
                f"sigma_1au={self.sigma_1au:.2e}, FFUV_G0={self.FFUV_G0})")


class DiskModel2D:
    """
    Calculate 2D (r, z) density structure for an externally irradiated protoplanetary disk.
    
    This class computes the gas density profile for a disk undergoing external 
    photoevaporation. The model solves for vertical hydrostatic equilibrium iteratively 
    and includes multiple components: a hydrostatic disk, a spherically diverging 
    photoevaporative wind launched from the disk surface, and a transition 
    region that smoothly connects the disk to the wind.
    
    Parameters
    ----------
    m_star : float
        Stellar mass in solar masses.
    r_d : float
        Characteristic disk radius (gravitational radius) in AU.
    sigma_1au : float
        Surface density at 1 AU in g/cm^2.
    FFUV_G0 : float
        External FUV field strength in Habing units (G0).
    n_points : int, optional
        Number of grid points in both radial and vertical directions (default: 1000).
        High resolution (>=1000) is recommended for accurate HSE solutions.
    gridsize : float, optional
        Outer radius of the computational grid in AU (default: 8 * r_d).
    m_dot : float, optional
        Mass loss rate in M_sun/yr. If None, interpolated from FRIED grid
        with a factor of 2 applied for 2D geometry.
    gamma : float, optional
        Exponential cutoff parameter for surface density profile. If None,
        calculated from scaling relation.
    p : float, optional
        Power law index for bowl transition (r_d component). Default is 0.2.
    q : float, optional
        Power law index for bowl transition (FFUV component). Default is 0.4.
    N_ITER : int, optional
        Number of hydrostatic equilibrium iterations (default: 20).
    
    Attributes
    ----------
    r_array : ndarray
        Radial grid in AU (logarithmically spaced).
    z_array : ndarray
        Vertical grid in AU (logarithmically spaced, includes midplane z=0).
    density : ndarray
        Total gas density in g/cm^3 at each (z, r) point. Shape is (nz, nr).
    temperature : ndarray
        Gas temperature in K at each (z, r) point. Shape is (nz, nr).
    rho_disk : ndarray
        Disk component density in g/cm^3. Shape is (nz, nr).
    rho_wind : ndarray
        Combined wind component density in g/cm^3. Shape is (nz, nr).
    model_name : str
        Model identifier string.
        
    Notes
    -----
    The density structure consists of:
    
    **Disk component:**
    - Power-law surface density with exponential truncation at r_d
    - Vertical structure solved iteratively assuming hydrostatic equilibrium
    - Temperature profile smoothly transitions from cool midplane to PDR temperature
    - Optical depth determines the disk surface (tau=1 surface)
    
    **Wind components:**
    - Spherical wind: Diverging outflow launched from tau=1 surface with density 
      proportional to 1/r^2, set by mass loss rate and sound speed
    - Transition region: Smooth transition from disk to spherical wind using exponential 
      blending function, tapered inside r_d
    - Wind is smoothly tapered inside r_d using cosine function
    
    The final density at each position is the maximum of the disk and wind components.
    
    **Hydrostatic Equilibrium:**
    The code iteratively solves the vertical structure by:
    1. Computing optical depths (vertical, radial inward/outward)
    2. Calculating FUV attenuation
    3. Determining temperature structure (disk + PDR)
    4. Solving d ln rho / dz = -Omega^2 z / c_s^2 - (1/T) dT/dz
    5. Renormalizing to match Sigma(r) at each radius
    
    The model aborts if:
    - Disk mass < 1e-5 M_sun
    - Disk is optically thin at r_d (tau < 1)
    """
    
    # Physical constants
    MU_DISK = 2.3
    MU_PDR  = 1.3
    BIG_G   = 6.67259e-8       # cm^3 g^-1 s^-2
    K_B     = 1.380626e-16     # erg K^-1
    M_SUN   = 1.9891e33        # g
    AU_CM   = 1.495979e13      # cm
    M_P     = 1.6726219e-24    # g
    YR_SEC  = 3.15576e7        # s
    
    # FUV opacity parameters
    SIGMA_FUV_WIND = 2.7e-23   # cm^2 per particle (wind)
    SIGMA_FUV_DISK = 8.0e-22   # cm^2 per particle (disk)
    
    
    def __init__(self, m_star, r_d, sigma_1au, FFUV_G0, n_points=1000, gridsize=None,
                 m_dot=None, gamma=None, p=0.2, q=0.4, k=1.75, N_ITER=20, verbose=True):
        
        # Initialise model
        self.verbose = verbose
        self._validate_inputs(m_star, r_d, sigma_1au, FFUV_G0, n_points, m_dot, gamma, p, q, k, gridsize)
        
        # Store input parameters
        self.m_star    = m_star
        self.r_d       = r_d
        self.sigma_1au = sigma_1au
        self.FFUV_G0   = FFUV_G0
        self.n_points  = n_points
        self.gridsize  = gridsize if gridsize is not None else r_d * 8
        self.m_dot     = m_dot
        self.gamma     = gamma
        self.p         = p
        self.q         = q
        self.k_smooth  = k
        self.N_ITER    = N_ITER
        
        # Derived quantities
        self.R_focal        = 0.5 * self.r_d     # focal point of the wind
        self.tau_surface    = 1.0                # tau value defining the disk surface
        self.kappa_FUV_disk = self.SIGMA_FUV_DISK / (self.MU_DISK * self.M_P)
        self.kappa_FUV_wind = self.SIGMA_FUV_WIND / (self.MU_PDR  * self.M_P)
        self.m_disk         = helpers.calculate_disk_mass(self.sigma_1au, r_max=self.r_d)
        self.model_name     = f'{m_star}_{FFUV_G0}_{r_d}_{sigma_1au}'
          
        # Initialize result arrays (computed when compute() is called)
        self.r_array     = None
        self.z_array     = None
        self.density     = None
        self.temperature = None
        self.rho_disk    = None
        self.rho_wind    = None
        
        
    def _validate_inputs(self, m_star, r_d, sigma_1au, FFUV_G0, n_points,
                         m_dot, gamma, p, q, k, gridsize):
        """Validate input parameters and raise errors for invalid values."""
        
        if m_star < 0.3 or m_star > 3.0:
            raise ValueError("Stellar mass must be between 0.3 to 3.0 M_sun")
        
        if r_d < 20 or r_d > 150:
            raise ValueError("Disk radius must be between 20 to 150 AU")
        
        if sigma_1au < 10 or sigma_1au > 1e4:
            raise ValueError("Surface density must be between 10 to 10,000 g/cm^2")
        
        if FFUV_G0 < 1e2 or FFUV_G0 > 1e5:
            raise ValueError("FUV field strength must be between 100 to 100,000 G0")
        
        if n_points < 1000:
            raise ValueError(f"Need at least 1000 grid points for reliable computation. Regrid afterwards if desired.")
        
        if gridsize is not None and gridsize <= r_d:
            raise ValueError(f"Grid size ({gridsize} AU) must be larger than disk radius ({r_d} AU)")
        
        if m_dot is not None:
            if m_dot <= 0:
                raise ValueError(f"Mass loss rate must be positive (got m_dot={m_dot})")
            if m_dot > 1e-5:
                print(f"WARNING: Unusually high mass loss rate (M_dot = {m_dot:.2e} M_sun/yr). "
                      f"Typical values are < 1e-5 M_sun/yr")
        
        if gamma is not None and gamma <= 0:
            raise ValueError(f"Taper parameter must be positive (got gamma={gamma})")
    
    
    def _setup_grid(self):
        """Create logarithmic radial and vertical grids."""
        self.r_array = np.logspace(np.log10(0.1), np.log10(self.gridsize), self.n_points)
        z_positive   = np.logspace(np.log10(0.01), np.log10(self.gridsize), self.n_points - 1)
        self.z_array = np.concatenate(([0.0], z_positive))  # include midplane as z=0
        
        self.nr = len(self.r_array)
        self.nz = len(self.z_array)
    
    
    def _interpolate_mass_loss_rate(self):
        """Interpolate mass loss rate from input parameters (with 2D geometry factor)."""
        if self.m_dot is None:
            factor_2d  = 2
            m_dot_log  = helpers.interpolate_mdot(self.m_star, self.r_d, self.sigma_1au, self.FFUV_G0)
            self.m_dot = factor_2d * 10**m_dot_log
    
    
    def _compute_gamma(self):
        """Compute exponential taper parameter from Keyte & Haworth (2026) prescription."""
        if self.gamma is None:
            self.gamma = (0.77 * self.m_star**0.10 * self.r_d**0.78 * (self.FFUV_G0 / 100)**0.07)
    
    
    def _compute_stellar_PDR_properties(self):
        """Compute stellar and PDR properties needed by the model."""
        if self.verbose:
            print(' > Computing stellar and PDR properties...')
        
        # Stellar properties
        L_star, R_star, T_eff = helpers.get_stellar_properties(self.m_star)
        fuv_frac              = helpers.fuv_fraction(T_eff)
        self.FUV_star         = fuv_frac * L_star * 3.828e33
        
        # Disk temperature normalization (80K from inspection of DALI models)
        self.t_1au = 80.0 * self.m_star**0.25
        
        # PDR temperature and sound speed
        T0_PDR      = 200.0
        self.t_PDR  = float(np.clip(T0_PDR * (self.FFUV_G0 / 1000.0)**0.2, 10.0, 3000.0))
        self.cs_PDR = np.sqrt((self.K_B * self.t_PDR) / (self.MU_PDR * self.M_P))
        
    
    def _compute_initial_disk(self): 
        """Set up initial disk density (Gaussian vertical profile) and FUV field."""
        if self.verbose:
            print(' > Computing initial disk structure...')
        
        nr, nz = self.nr, self.nz
        
        # Surface density profile
        self.sigma_r = (self.sigma_1au * self.r_array**(-1.0) * np.exp(-(self.r_array / (self.r_d))**self.gamma))
        
        # Allocate arrays
        self.rho_disk    = np.zeros((nz, nr))
        self.temperature = np.zeros((nz, nr))
        self.H           = np.zeros((nz, nr))
        FUV_field_g0     = np.zeros((nz, nr))
        
        # Initial Gaussian vertical profile
        for i in range(nr):
            Omega = np.sqrt((self.BIG_G * self.m_star * self.M_SUN) / (self.r_array[i] * self.AU_CM)**3)
            T_mid = max(self.t_1au * self.r_array[i]**(-0.5), 10.0)
            cs    = np.sqrt((self.K_B * T_mid) / (self.MU_DISK * self.M_P))
            H_i   = cs / Omega
            rho0  = self.sigma_r[i] / (np.sqrt(2.0 * np.pi) * H_i)
            
            for j in range(nz):
                zcm = self.z_array[j] * self.AU_CM
                self.H[j, i]        = H_i
                self.rho_disk[j, i] = rho0 * np.exp(-0.5 * (zcm / H_i)**2)
        
        # Unattenuated FUV field (stellar + external)
        for i in range(nr):
            for j in range(nz):
                zcm            = self.z_array[j] * self.AU_CM
                distance_cm    = np.sqrt((self.r_array[i] * self.AU_CM)**2 + zcm**2)
                FUV_star_local = self.FUV_star / (4.0 * np.pi * distance_cm**2)
                FUV_field_g0[j, i] = max(FUV_star_local / 1.6e-3 + self.FFUV_G0, 1e-30)
        
        self.FUV_field_g0 = FUV_field_g0
        
        # Index for r_d on the grid
        self.idx_rd = np.searchsorted(self.r_array, self.r_d, side="right") - 1
    
    
    def _compute_optical_depths(self, rho_field, kappa):
        """
        Compute optical depths from three directions and return the minimum.
        
        Parameters
        ----------
        rho_field : ndarray
            Density field, shape (nz, nr).
        kappa : float
            FUV opacity in cm^2/g.
        
        Returns
        -------
        tau_vert, tau_in, tau_out, tau_min : ndarray
            Optical depths from vertical, radial-inward, radial-outward directions,
            and their element-wise minimum. All shape (nz, nr).
        """
        nr, nz = self.nr, self.nz
        au = self.AU_CM
        
        tau_vert = np.zeros_like(rho_field)
        tau_in   = np.zeros_like(rho_field)
        tau_out  = np.zeros_like(rho_field)
        
        # Vertical (integrate from top down)
        for i in range(nr):
            cumulative = 0.0
            for j in range(nz - 2, -1, -1):
                dz_cm   = (self.z_array[j+1] - self.z_array[j]) * au
                rho_avg = 0.5 * (rho_field[j+1, i] + rho_field[j, i])
                cumulative += rho_avg * dz_cm * kappa
                tau_vert[j, i] = cumulative
            tau_vert[-1, i] = 0.0
        
        # Radial inward (integrate from outer edge inward)
        for j in range(nz):
            cumulative = 0.0
            for i in range(nr - 2, -1, -1):
                dr_cm   = (self.r_array[i+1] - self.r_array[i]) * au
                rho_avg = 0.5 * (rho_field[j, i+1] + rho_field[j, i])
                cumulative += rho_avg * dr_cm * kappa
                tau_in[j, i] = cumulative
            tau_in[j, -1] = 0.0
        
        # Radial outward (integrate from inner edge outward)
        for j in range(nz):
            cumulative = 0.0
            for i in range(1, nr):
                dr_cm   = (self.r_array[i] - self.r_array[i-1]) * au
                rho_avg = 0.5 * (rho_field[j, i] + rho_field[j, i-1])
                cumulative += rho_avg * dr_cm * kappa
                tau_out[j, i] = cumulative
            tau_out[j, 0] = 0.0
        
        tau_min = np.minimum.reduce([tau_vert, tau_in, tau_out])
        return tau_vert, tau_in, tau_out, tau_min
    
    
    def _iterate_hydrostatic_equilibrium(self):
        """
        Iteratively solve for vertical hydrostatic equilibrium.
        
        Each iteration:
        1. Computes optical depths from disk density
        2. Attenuates FUV field
        3. Sets PDR temperature
        4. Computes 2D temperature (midplane to PDR transition at tau=1 surface)
        5. Solves d ln(rho)/dz and renormalizes to match Sigma(r)
        
        Returns
        -------
        success : bool
            True if HSE converged without aborting.
        """
        nr, nz = self.nr, self.nz
        au     = self.AU_CM
        mid_j  = 0                   # midplane index
        
        T0_PDR   = 200.0
        TPDR_FUV = np.zeros((nz, nr))
        dT_dz    = np.zeros((nz, nr))
        
        for it in range(self.N_ITER):
            if self.verbose:
                print(f'\r > Computing hydrostatic equilibrium ({it+1}/{self.N_ITER})...', end='', flush=True)

            
            # 1. Optical depths
            tau_vert, tau_in, tau_out, tau_min = self._compute_optical_depths(
                self.rho_disk, self.kappa_FUV_disk)
            
            # -- Find tau=1 surface at midplane (scanning from outside in)
            tau_midplane = tau_min[mid_j, :]
            idx_outside_in = np.where(tau_midplane[::-1] >= 1)[0]
            
            if idx_outside_in.size > 0:
                self.idx_tau1_midplane = len(tau_midplane) - 1 - idx_outside_in[0]
            else:
                self.idx_tau1_midplane = None
            
            # -- Abort checks
            if self.idx_tau1_midplane is None:
                self._abort_reason = "ABORTED: Disk optically thin (tau_rd=0)"
                return False
            elif self.r_array[self.idx_tau1_midplane] < self.r_array[self.idx_rd]:
                self._abort_reason = "ABORTED: Disk optically thin (tau=1 inside r_d)"
                return False
            
            # 2. Attenuated FUV field 
            FUV_att_g0 = self.FUV_field_g0 * np.exp(-tau_min)
            
            # 3. PDR temperature field (external FUV only; stellar contribution not included)
            # -- Note: Replace self.FFUV_G0 with FUV_att_g0[j,i] to include stellar heating  (not supported)
            for i in range(nr):
                for j in range(nz):
                    TPDR_FUV[j, i] = float(np.clip(
                        T0_PDR * (self.FFUV_G0 / 1000.0)**0.2, 10.0, 3000.0))     
            
            #  4. 2D temperature structure
            for i in range(nr):
                # -- Find tau=1 surface height at this radius
                z_surface = None
                for j in range(nz - 1, -1, -1):
                    if tau_vert[j, i] >= self.tau_surface:
                        z_surface = self.z_array[j]
                        break
                
                T_midplane = max(self.t_1au * self.r_array[i]**(-0.5), 10.0)
                
                if z_surface is None or z_surface == 0:
                    for j in range(nz):
                        self.temperature[j, i] = TPDR_FUV[j, i]
                elif z_surface > 0.0:
                    for j in range(nz):
                        if self.z_array[j] <= z_surface:
                            z_norm = self.z_array[j] / z_surface
                            self.temperature[j, i] = helpers.smooth_temperature_profile(
                                z_norm, T_midplane, TPDR_FUV[j, i], "tanh", k=self.k_smooth)
                        else:
                            self.temperature[j, i] = TPDR_FUV[j, i]
                else:
                    for j in range(nz):
                        print(f'No temperature set at r={self.r_array[i]}au z={self.z_array[j]}au')
            
            # 5. Solve HSE
            
            # -- Temperature gradients dT/dz (K/AU)
            for i in range(nr):
                for j in range(nz):
                    if j == 0:
                        dT_dz[j, i] = ((self.temperature[j+1, i] - self.temperature[j, i]) /
                                        (self.z_array[j+1] - self.z_array[j]))
                    elif j == nz - 1:
                        dT_dz[j, i] = ((self.temperature[j, i] - self.temperature[j-1, i]) /
                                        (self.z_array[j] - self.z_array[j-1]))
                    else:
                        dT_dz[j, i] = ((self.temperature[j+1, i] - self.temperature[j-1, i]) /
                                        (self.z_array[j+1] - self.z_array[j-1]))
            
            # -- Integrate d ln(rho)/dz upward, then renormalize to Sigma(r)
            for i in range(nr):
                Omega = np.sqrt((self.BIG_G * self.m_star * self.M_SUN) /
                               (self.r_array[i] * au)**3)
                
                rho_mid = max(self.rho_disk[mid_j, i], 1e-40)
                ln_rho  = np.log(rho_mid)
                self.rho_disk[mid_j, i] = rho_mid
                
                for j in range(1, nz):
                    dz_cm   = (self.z_array[j] - self.z_array[j-1]) * au
                    z_cm    = self.z_array[j-1] * au
                    T_loc   = max(self.temperature[j-1, i], 10.0)
                    cs2     = (self.K_B * T_loc) / (self.MU_DISK * self.M_P)
                    dTdz_cm = dT_dz[j-1, i] / au  # K/AU -> K/cm
                    
                    rhs = -(Omega**2 * z_cm / cs2) - (dTdz_cm / T_loc)
                    ln_rho += rhs * dz_cm
                    self.rho_disk[j, i] = np.exp(ln_rho)
                
                # -- Renormalize to match surface density
                col = np.trapz(self.rho_disk[:, i], self.z_array * au)
                if col > 0.0:
                    self.rho_disk[:, i] *= (self.sigma_r[i] / col)
        
        if self.verbose:
            print()       # newline after progress counter
        
        return True
    
    
    def _compute_spherical_wind(self):
        """Compute spherically diverging wind launched from the tau=1 surface."""
        if self.verbose:
            print(' > Computing spherical wind...')
        
        nr, nz = self.nr, self.nz
        au     = self.AU_CM
        
        rho_wind_sph = np.zeros((nz, nr))
        rho_tau1_midplane = self.rho_disk[0, self.idx_tau1_midplane]
        
        # Use exponential wind where available, else disk density for tau surface
        rho_for_tau = self.rho_disk.copy()
        
        # Identify region interior to tau=1 surface
        base_mask = (self.rho_disk >= rho_tau1_midplane)
        surface_indices = np.argwhere(base_mask)
        
        if surface_indices.size == 0:
            self.rho_wind_sph = rho_wind_sph
            self.surface_rhos = np.array([])
            return
        
        surface_coords = np.column_stack((
            self.r_array[surface_indices[:, 1]],
            self.z_array[surface_indices[:, 0]]))
        surface_rhos = rho_for_tau[surface_indices[:, 0], surface_indices[:, 1]]
        
        # KD-tree for nearest distance to base region
        tree = cKDTree(surface_coords)
        
        for i in range(nr):
            for j in range(nz):
                dist, idx = tree.query([self.r_array[i], self.z_array[j]])
                
                if dist > 0:
                    total_dist = dist + self.R_focal
                    rho_wind_sph[j, i] = ((self.m_dot * self.M_SUN) /
                        (4.0 * np.pi * total_dist**2 * au**2 * self.YR_SEC * self.cs_PDR))
                else:
                    rho_wind_sph[j, i] = 1e-30
        
        # Rescale wind if density at wind base exceeds disk density in the same region
        surface_rho_min = np.nanmin(surface_rhos[surface_rhos > 1e-30])
        max_wind = np.nanmax(rho_wind_sph)
        wind_scaling = max_wind / surface_rho_min
        
        if wind_scaling > 1:
            rho_wind_sph /= wind_scaling
            if self.verbose:
                print(f'   - Wind scaled by x{wind_scaling:.2f}')
        
        # Store for transition region computation
        self.rho_wind_sph = rho_wind_sph
        self.surface_rhos = surface_rhos
    
    
    def _compute_disk_wind_transition(self):
        """Compute transition region between disk and wind."""
        if self.verbose:
            print(' > Computing disk-wind transition...')
        
        nr, nz = self.nr, self.nz
        rho_tau1_midplane = self.rho_disk[0, self.idx_tau1_midplane]
        
        rho_bowl = np.zeros((nz, nr))
        
        # Identify base region (same mask as spherical wind)
        base_mask = (self.rho_disk >= rho_tau1_midplane)
        surface_indices = np.argwhere(base_mask)
        
        if surface_indices.size == 0:
            self.rho_bowl = rho_bowl
            return
        
        surface_coords = np.column_stack((
            self.r_array[surface_indices[:, 1]],
            self.z_array[surface_indices[:, 0]]))
        
        tree = cKDTree(surface_coords)
        
        for i in range(nr):
            for j in range(nz):
                dist, idx = tree.query([self.r_array[i], self.z_array[j]])
                rho_base = self.surface_rhos[idx]
                
                if dist > 0:
                    # Transition function (2D analogue of 1D plateau)
                    x_norm = np.log10(1 + dist / self.r_d)
                    lamda = self.r_d**self.p + (self.FFUV_G0 / 100)**self.q
                    f = (1 - np.exp(-lamda * x_norm)) / (1 - np.exp(-lamda))
                    
                    rho_bowl[j, i] = rho_base * (self.rho_wind_sph[j, i] / rho_base)**f
                else:
                    rho_bowl[j, i] = 1e-30
                
                # Taper inside r_d
                if self.r_array[i] < self.r_d:
                    rho_bowl[j, i] *= (self.r_array[i] / self.r_d)
        
        self.rho_bowl = rho_bowl
    
    
    def _taper_and_assemble_wind(self):
        """Combine wind components, apply spatial taper, and assemble final density."""
        if self.verbose:
            print(' > Tapering wind and assembling total density...')
        
        nr, nz = self.nr, self.nz
        rho_tau1_midplane = self.rho_disk[0, self.idx_tau1_midplane]
        
        # Combined wind = max(spherical, transition)
        self.rho_wind = np.maximum(self.rho_wind_sph, self.rho_bowl)
        
        # Apply spatial taper (cosine tapering near inner edge)
        for i in range(nr):
            for j in range(nz):
                if self.z_array[j] >= self.r_d / 2:
                    edge = self.r_d
                else:
                    edge = self.r_d * (self.z_array[j] / (self.r_d / 2))
                
                f_scale = helpers.smooth_cosine(self.r_array[i], edge, steepness=2)
                self.rho_wind[j, i] *= f_scale
        
        # Floor disk density below tau=1 threshold
        self.rho_disk[self.rho_disk < rho_tau1_midplane] = 1e-30
        
        # Final density = max(disk, wind)
        self.density = np.maximum(self.rho_disk, self.rho_wind)
    
    
    def compute(self):
        """
        Run the 2D disk model calculation.
        
        Returns
        -------
        r_array : ndarray
            Radial grid in AU.
        z_array : ndarray
            Vertical grid in AU.
        density : ndarray
            Total gas density in g/cm^3, shape (nz, nr).
            
        Returns (on abort)
        ------------------
        result : str
            Error message if model fails.
        model_name : str
            Model identifier string.
        
        Notes
        -----
        Results are also stored as instance attributes:
        - self.r_array, self.z_array, self.density
        - self.temperature
        - self.rho_disk, self.rho_wind
        
        Examples
        --------
        import puffin_disk
        disk = puffin_disk.DiskModel2D(m_star=1.0, r_d=100, sigma_1au=100, FFUV_G0=1000)
        r, z, rho = disk.compute()
        plt.contourf(r, z, np.log10(rho), levels=np.arange(-20,-11,0.2), cmap='Spectral_r', extend='both')
        """
        
        # Setup
        self._setup_grid()
        self._interpolate_mass_loss_rate()
        self._compute_gamma()
        
        # Check disk mass
        if self.m_disk < 1e-5:
            return 'ABORTED: Disk mass < 1e-5 M_sun', self.model_name
        
        if self.verbose:
            helpers.print_initialisation_2D(self.m_star, self.r_d, self.sigma_1au, self.FFUV_G0,
                                          self.n_points, self.gridsize, self.m_dot, self.gamma,
                                          self.p, self.q, self.k_smooth)
            helpers.log_section('COMPUTING 2D MODEL')
        
        # Star, PDR, and disk
        self._compute_stellar_PDR_properties()
        self._compute_initial_disk()
        
        # Iterate hydrostatic equilibrium
        success = self._iterate_hydrostatic_equilibrium()
        if not success:
            return self._abort_reason, self.model_name
        
        # Wind components
        self._compute_spherical_wind()
        self._compute_disk_wind_transition()
        self._taper_and_assemble_wind()
        
        if self.verbose:
            print(' > MODEL COMPLETE')
        
        return self.r_array, self.z_array, self.density
    
    
    def __repr__(self):
        return (f"DiskModel2D(m_star={self.m_star}, r_d={self.r_d}, "
                f"sigma_1au={self.sigma_1au:.2e}, FFUV_G0={self.FFUV_G0})")