"""
PUFFIN Helper Functions
========================
Supporting utilities for the PUFFIN (Python Utility For FUV Irradiated disk 
deNsities) parametric disk modeling framework.

This module provides essential computational and I/O utilities for modeling
externally FUV-irradiated protoplanetary disks, including stellar property
calculations, temperature profiles, mass-loss rate interpolation, and formatted
console output.

Module Contents
---------------
Temperature Profiles:
    smooth_temperature_profile : Create smooth vertical temperature transitions

Stellar Properties (Eker et al. 2018 empirical relations):
    eker_mlr : Six-piece mass-luminosity relation
    get_stellar_properties : Compute luminosity, radius, and effective temperature
    planck_wavelength : Planck function in wavelength form
    fuv_fraction : Calculate FUV (912-2000 Å) luminosity fraction

Disk Properties:
    calculate_disk_mass : Integrate total disk mass from surface density profile

Mass-Loss Rate Interpolation:
    interpolate_mdot : Query FRIED grid lookup table for photoevaporation rates

Mathematical Utilities:
    smooth_cosine : Adjustable cosine taper function for smooth transitions

Console Output:
    print_initialisation_1D/2D : Formatted model initialization displays
    log_section/table_* : Structured table output with Unicode formatting
    format_scientific : Format numbers with superscript exponents

Dependencies
------------
numpy : Array operations and mathematical functions
scipy.integrate : Numerical integration (Planck function, disk mass)
scipy.interpolate : N-dimensional interpolation for mass-loss rates
scipy.constants : Physical constants (h, c, k)
matplotlib.pyplot : Plotting utilities (imported but not actively used here)
pathlib : Cross-platform file path handling
platform : Operating system detection for Unicode box-drawing

Data Files
----------
FRIEDV2_ALL_fPAH1p0_growth.dat : Pre-computed mass-loss rate lookup table
    from the FRIED photoevaporation grid (Haworth et al. 2018, 2023)

References
----------
Keyte & Haworth (2026) : PUFFIN overview paper
Eker et al. (2018) : Empirical mass-luminosity, mass-radius, and mass-temperature relations
Haworth et al. (2018, 2023) : FRIED photoevaporation grid

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.constants import h, c, k
from scipy.interpolate import LinearNDInterpolator
from pathlib import Path
import platform

#######################
# TEMPERATURE PROFILE #
#######################

def smooth_temperature_profile(z_norm, T_midplane, T_surface, profile_type="tanh", k=2):
    """
    Create a smooth temperature profile from midplane to surface.
    
    Parameters:
    z_norm: normalized height (0 at midplane, 1 at tau surface)
    T_midplane: temperature at the midplane
    T_surface: temperature at the tau surface (TPDR)
    profile_type: type of smooth function to use
    
    Currently only the 'tanh' option has been validated.
    """
    
    if profile_type == "tanh":
        # Hyperbolic tangent provides very smooth transition
        # The factor k controls the steepness - larger values make it more linear
        smooth_factor = 0.5 * (1 + np.tanh(k * (z_norm - 0.5)))
        return T_midplane + (T_surface - T_midplane) * smooth_factor
    
    elif profile_type == "exponential":
        if z_norm == 0:
            return T_midplane
        else:
            alpha = 2.0  # Controls how quickly temperature rises with height
            smooth_factor = (1 - np.exp(-alpha * z_norm)) / (1 - np.exp(-alpha))
            return T_midplane + (T_surface - T_midplane) * smooth_factor
    
    elif profile_type == "cubic":
        smooth_factor = 3 * z_norm**2 - 2 * z_norm**3
        return T_midplane + (T_surface - T_midplane) * smooth_factor
    
    elif profile_type == "linear":
        return T_midplane + (T_surface - T_midplane) * z_norm
    
    else:  # default to linear if unknown type specified
        return T_midplane + (T_surface - T_midplane) * z_norm
    
    

##################################################
# STELLAR FUV CALCULATION FUNCTIONS              #
# Uses MLR, MRR, and MTR from Eker et al. (2018) #
##################################################

# Mass-luminosity
def eker_mlr(mass):
    """
    Eker et al. (2018) six-piece mass-luminosity relation.
    
    Parameters:
    -----------
    mass : float
        Stellar mass in solar masses (M/M☉)
        Valid range: 0.179 ≤ M ≤ 31 M☉
    
    Returns:
    --------
    log_luminosity : float
        log10(L/L☉) - logarithm of luminosity in solar units
    """
    if mass < 0.179:
        mass = 0.179
        print(f"Warning: Mass {mass} below valid range (<0.179 Msun). Using M=0.179 for luminosity calculation.")
    elif mass > 31.0:
        mass = 31.0
        print(f"Warning: Mass {mass} above valid range (>31.0 Msun). Using M=31.0 for luminosity calculation.")
    
    log_mass = np.log10(mass)
    
    # Simple if/elif statements for the six domains
    if mass <= 0.45:
        # Ultra low mass: 0.179 < M* < 0.45
        log_luminosity = 2.028 * log_mass - 0.976
    elif mass <= 0.72:
        # Very low mass: 0.45 < M* < 0.72
        log_luminosity = 4.572 * log_mass - 0.102
    elif mass <= 1.05:
        # Low mass: 0.72 < M* < 1.05
        log_luminosity = 5.743 * log_mass - 0.007
    elif mass <= 2.40:
        # Intermediate mass: 1.05 < M* < 2.40
        log_luminosity = 4.329 * log_mass + 0.010
    elif mass <= 7:
        # High mass: 2.4 < M* < 7.0
        log_luminosity = 3.967 * log_mass + 0.093
    else:
        # Very high mass: 7.0 < M* < 31.0
        log_luminosity = 2.865 * log_mass + 1.105
    
    return log_luminosity


# MLR, MRR, MTR
def get_stellar_properties(mass):
    """
    Get stellar properties empirical relations from Eker aet al. (2018)
    
    Parameters:
    -----------
    mass : float
        Stellar mass in solar masses (M/M☉)
    
    Returns:
    --------
    L : float
        Luminosity in solar units (L/L☉)
    R : float  
        Radius in solar units (R/R☉)
    T_eff : float
        Effective temperature in Kelvin
    """
    # Physical constants (SI units)
    sigma_sb = 5.67e-8  # Stefan-Boltzmann constant [W m^-2 K^-4]
    L_sun = 3.828e26    # Solar luminosity [W]
    R_sun = 6.96e8      # Solar radius [m]
    
    # Get luminosity from MLR (in solar units)
    log_L = eker_mlr(mass)
    if log_L is None:
        return None, None, None
    
    L_solar = 10**log_L  # Luminosity in L☉
    L_watts = L_solar * L_sun  # Convert to watts
    
    if mass <= 1.5:
        # Use empirical MRR to get radius in solar units
        R_solar = 0.438 * mass**2 + 0.479 * mass + 0.075
        R_meters = R_solar * R_sun  # Convert to meters
        
        # L = 4πR^2sigmaT^4, solve for T.
        T_eff = (L_watts / (4 * np.pi * R_meters**2 * sigma_sb))**(1/4)
        
    else:
        # Use empirical MTR to get temperature
        log_mass = np.log10(mass)
        log_T_eff = -0.170 * log_mass**2 + 0.888 * log_mass + 3.671
        T_eff = 10**log_T_eff
        
        # Stefan-Boltzmann law: L = 4πR²σT⁴
        # Solve for R: R = sqrt(L / (4πσT⁴))
        R_meters = (L_watts / (4 * np.pi * sigma_sb * T_eff**4))**0.5
        R_solar = R_meters / R_sun  
    
    return L_solar, R_solar, T_eff


# Planck function
def planck_wavelength(wavelength, T):
    """
    Planck function in wavelength form
    wavelength: in meters
    T: temperature in K
    Returns: spectral radiance per unit wavelength
    """
    return (2 * h * c**2 / wavelength**5) / (np.exp(h * c / (wavelength * k * T)) - 1)


# Calculate FUV fraction
def fuv_fraction(T_eff):
    """
    Calculate FUV fraction by integrating Planck function
    T_eff: effective temperature in K
    Returns: fraction of luminosity in FUV band (912-2000 Å)
    """
    # Convert Angstroms to meters
    lambda_min = 912e-10   # 912 Å
    lambda_max = 2000e-10  # 2000 Å
    
    # Integrate FUV portion
    fuv_integral, _ = quad(lambda lam: planck_wavelength(lam, T_eff), 
                          lambda_min, lambda_max)
    
    # Integrate total (approximate with reasonable bounds)
    # For practical purposes, integrate from 100 Å to 100 μm
    total_integral, _ = quad(lambda lam: planck_wavelength(lam, T_eff), 
                            100e-10, 100e-5)
    
    return fuv_integral / total_integral


#######################
# CALCULATE DISK MASS #
#######################

def calculate_disk_mass(sigma_1au, r_min=0.001, r_max=100.0, n_points=1000):
    """
    Calculate total disk mass for protoplanetary disk with power-law surface density.
    
    Parameters:
    sigma_1au: reference surface density at 1 AU (g/cm²)
    r_min: inner radius (AU) - default 0.1 to avoid divergence at r=0
    r_max: outer radius (AU) - default 100
    n_points: number of integration points
    
    Returns:
    mass in solar masses
    """
    # Create radial array
    r_array = np.linspace(r_min, r_max, n_points)
    
    # Surface density: sigma = sigma_1au * r^(-1)
    sigma = sigma_1au * (r_array**(-1.0))
    
    # Mass element: dm = sigma * 2π * r * dr
    dr = r_array[1] - r_array[0]  
    mass_elements = sigma * 2 * np.pi * r_array * dr
    
    # Total mass in grams
    total_mass_g = np.sum(mass_elements)
    
    # Convert AU² to cm² and then to solar masses
    au_to_cm = 1.496e13
    solar_mass_g = 1.989e33
    
    total_mass_g *= (au_to_cm**2) 
    total_mass_solar = total_mass_g / solar_mass_g
    
    return total_mass_solar


#################################
# INTERPOLATE MASS-LOSS RATE    #
# Uses data from the FRIED grid #
#################################

def interpolate_mdot(m_star, r_d, sigma_1au, F_FUV):
    """
    Interpolate log10(m_dot) from the lookup table.
    
    Parameters:
    -----------
    m_star : float or array
        Stellar mass (0.3 - 3.0)
    r_d : float or array
        Disk radius (10 - 150)
    sigma_1au : float or array
        Surface density at 1 AU (1e1 - 1e4)
    F_FUV : float or array
        FUV flux (1e2 - 1e5)
    
    Returns:
    --------
    log10(m_dot) : float or array
    """
    
    # Get data file
    helpers_dir = Path(__file__).parent
    data_file = helpers_dir / 'data' / 'FRIEDV2_ALL_fPAH1p0_growth.dat'
    data = np.loadtxt(data_file, skiprows=1, delimiter=',')
    
    # Extract columns (excluding sigma_rd)
    m_star_list = data[:, 0]
    r_d_list = data[:, 1]
    sigma_1au_list = data[:, 2]
    # sigma_rd_list = data[:, 3]  # IGNORED
    F_FUV_list = data[:, 4]
    log_mdot_list = data[:, 5]
    
    # Transform to log space for parameters spanning orders of magnitude
    points = np.column_stack([
        m_star_list,              # 0.3 - 3.0
        r_d_list,                 # 10 - 150
        np.log10(sigma_1au_list), # 1e1 - 1e4 (log space)
        np.log10(F_FUV_list)      # 1e2 - 1e5 (log space)
    ])
    
    query = np.array([[
        m_star,
        r_d,
        np.log10(sigma_1au),
        np.log10(F_FUV)
    ]])
    
    interpolator = LinearNDInterpolator(points, log_mdot_list)
    
    result = interpolator(query)[0]
    
    if np.isnan(result):
        print("Warning: Requested point is outside the interpolation domain")
    
    return result


################
# COSINE TAPER #
################

def smooth_cosine(x, x1, steepness):
    # Handle edge case: if x1=0, return 1.0 (no tapering)
    if x1 == 0:
        return 1.0
    # Clamp to [0, x1]
    x = np.clip(x, 0, x1)
    # Apply adjustable exponent to control curve shape
    t = (x / x1) ** steepness
    return (1 - np.cos(np.pi * t)) / 2


############################
# CONSOLE OUTPUT FUNCTIONS #
############################

# Safe message formatting for Windows/Mac
def safe_logger(message):
    """
    Safe logging with ASCII replacements on Windows.
    """
    if platform.system() == 'Windows':
        replacements = {
            '▓': '#', '►': '>', '│':'|',
            '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5',
            '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9', '⁰': '0',
            '⁻': '-', '×': 'x', '•': '-'
        }

        for unicode_char, ascii_char in replacements.items():
            message = message.replace(unicode_char, ascii_char)
        
        # Handle any remaining Unicode
        safe_message = message.encode('ascii', 'replace').decode('ascii')
        print(safe_message)
    else:
        print(message)

# 1D code header
def print_initialisation_1D(m_star, r_d, sigma_1au, FFUV_G0, n_points, gridsize, m_dot, gamma, p, q):
    if platform.system() == 'Windows':
        print("|" + ("-" * 68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + "PUFFIN".center(68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + "Python Utility For Fuv Irradiated disk deNsities".center(68) + "|")
        print("|" + "by Luke Keyte".center(68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + "Version 1.0.0 | 2026".center(68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + ("-" * 68) + "|")
        print("\n")
    else:
        print("┏" + ("━" * 68) + "┓")
        print("┃" + "".center(68) + "┃")
        print("┃" + "      ▗▀▖▗▀▖▗    ".center(68) + "┃")
        print("┃" + "▛▀▖▌ ▌▐  ▐  ▄ ▛▀▖".center(68) + "┃")
        print("┃" + "▙▄▘▌ ▌▜▀ ▜▀ ▐ ▌ ▌".center(68) + "┃")
        print("┃" + "▌  ▝▀▘▐  ▐  ▀▘▘ ▘".center(68) + "┃")
        print("┃" + "".center(68) + "┃")
        print("┃" + "Python Utility For Fuv Irradiated disk deNsities".center(68) + "┃")
        print("┃" + "by Luke Keyte".center(68) + "┃")
        print("┃" + "".center(68) + "┃")
        print("┃" + "Version 1.0.0 | 2026".center(68) + "┃")
        print("┃" + "".center(68) + "┃")
        print("┗" + ("━" * 68) + "┛")
        print("\n")
    # Print model properties
    log_section('INITIALIZATION')
    log_table_header(' > Model initialized')
    log_table_row("Stellar mass", m_star, 'M_sun')
    log_table_row("Disk radius", r_d, 'AU')
    log_table_row("Sigma (1AU)", sigma_1au, 'g cm^-2')
    log_table_row("FUV field", FFUV_G0, 'G0')
    log_table_row("N grid cells", n_points)
    log_table_row("Grid size", gridsize, 'AU')
    log_table_row("Mass-loss rate", m_dot, 'M_sun/yr')
    log_table_row("Gamma", gamma)
    log_table_row("p", p)
    log_table_row("q", q)
    log_table_footer()
    
# 2D code header
def print_initialisation_2D(m_star, r_d, sigma_1au, FFUV_G0, n_points, gridsize, m_dot, gamma, p, q, k):
    if platform.system() == 'Windows':
        print("|" + ("-" * 68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + "PUFFIN".center(68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + "Python Utility For Fuv Irradiated disk deNsities".center(68) + "|")
        print("|" + "by Luke Keyte".center(68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + "Version 1.0.0 | 2026".center(68) + "|")
        print("|" + "".center(68) + "|")
        print("|" + ("-" * 68) + "|")
        print("\n")
    else:
        print("┏" + ("━" * 68) + "┓")
        print("┃" + "".center(68) + "┃")
        print("┃" + "      ▗▀▖▗▀▖▗    ".center(68) + "┃")
        print("┃" + "▛▀▖▌ ▌▐  ▐  ▄ ▛▀▖".center(68) + "┃")
        print("┃" + "▙▄▘▌ ▌▜▀ ▜▀ ▐ ▌ ▌".center(68) + "┃")
        print("┃" + "▌  ▝▀▘▐  ▐  ▀▘▘ ▘".center(68) + "┃")
        print("┃" + "".center(68) + "┃")
        print("┃" + "Python Utility For Fuv Irradiated disk deNsities".center(68) + "┃")
        print("┃" + "by Luke Keyte".center(68) + "┃")
        print("┃" + "".center(68) + "┃")
        print("┃" + "Version 1.0.0 | 2026".center(68) + "┃")
        print("┃" + "".center(68) + "┃")
        print("┗" + ("━" * 68) + "┛")
        print("\n")
    # Print model properties
    log_section('INITIALIZATION')
    log_table_header(' > Model initialized')
    log_table_row("Stellar mass", m_star, 'M_sun')
    log_table_row("Disk radius", r_d, 'AU')
    log_table_row("Sigma (1AU)", sigma_1au, 'g cm^-2')
    log_table_row("FUV field", FFUV_G0, 'G0')
    log_table_row("N grid cells", n_points)
    log_table_row("Grid size", gridsize, 'AU')
    log_table_row("Mass-loss rate", m_dot, 'M_sun/yr')
    log_table_row("Gamma", gamma)
    log_table_row("p", p)
    log_table_row("q", q)
    log_table_row("k", k)
    log_table_footer()

def log_section(title):
    safe_logger("▓" * 24 + title.upper().center(22) + "▓" * 24 +"\n")
    
def log_table_header(title):
    """Log a table section header."""
    safe_logger(f" {title}")
    print(" ┌────────────────────┬───────────────┬──────────┐")
    print(" │ Parameter          │ Value         │ Unit     │")
    print(" ├────────────────────┼───────────────┼──────────┤")
    
def log_table_row(parameter, value, unit=""):
    """Format and log a parameter as a table row."""
    formatted_value = format_scientific(value) if isinstance(value, float) else value
    # Handle boolean values
    if isinstance(value, bool):
        formatted_value = str(value)
    
    if unit:
        # Convert standard units to proper Unicode
        unit = (unit.replace("^-1", "⁻¹")
                   .replace("^-2", "⁻²")
                   .replace("^-3", "⁻³")
                   .replace("^2", "²")
                   .replace("^3", "³"))
        safe_logger(f" │ {parameter:<18} │ {formatted_value:<13} │ {unit:<8} │")
    else:
        safe_logger(f" │ {parameter:<18} │ {formatted_value:<13} │          │")

def log_table_footer():
    """Log the table footer."""
    print(" └────────────────────┴───────────────┴──────────┘\n")
    
def format_scientific(number):
    """Format a number in scientific notation with proper symbols."""
    # For small numbers close to zero, just return the formatted number
    if abs(number) < 0.1 and abs(number) > 0.0001:
        return f"{number:.4f}"
    
    # For numbers that don't need scientific notation
    if abs(number) >= 0.1 and abs(number) < 1000:
        # Format with appropriate decimal places
        if abs(number) >= 100:
            return f"{number:.1f}"
        elif abs(number) >= 10:
            return f"{number:.2f}"
        else:
            return f"{number:.3f}"
    
    # Convert to scientific notation
    sci_notation = f"{number:.2e}"
    parts = sci_notation.split('e')
    mantissa = float(parts[0])
    exponent = int(parts[1])
    
    # Format exponent with superscripts
    exponent_map = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', 
                    '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', 
                    '-': '⁻', '+': ''}
    exponent_str = ''.join(exponent_map[char] for char in str(exponent))
    
    # Return formatted string
    if platform.system() == 'Windows':
        return f"{mantissa:.2f}e{exponent_str}"
    else:
        return f"{mantissa:.2f} × 10{exponent_str}"



        