# This module deals with input/output for the parametric disk + wind model
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.constants import h, c, k
from scipy.interpolate import LinearNDInterpolator



# Read input file and set up parameters
def GetInputParameters():
    paramFile = "inputparams.dat"
    parameters = np.genfromtxt(paramFile, delimiter=",")
    Mstar=parameters[0,1]
    Rd=parameters[1,1]
    Sigma1au=parameters[2,1]
    T1au=parameters[3,1]
    FUVG0=parameters[4,1]
    nR=parameters[5,1]
    Rout=parameters[6,1]
    print("=========================")
    print("Input parameters:")
    print("=========================")
    print("M_star      :", Mstar, "M_sun")
    print("R_d         :", Rd ,"au")
    print("R_out       :", Rout ,"au")
    print("Sigma_1au   :", Sigma1au , "g/cm^2")
    print("T_1au       :", T1au, "K")
    print("Ambient FUV :", FUVG0, "G0")
    print("n_cells (r) :", int(nR))
    print("=========================")
    print(" ")
    return Mstar, Rd, Sigma1au, T1au, FUVG0, nR, Rout


# Save results of 1D model to file 
def dump1Dmodel(filename, R, rho, T, vR, vPhi, sigmaFUV, dustToGas):
    np.savetxt(filename, np.c_[R, rho, T, vR, vPhi, sigmaFUV, dustToGas], delimiter=",", fmt='%f %e %f %f %f %f %f')


# Read FRIED 'radial' file
def readFriedRadial(filename):
    param_list = filename.split('_')
    m_star     = float(param_list[0].replace('Msol', '').replace('p', '.'))
    g0         = float(param_list[1].replace('G0', ''))
    sigma      = float(param_list[2].replace('Sigma', ''))
    rd         = float(param_list[3].replace('au', ''))
    return m_star, g0, sigma, rd


# Extract FRIED log_mdot
def readFriedMdot(filename, Mstar, Rd, Sigma1au, FUVG0): 
    with open(filename, 'r') as file:
        for line in file:
            current_line = line.split()

            if float(current_line[0].replace(',', '')) == Mstar:
                if float(current_line[1].replace(',', '')) == Rd:
                    if float(current_line[2].replace(',', '')) == Sigma1au:
                            if float(current_line[4].replace(',', '')) == FUVG0:
                                log_mdot = float(current_line[-1])
    return log_mdot


# Convert error values to 0
def convert_to_zero(x):
    try:
        return float(x)
    except ValueError:
        return 0
    
    
# Plot results
def plot_all(figname, xlim, ylim, r_array, z_array, rho_disk, rho_wind, rho, tau_FUV, FUV_field_attenuated, T, tau_surface, savefig=False):
    
    fig, ax = plt.subplots(2,3, figsize=(14,6.75))
    ax = ax.flatten()
    
    cmap = 'Spectral_r'
    
    levels_rho = np.arange(-25, -11, 0.2)
    
    rho_disk_plot = ax[0].contourf(r_array, z_array, np.log10(rho_disk), cmap=cmap, levels=levels_rho, extend='both')
    cbar_rho_disk = fig.colorbar(rho_disk_plot, ax=ax[0], orientation='vertical')
    ax[0].set_title('$\\rho_\mathrm{disk}$', fontsize=14)
    
    rho_wind_plot = ax[1].contourf(r_array, z_array, np.log10(rho_wind), cmap=cmap, levels=levels_rho, extend='both')
    cbar_rho_wind = fig.colorbar(rho_wind_plot, ax=ax[1], orientation='vertical')
    ax[1].set_title('$\\rho_\mathrm{wind}$', fontsize=14)
    
    rho_plot = ax[2].contourf(r_array, z_array, np.log10(rho), cmap=cmap, levels=levels_rho, extend='both')
    cbar_rho = fig.colorbar(rho_plot, ax=ax[2], orientation='vertical')
    ax[2].set_title('$\\rho$', fontsize=14)
    
    tau_plot = ax[3].contourf(r_array, z_array, np.log10(tau_FUV), cmap=cmap, levels=np.arange(-2, 4.1, 0.1), extend='both')
    cbar_tau = fig.colorbar(tau_plot, ax=ax[3], orientation='vertical')
    ax[3].set_title('$\\tau_\mathrm{FUV}$', fontsize=14)

    tau1_contour = ax[3].contour(r_array, z_array, tau_FUV, levels=[tau_surface*3.02], colors='red', linewidths=1.5, linestyles=':')
        
    fuv_plot = ax[4].contourf(r_array, z_array, np.log10(FUV_field_attenuated), cmap=cmap, levels=np.arange(0.1, 2, 0.01), extend='both')
    cbar_fuv = fig.colorbar(fuv_plot, ax=ax[4], orientation='vertical')
    ax[4].set_title('FUV', fontsize=14)
    
    fuv1_contour = ax[4].contour(r_array, z_array, tau_FUV, levels=[1.0], colors='white', linewidths=1.5, linestyles=':')
    
    T_plot = ax[5].contourf(r_array, z_array, T, cmap=cmap, levels=np.logspace(0,np.log10(40),25), extend='both')
    cbar_T = fig.colorbar(T_plot, ax=ax[5], orientation='vertical')
    ax[5].set_title('T$_\mathrm{gas}$', fontsize=14)
    
    for i in range(0,6):
        ax[i].set_xlim(0, xlim)
        ax[i].set_ylim(0, ylim)
        
    plt.tight_layout()
    if savefig == True:
        plt.savefig(f'/Users/luke/Documents/QMUL/{figname}.png', dpi=150)
        
        
        
def plot_density(figname, anno, r_array, z_array, rho_disk, xlim, ylim, savefig=True):

    cmap = 'Spectral_r'
        
    levels_rho = np.arange(-20, -11, 0.2)

    plt.contourf(r_array, z_array, np.log10(rho_disk), cmap=cmap, levels=levels_rho, extend='both')
    plt.colorbar()

    plt.annotate(f'Iteration {anno}', (500, 600), c='white', fontsize=14)

    plt.xlim(0, 700)
    plt.ylim(0, 700)

    plt.tight_layout()
    
    if savefig == True:
        plt.savefig(f'/Users/luke/Documents/QMUL/{figname}.png', dpi=150)
    plt.show()
    
    
def smooth_temperature_profile(z_norm, T_midplane, T_surface, profile_type="tanh", k=2):
    """
    Create a smooth temperature profile from midplane to surface.
    
    Parameters:
    z_norm: normalized height (0 at midplane, 1 at tau surface)
    T_midplane: temperature at the midplane
    T_surface: temperature at the tau surface (TPDR)
    profile_type: type of smooth function to use
    """
    
    if profile_type == "tanh":
        # Hyperbolic tangent provides very smooth transition
        # The factor 2 controls the steepness - larger values make it more linear
        smooth_factor = 0.5 * (1 + np.tanh(k * (z_norm - 0.5)))
        return T_midplane + (T_surface - T_midplane) * smooth_factor
    
    elif profile_type == "exponential":
        # Exponential approach to surface temperature
        # This mimics physical heating that becomes more efficient at lower optical depths
        if z_norm == 0:
            return T_midplane
        else:
            alpha = 2.0  # Controls how quickly temperature rises with height
            smooth_factor = (1 - np.exp(-alpha * z_norm)) / (1 - np.exp(-alpha))
            return T_midplane + (T_surface - T_midplane) * smooth_factor
    
    elif profile_type == "cubic":
        # Cubic polynomial with zero derivatives at endpoints
        # This ensures smooth connection to constant regions above/below
        smooth_factor = 3 * z_norm**2 - 2 * z_norm**3
        return T_midplane + (T_surface - T_midplane) * smooth_factor
    
    elif profile_type == "linear":
        # Linear interpolation (your original approach)
        # Simple straight-line transition from midplane to surface
        return T_midplane + (T_surface - T_midplane) * z_norm
    
    else:  # default to linear if unknown type specified
        return T_midplane + (T_surface - T_midplane) * z_norm
    
    


####################################################
# DETERMINE FRACTION OF TOTAL LUMINOSITY IN FUV    #
# - Uses MLR, MRR, and MTR from Eker et al. (2018) #
####################################################


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
    
    
    data = np.loadtxt('FRIEDV2_ALL_fPAH1p0_growth.dat', skiprows=1, delimiter=',')
    
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

def smooth_cosine(x, x1, steepness):
    # Handle edge case: if x1=0, return 1.0 (no tapering)
    if x1 == 0:
        return 1.0
    # Clamp to [0, x1]
    x = np.clip(x, 0, x1)
    # Apply adjustable exponent to control curve shape
    t = (x / x1) ** steepness
    return (1 - np.cos(np.pi * t)) / 2