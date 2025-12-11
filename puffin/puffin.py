import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import helpers

def DiskModel1D(m_star, r_d, sigma_1au, FFUV_G0, n_points=200, gridsize=None, m_dot=None, gamma=None, p=0.2, q=0.4):
    
    """
    Calculate 1D radial density structure for externally FUV-irradiated protoplanetary disks.
    
    This function computes the gas density profile for a disk undergoing external 
    photoevaporation. The model includes three components: the underlying disk structure,
    an outflowing wind, and a transitional plateau region that smoothly connects them.
    
    Parameters
    ----------
    m_star : float
        Stellar mass in solar masses.
    r_d : float
        Characteristic disk radius in AU.
    sigma_1au : float
        Surface density at 1 AU in g/cm^2.
    FFUV_G0 : float
        External FUV field strength in Habing units (G0).
    n_points : int, optional
        Number of radial grid points. Default is 200.
    gridsize : float, optional
        Outer radius of the computational grid in AU. Default is 8 * r_d.
    m_dot : float, optional
        Mass loss rate in M_sun/yr. If None, interpolated from FRIED grid.
    gamma : float, optional
        Exponential cutoff parameter for surface density profile. If None, 
        calculated from scaling relation.
    p : float, optional
        Power law index for plateau transition (r_d**p). Default is 0.2.
    q : float, optional
        Power law index for plateau transition ((FFUV_G0/100)**q). Default is 0.4.
    
    Returns
    -------
    r_array : numpy.ndarray
        Radial grid in AU (logarithmically spaced).
    rho : numpy.ndarray
        Total gas density in g/cm^3 at each radial point.
    
    Notes
    -----
    The density structure consists of:
    - Disk component: Power-law surface density with exponential taper
    - Wind component: Spherical outflow driven by photoevaporation (r >= r_d)
    - Plateau component: Smooth transition between disk and wind (r > r_d)
    
    The final density at each radius is the maximum of the three components.
    """
    
    ###############
    # Setup model #
    ###############
    
    # Grid
    if gridsize is None:
        gridsize = r_d * 8
    r_array = np.logspace(np.log10(0.1), np.log10(gridsize), n_points)
    
    # Interpolate M_dot
    if m_dot is None:
        m_dot_log = helpers.interpolate_mdot(m_star, r_d, sigma_1au, FFUV_G0)
        m_dot = 10.0**(m_dot_log)
    
    # Constants
    mu_d   = 2.3
    mu_PDR = 1.3
    big_G  = 6.67259e-8
    k_B    = 1.380626e-16
    m_sun  = 1.9891e33
    au_cm  = 1.495979e13
    m_p    = 1.67e-24
    yr_sec = 3.15576e+07
    
    # Disk temperature
    t_1au = 150.0 * (m_star)**0.25 
    
    # Wind parameters
    t0_PDR = 200.0
    t_PDR  = min(max(t0_PDR * (FFUV_G0/1000.0)**0.2, 10), 3000.0)
    cs_PDR = np.sqrt((k_B * t_PDR) / (mu_PDR * m_p))
    
    # Calculate mathcal_F at disk edge
    t_d       = max(t_1au * r_d**(-0.5), 10)
    cs_d      = np.sqrt((k_B * t_d) / (mu_d * m_p))
    omega     = np.sqrt((big_G * m_star * m_sun) / (r_d*au_cm)**3)
    h_d       = cs_d / omega
    mathcal_F = h_d / np.sqrt((r_d*au_cm)**2 + h_d**2)
    
    # Initialise arrays
    sigma       = np.zeros(len(r_array))
    temp        = np.zeros(len(r_array))
    rho_disk    = np.zeros(len(r_array))
    rho_wind    = np.zeros(len(r_array))
    rho_plateau = np.zeros(len(r_array))
    rho         = np.zeros(len(r_array))
    
    # Index for r_d
    idx_rd = (np.abs(r_array - r_d)).argmin()
    
    
    ###############################
    # Calculate density structure #
    ###############################
    
    if gamma is None:
        gamma = (0.77) * (m_star**0.10) * (r_d)**0.78 * ((FFUV_G0)/100)**0.07
        
    for i in range(len(r_array)):
        # Disk component
        sigma[i]    = sigma_1au * (r_array[i]**(-1.0)) * np.exp(-((r_array[i])/(r_d*1.1))**(gamma)) 
        temp[i]     = max(t_1au * (r_array[i]**(-0.5)), 10)
        cs          = np.sqrt((k_B * temp[i]) / (mu_d * m_p))
        omega       = np.sqrt((big_G * m_star * m_sun) / (r_array[i]*au_cm)**3)
        H           = cs / omega
        rho_disk[i] = sigma[i] / H
        
        # Wind component
        if r_array[i] >= r_d:
            rho_wind[i] = (m_dot * m_sun) / (4.0 * np.pi * r_array[i]**2 * au_cm**2 * yr_sec * mathcal_F * cs_PDR)
        
        # Plateau component
        if r_array[i] > r_d:
            x_norm         = np.log10(r_array[i]/r_d)
            lamda          = r_d**p + (FFUV_G0/100)**q
            f              = (1 - np.exp(-lamda * x_norm)) / (1 - np.exp(-lamda))
            rho_plateau[i] = rho_disk[idx_rd] * (rho_wind[i] / rho_disk[idx_rd])**f
    
    # Total density (maximum of all components)
    for j in range(len(r_array)):
        if r_array[j] >= r_d:
            rho[j] = np.maximum.reduce([rho_disk[j], rho_wind[j], rho_plateau[j]])
        else:
            rho[j] = rho_disk[j]
    
    return r_array, rho




def DiskModel2D(m_star, r_d, sigma_1au, FFUV_G0, n_points=1000, gridsize=None, m_dot=None, gamma=None, p=0.2, q=0.4, N_ITER=20):
    
    """
    Calculate 2D density structure for externally FUV-irradiated protoplanetary disks.
    
    This function computes the gas density profile for a disk undergoing external 
    photoevaporation. The model solves for vertical hydrostatic equilibrium iteratively 
    and includes multiple components: a hydrostatic disk, a spherically diverging 
    photoevaporative wind launched from the disk surface, and a 'bowl'-shaped transition 
    region that smoothly connects the disk to the wind.
    
    Parameters
    ----------
    m_star : float
        Stellar mass in solar masses.
    FFUV_G0 : float
        External FUV field strength in Habing units (G0).
    r_d : float
        Characteristic disk radius (gravitational radius) in AU.
    sigma_1au : float
        Surface density at 1 AU in g/cm^2.
    n_points : int, optional
        Number of grid points in both radial and vertical directions. Default is 1000.
        High resolution (≥1000) is recommended for accurate HSE solutions.
    gridsize : float, optional
        Outer radius of the computational grid in AU. Default is min(8 * r_d, 800).
    m_dot : float, optional
        Mass loss rate in M_sun/yr. If None, interpolated from FRIED grid with 
        a factor of 2 applied for 2D geometry.
    gamma : float, optional
        Exponential cutoff parameter for surface density profile. If None, 
        calculated from scaling relation.
    p : float, optional
        Power law index for bowl transition (r_d component). Default is 0.2.
    q : float, optional
        Power law index for bowl transition (FFUV component). Default is 0.4.
    N_ITER : int, optional
        Number of hydrostatic equilibrium iterations. Default is 20.
    
    Returns
    -------
    r_array : numpy.ndarray
        Radial grid in AU (logarithmically spaced).
    z_array : numpy.ndarray
        Vertical grid in AU (logarithmically spaced, includes midplane z=0).
    rho_total : numpy.ndarray
        Total gas density in g/cm^3 at each (z, r) point. Shape is (nz, nr).
        
    Returns (on error)
    ------------------
    result : str
        Error message if model fails (e.g., "ABORTED: Disk mass < 1e-5 M_sun").
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
      ∝ 1/r² set by mass loss rate and sound speed
    - Bowl region: Smooth transition from disk to spherical wind using exponential 
      blending function, tapered inside r_d
    - Wind is smoothly tapered at edges using cosine function
    
    The final density at each position is the maximum of the disk and wind components.
    
    **Hydrostatic Equilibrium:**
    The code iteratively solves the vertical structure by:
    1. Computing optical depths (vertical, radial inward/outward)
    2. Calculating FUV attenuation
    3. Determining temperature structure (disk + PDR)
    4. Solving d ln ρ / dz = -Ω²z/c_s² - (1/T) dT/dz
    5. Renormalizing to match Σ(r) at each radius
    
    The model aborts if:
    - Disk mass < 1e-5 M_sun
    - Disk is optically thin at r_d (tau < 1)
    """
    
    # Model name
    model_name = f'{m_star}_{FFUV_G0}_{r_d}_{sigma_1au}'
    
    # Interpolate M_dot
    if m_dot is None:
        factor_2d = 2
        m_dot_log = helpers.interpolate_mdot(m_star, r_d, sigma_1au, FFUV_G0)
        m_dot     = factor_2d * 10**m_dot_log          
    
    # Star/disk properties
    L_star, R_star, T_eff = helpers.get_stellar_properties(m_star)
    fuv_fraction          = helpers.fuv_fraction(T_eff)
    FUV_star              = fuv_fraction * L_star * 3.828e33   

    t_1au  = 80.0 * (m_star)**0.25    # from DALI models
    m_disk = helpers.calculate_disk_mass(sigma_1au, r_max=r_d)  
    mid_j  = 0

    if m_disk < 1e-5:
        result = 'ABORTED: Disk mass < 1e-5 M_sun'
        return result, model_name


    ########
    # GRID #
    ########

    if gridsize is None:
        gridsize = r_d * 8
    
    if n_points < 1000:
        print('WARNING: n_points < 1000')
        print('HSE solve requires higher resolution for best results')
        print('Regrid your model afterwards if desired')

    r_array    = np.logspace(np.log10(0.1), np.log10(gridsize), n_points)     # AU, ascending
    z_positive = np.logspace(np.log10(0.01), np.log10(gridsize), n_points-1)  # AU (no zero)
    z_array    = np.concatenate(([0.0], z_positive))                          # include midplane as z=0

    nr = len(r_array)
    nz = len(z_array)


    #############
    # CONSTANTS #
    #############

    mu_d   = 2.3    
    mu_PDR = 1.3        
    big_G  = 6.67259e-8
    k_B    = 1.380626e-16
    m_sun  = 1.9891e33
    au     = 1.495979e13
    m_p    = 1.6726219e-24
    yr_sec = 3.15576e7

    sigma_FUV_wind = 2.7e-23                          # cm^2 per *particle* 
    sigma_FUV_disk = 8.0e-22                          # cm^2 per *particle* 
    kappa_FUV_wind = sigma_FUV_wind / (mu_PDR * m_p)  # cm^2/g
    kappa_FUV_disk = sigma_FUV_disk / (mu_d   * m_p)  # cm^2/g
    
    Rfocal = 0.5 * r_d


    ##################
    # PDR PROPERTIES #
    ##################

    T0_PDR = 200.0
    TPDR   = float(np.clip(T0_PDR * (FFUV_G0/1000.0)**0.2, 10.0, 3000.0))
    csPDR  = np.sqrt((k_B * TPDR) / (mu_PDR * m_p)) 

    tau_surface = 1.0      # value of tau to define disk surface
    k_smooth    = 1.75     # determines steepness of vertical temperature gradient 


    ###################
    # ALLOCATE ARRAYS #
    ###################

    sigma_r             = np.zeros(nr)          
    rho_disk            = np.zeros((nz, nr))
    rho_disk_pre_hse    = np.zeros((nz, nr))
    rho_wind_exp        = np.zeros((nz, nr))
    rho_wind_sph        = np.zeros((nz, nr))
    rho_bowl            = np.zeros((nz, nr))
    T                   = np.zeros((nz, nr))
    H                   = np.zeros((nz, nr))    
    FUV_field_g0        = np.zeros((nz, nr))
    tau_FUV_vertical    = np.zeros((nz, nr))
    tau_FUV_radial_in   = np.zeros((nz, nr))
    tau_FUV_radial_out  = np.zeros((nz, nr))
    tau_FUV_min         = np.zeros((nz, nr))
    dT_dz               = np.zeros((nz, nr))
    FUV_field_att_g0    = np.zeros((nz, nr))   
    dist_array_tau      = np.zeros((nz, nr)) 
    dist_array_bowl     = np.zeros((nz, nr))
    TPDR_FUV            = np.zeros((nz, nr)) 
    f_scale             = np.zeros((nz, nr)) 


    #############
    # FUNCTIONS #
    #############

    # Compute optical depths
    def compute_optical_depths(rho_field, kappa):
        """Return (tau_vert, tau_in, tau_out, tau_min) with shapes (nz, nr)."""
        tau_v  = np.zeros_like(rho_field)
        tau_in = np.zeros_like(rho_field)
        tau_out= np.zeros_like(rho_field)

        # 1) Vertical (integrate from top down)
        for i in range(nr):
            this_col = 0.0
            for j in range(nz-2, -1, -1):  # from top-1 down to 0
                dz_cm   = (z_array[j+1] - z_array[j]) * au
                rho_avg = 0.5 * (rho_field[j+1, i] + rho_field[j, i])
                this_col += rho_avg * dz_cm * kappa
                tau_v[j, i] = this_col
            tau_v[-1, i] = 0.0  # top boundary

        # 2) Radial inward (integrate from outer edge inward)
        for j in range(nz):
            this_col = 0.0
            for i in range(nr-2, -1, -1):
                dr_cm  = (r_array[i+1] - r_array[i]) * au
                rho_avg= 0.5 * (rho_field[j, i+1] + rho_field[j, i])
                this_col += rho_avg * dr_cm * kappa
                tau_in[j, i] = this_col
            tau_in[j, -1] = 0.0

        # 3) Radial outward (integrate from inner edge outward)
        for j in range(nz):
            this_col = 0.0
            for i in range(1, nr):
                dr_cm  = (r_array[i] - r_array[i-1]) * au
                rho_avg= 0.5 * (rho_field[j, i] + rho_field[j, i-1])
                this_col += rho_avg * dr_cm * kappa
                tau_out[j, i] = this_col
            tau_out[j, 0] = 0.0

        tau_min = np.minimum.reduce([tau_v, tau_in, tau_out])
        return tau_v, tau_in, tau_out, tau_min





    #####################
    # INTIAL CONDITIONS #
    #####################

    if gamma is None:
        gamma = (0.77) * (m_star**0.10) * (r_d)**0.78 * ((FFUV_G0)/100)**0.07

    # Disk component
    sigma_r = sigma_1au * (r_array**-1.0) * np.exp(-(r_array/r_d)**gamma)
    for i in range(nr):
        Omega = np.sqrt((big_G * m_star * m_sun) / (r_array[i]*au)**3)
        T_mid = max(t_1au * (r_array[i]**-0.5), 10.0)
        cs    = np.sqrt((k_B * T_mid) / (mu_d * m_p))
        H_i   = cs / Omega  
        rho0  = sigma_r[i] / (np.sqrt(2.0*np.pi) * H_i)

        for j in range(nz):
            zcm = z_array[j] * au
            H[j, i] = H_i  
            rho_disk[j, i] = rho0 * np.exp(-0.5 * (zcm/H_i)**2)

    # Stellar & external FUV field (unattenuated)
    for i in range(nr):
        for j in range(nz):
            zcm = z_array[j] * au
            distance_cm = np.sqrt((r_array[i]*au)**2 + zcm**2)
            FUV_star_local = FUV_star / (4.0 * np.pi * distance_cm**2)  
            FUV_field_g0[j, i] = max(FUV_star_local / 1.6e-3 + FFUV_G0, 1e-30)


    # Properties at r_d
    idx_rd = np.searchsorted(r_array, r_d, side="right") - 1
    rho_rd = rho_disk[0, idx_rd]
    


    ###############
    # ITERATE HSE #
    ###############
    
    for it in range(N_ITER):
        
        # ----------------- #
        # 1. OPTICAL DEPTH  #
        # ----------------- #
        
        # Compute optical depths
        tau_FUV_vertical, tau_FUV_radial_in, tau_FUV_radial_out, tau_FUV_min = [val for val in compute_optical_depths(rho_disk, kappa_FUV_disk)]

        # Determine surface
        tau_midplane = tau_FUV_min[mid_j, :]  
        idx_tau_outside_in = np.where(tau_midplane[::-1] >= 1)[0]
        if idx_tau_outside_in.size > 0:
            idx_tau1_midplane = len(tau_midplane) - 1 - idx_tau_outside_in[0]
        else:
            idx_tau1_midplane = None
        
        if idx_tau1_midplane is None: 
            result = "ABORTED: Disk optically thin (tau_rd=0)"
            return result, model_name
        elif r_array[idx_tau1_midplane] < r_array[idx_rd]:
            result = "ABORTED: Disk optically thin (tau=1 inside r_d)"
            return result, model_name
        else:
            surface_choice = 'Tau'


        # ------------- #
        # 2. FUV FIELD  #
        # ------------- #
        
        FUV_field_att_g0[:, :] = FUV_field_g0 * np.exp(-tau_FUV_min)
        
        
        # ------------------- #
        # 3. TEMPERATURE PDR  #
        # ------------------- #
        
        for i in range(nr):
            for j in range(nz):
                # TPDR_FUV[j,i] = float(np.clip(T0_PDR * (FUV_field_att_g0[j,i]/1000.0)**0.2, 10.0, 3000.0))    
                TPDR_FUV[j,i] = float(np.clip(T0_PDR * (FFUV_G0/1000.0)**0.2, 10.0, 3000.0))                   # no stellar heating

        # -------------------- #
        # 4. TEMPERATURE DISK  #
        # -------------------- #
        
        # Recalculate rho_rd for this iteration
        rho_rd = rho_disk[0, idx_rd]
        
        # Compute temperature structure
        for i in range(nr):
            
            # Find tau=1 surface
            z_surface_tau = None
            for j in range(nz-1, -1, -1):
                if tau_FUV_vertical[j, i] >= tau_surface:
                    z_surface_tau = z_array[j]
                    break

            z_surface = z_surface_tau

            # 2D temperature profile
            
            T_midplane = max(t_1au * (r_array[i]**-0.5), 10.0)

            if z_surface is None:
                for j in range(nz):
                    T[j,i] = TPDR_FUV[j,i]
            elif z_surface > 0.0:
                for j in range(nz):
                    if z_array[j] <= z_surface:
                        z_norm = z_array[j] / z_surface
                        T[j,i] = helpers.smooth_temperature_profile(z_norm, T_midplane, TPDR_FUV[j,i], "tanh", k=k_smooth)
                    else:
                        T[j,i] = TPDR_FUV[j,i]
            elif z_surface == 0:
                for j in range(nz):
                    T[j,i] = TPDR_FUV[j,i]
            else:
                for j in range(nz):
                    print(f'No temperature set at r={r_array[i]}au z={z_array[j]}au')


        # --------------------------------- #
        # 5. SOLVE HYDROSTATIC EQUILIBRIUM  #
        # --------------------------------- #

        # Gradients dT/dz (in K/AU)
        
        for i in range(nr):
            for j in range(nz):
                if j == 0:
                    dT_dz[j, i] = (T[j+1, i] - T[j, i]) / (z_array[j+1] - z_array[j])
                elif j == nz-1:
                    dT_dz[j, i] = (T[j, i] - T[j-1, i]) / (z_array[j] - z_array[j-1])
                else:
                    dT_dz[j, i] = (T[j+1, i] - T[j-1, i]) / (z_array[j+1] - z_array[j-1])


        # Solve HSE upward for ln(rho) then renormalize to match Sigma(r)

        for i in range(nr):
        
            Omega = np.sqrt((big_G * m_star * m_sun) / (r_array[i]*au)**3)

            # Midplane boundary condition: keep current rho at z=0
            rho_mid = max(rho_disk[mid_j, i], 1e-40)
            ln_rho  = np.log(rho_mid)

            rho_disk[mid_j, i] = rho_mid  # store midplane explicitly

            for j in range(1, nz):
                dz_cm   = (z_array[j] - z_array[j-1]) * au
                z_cm    = (z_array[j-1]) * au
                T_loc   = max(T[j-1, i], 10.0)
                cs2     = (k_B * T_loc) / (mu_d * m_p)
                dTdz_cm = (dT_dz[j-1, i] / au)  # convert K/AU -> K/cm

                # d ln ρ / dz = - Ω^2 z / c_s^2 - (1/T) dT/dz
                rhs = -(Omega**2 * z_cm / cs2) - (dTdz_cm / T_loc)

                ln_rho += rhs * dz_cm
                rho_disk[j, i] = np.exp(ln_rho)

            # Renormalize column to match Sigma(r) at this radius
            col = np.trapz(rho_disk[:, i], z_array*au)  # g/cm^2
            if col > 0.0:
                rho_disk[:, i] *= (sigma_r[i] / col)


    rho_tau1_midplane = rho_disk[0, idx_tau1_midplane]



    ##################
    # SPHERICAL WIND #
    ##################
        
    # Recalculate tau (disk + expo)
    rho_new_tau = np.where(rho_wind_exp > 1e-30, rho_wind_exp, rho_disk)

    # base_mask_tau = (tau_FUV_min >= 1)
    base_mask_tau = (rho_disk >= rho_tau1_midplane)


    # Step 1: coordinates of wind-base "solid" region (interior to the chosen surface)
    surface_indices_tau = np.argwhere(base_mask_tau)
    if surface_indices_tau.size > 0:
        surface_coords_tau = np.column_stack((r_array[surface_indices_tau[:,1]], z_array[surface_indices_tau[:,0]])) # z
        surface_rhos_tau   = rho_new_tau[surface_indices_tau[:,0], surface_indices_tau[:,1]]

        # Step 2: KD-tree for nearest distance to the base region
        tree_tau = cKDTree(surface_coords_tau)

        # Step 3: For each cell, distance to base & exponential decay
        for i in range(len(r_array)):
            for j in range(len(z_array)):
                r = r_array[i]
                z = z_array[j]

                dist_tau, idx_tau = tree_tau.query([r, z])
                dist_array_tau[j, i] = dist_tau
                rho_base_tau = surface_rhos_tau[idx_tau]
                
                total_dist_tau = dist_array_tau[j,i] + Rfocal
                
                if dist_tau > 0:
                    rho_wind_sph[j,i] =  ((m_dot * m_sun) / (4.0 * np.pi * (total_dist_tau)**2 * au**2 * yr_sec * csPDR))
                else:
                    rho_wind_sph[j, i] = 1e-30   # ie. use rho_disk inside tau=1 surface



    ################
    # Rescale Wind #
    ################

    surface_rho_min = np.nanmin(surface_rhos_tau[surface_rhos_tau > 1e-30])
    max_wind_sph = np.nanmax(rho_wind_sph)

    wind_scaling = max_wind_sph/surface_rho_min

    if wind_scaling > 1:
        rho_wind_sph = rho_wind_sph / wind_scaling
        print(f'WIND SCALED BY {wind_scaling}')
    
   
    
    
    ########
    # BOWL #
    ########

    base_mask_bowl = (rho_disk >= rho_tau1_midplane)

    surface_indices_bowl = np.argwhere(base_mask_bowl)
    if surface_indices_bowl.size > 0:
        surface_coords_bowl = np.column_stack((r_array[surface_indices_bowl[:,1]], 
                                            z_array[surface_indices_bowl[:,0]]))
        surface_rhos_bowl = rho_disk[surface_indices_bowl[:,0], surface_indices_bowl[:,1]]
        
        tree_bowl = cKDTree(surface_coords_bowl)
        
        for i in range(len(r_array)):
            for j in range(len(z_array)):
                
                r = r_array[i]
                z = z_array[j]
                
                dist_bowl, idx_bowl = tree_bowl.query([r_array[i], z_array[j]])
                dist_array_bowl[j, i] = dist_bowl
                
                rho_base_bowl = surface_rhos_tau[idx_bowl]
                
                if dist_bowl > 0:
                    
                    # Bowl smoothing function (1D equivalent)
                    x_norm = np.log10(1+dist_bowl / r_d)
                    k = r_d**p + (FFUV_G0/100)**q
                    f = (1 - np.exp(-k * x_norm)) / (1 - np.exp(-k))
                    
                    # Blend surface density → spherical wind
                    rho_bowl[j, i] = rho_base_bowl * (rho_wind_sph[j, i] / rho_base_bowl)**f
                else:
                    rho_bowl[j, i] = 1e-30  
                    
                # Taper inside r_d
                if r_array[i] < r_d:                 
                    scaling = (r_array[i]/r_d)
                    rho_bowl[j, i] = rho_bowl[j, i] * scaling


    ##############
    # Taper Wind #
    ##############
    
    rho_wind = np.maximum(rho_wind_sph, rho_bowl)
      
    for i in range(len(r_array)):
        for j in range(len(z_array)):
            
            if z_array[j] >= r_d/2:
                edge = r_d
            else:
                edge = r_d * (z_array[j]/(r_d/2))
            
            f_scale[j,i]  = helpers.smooth_cosine(r_array[i], edge, steepness=2)
            rho_wind[j,i] = rho_wind[j,i] * f_scale[j,i]

    
    #############
    # Final rho #
    #############
    
    rho_disk[rho_disk < rho_tau1_midplane] = 1e-30
    
    rho_total = np.maximum(rho_disk, rho_wind)

    return r_array, z_array, rho_total