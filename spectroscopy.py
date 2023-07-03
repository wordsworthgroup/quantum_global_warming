import numpy as np

# fundamental constants
N_avo     = 6.0221409e23     # Avogadro's number []
hplanck   = 6.62607015e-34   # Planck's constant [J.s]
kboltz    = 1.380649e-23     # Boltzmann's constant [J K^-2]
c_light   = 299792458.       # speed of light in vacuum [m s^-1]
epsilon_0 = 8.85418e-12      # permittivity of free space [F m^-1]
stefan    = 5.67e-8          # Stefan-Boltzmann constant [W m^-2 K^-4]

# atomic masses
u_mas    = 1.660539e-27      # atomic mass [kg]
m_oxygen = 2.656e-26         # mass of oxygen atom [kg]
m_carbon = 1.994e-26         # mass of carbon atom [kg]

# CO2 molecule parameters
a_e          = 116e-12                     # C-O bond length [m]
I_co2        = 2*m_oxygen*a_e**2           # rotational moment of inertia of CO2 [kg.m2]
B_co2        = hplanck/(8*np.pi**2*I_co2)  # rotational constant [Hz]
trans_dipole = 3.35e-31                    # transition dipole [D]
k_delta      = 7.7e-19                     # bending force constant of CO2 molecule [N.m]
b            = 2.14e12                     # nonlinear coupling constant [Hz]
Delta_0      = 0.5e12                      # Zeroth order level separation [Hz]
Delta_F      = np.sqrt(Delta_0**2+2*b**2)  # Fermi resonance frequency [Hz]

# thermodynamic properties of atmosphere
Ts        = 288    # surface temperature [K]
Tt        = 216.65 # tropopause temperature [K]
T0        = Ts     # reference temperature [K]
p0        = 1e5    # reference pressure [Pa]
T0_hitran = 296.   # HITRAN reference temperature [K]
T_spec    = T0     # temperature to calculate spectra at [K]
p_spec    = p0     # pressure to calculate spectra at [Pa]

# spectral grid 
nS       = 500000                     # number of spectral intervals []
nu_a     = np.linspace(15,25,nS)*1e12 # frequency [Hz]
nu_THz_a = nu_a/1e12                  # frequency [THz]

# plotting variables
nu_minA  = 15   # minimum frequency in plot A [THz] 
nu_maxA  = 25   # maximum frequency in plot A [THz]
sig_minA = 1e-29 # minimum cross-section in plot A [m2 molecule-1]
sig_maxA = 1e-20  # maximum cross-section in plot A [m2 molecule-1]
nu_minB  = 19   # minimum frequency in plot B [THz] 
nu_maxB  = 21   # maximum frequency in plot B [THz]
sig_minB = 1e-26 # minimum cross-section in plot B [m2 molecule-1]
sig_maxB = 1e-20  # maximum cross-section in plot B [m2 molecule-1]

# functions to interconvert between frequency and wavenumber

def hertz_to_wavenum(Hz):
    "Convert frequency [Hz] to wavenumber [cm^-1]"
    return(1.e-2*Hz/c_light)

def wavenum_to_hertz(wavenum):
    "Convert wavenumber [cm^-1] to frequency [Hz]"
    return(wavenum*c_light*1.e2)

def terahertz_to_wavenum(Hz):
    "Convert frequency [THz] to wavenumber [cm^-1]"
    return(1e12*1.e-2*Hz/c_light)

def normal_mode_2():
    "Calculate frequency of the second normal mode of CO2 (bending mode) [Hz]"
    nu2 = (1/(2*np.pi)*np.sqrt((2/m_oxygen)*(1+2*m_oxygen/m_carbon)*(k_delta/a_e**2)))
    return nu2

def planck_spec_irradiance(nu, T):
    """
    Calculate Planck spectral irradiance [W m-2 sr^-1 Hz-1].
    
    Input
        nu: frequency [Hz]
        T: temperature [K]    
    """
    x = hplanck*nu/(kboltz*T)
    A = (2*hplanck*nu**3)/c_light**2
    B_nu = A/(np.exp(x) - 1)
    return B_nu
