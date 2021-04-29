import numpy as np
import scipy.constants as cte
import scipy.optimize as opt

# scale variables
def rho_s( B_0, m_i, T_e) :
    """ Spatial scale [m]: ion Larmor radius at electron temperature"""
    return  np.sqrt( T_e *cte.eV * m_i)  / cte.e / B_0

def omega_0_inv( B_0, m_i):
    """ Time scale [s]: inverse ion Larmor frequency"""
    return m_i / cte.e / B_0

def c_s( m_i, T_e) :
    """ Velocity scale [m/s]: Ion sound speed at electron temperature"""
    return np.sqrt( T_e*cte.eV/m_i)

def Phi_0( T_e) :
    """Scale of potential [eV]"""
    return T_e*cte.eV/cte.e

def Psip_0(B_0, m_i, T_e, R):
    """ Scale of magnetic flux []"""
    return B_0*rho_s( B_0, m_i, T_e)*R

def lambda_D( T_e, n_0):
    """Debye length [m]"""
    return np.sqrt(cte.epsilon_0*T_e*cte.eV/n_0/1e19/cte.e**2)

# independent numerical variables
def mu_e ( m_i) : # proton mass, deuteron mass, triton mass
    """ electron mass ratio"""
    return -cte.m_e/ m_i

def tau_i ( T_e, T_i) :
    """ ion to electron temperature ratio"""
    return T_i / T_e

def beta( n_0, T_e, B_0) :  # n_0 in 1e-19, T_e in eV, B_0 in T
    """ electron plasma beta"""
    return n_0*1e19*T_e*cte.eV / (B_0**2 / cte.mu_0)

def resistivity( n_0, T_e, B_0):
    """plasma resistivity"""
    return 0.51*np.sqrt( 2*cte.m_e)*cte.e**3*10/ 12/np.pi**(3/2)/\
            cte.epsilon_0**2 * n_0*1e19 / B_0 / (T_e*cte.eV)**(3/2)

def R_0 ( B_0, m_i, T_e, R) :
    """machine radius relative to Larmor radius"""
    return R / rho_s(B_0, m_i, T_e)

def epsilon_D ( n_0, B_0, m_i) :
    """Square Debye length relative to Larmor radius"""
    return cte.epsilon_0/n_0/1e19/m_i*B_0**2

# dependent numerical variables

def viscosity_e( resistivity) :
    """ parallel electron viscosity"""
    return 0.37/resistivity

def viscosity_i( resistivity, mu_e) :
    """ parallel ion viscosity"""
    #return 0.96/np.sqrt( m_i)/cte.e**3/10* 12*np.pi**(3/2)*\
    #        cte.epsilon_0**2 * B_0 * (T_i*cte.eV)**(3/2) /(n_0*1e19)
    return 0.69/resistivity*np.sqrt( np.abs(mu_e))

# Miscellaneous ( epsilon_a_inv, q)
def minor_radius( epsilon_a_inv, R_0) :
    """as a function of inverse aspect ratio and R_0 (any unit)"""
    return epsilon_a_inv*R_0

def lparallel( q, R_0):
    return q*R_0

#def alpha_MHD( q, epsilon_a_inv, beta, DeltaN):
#    """From Eich2020 paper Eq. (3)"""
#    return q**2/epsilon_a_inv*beta*DeltaN

def Myra (resistivity, lparallel):
    return resistivity * lparallel

def Knudsen (mu_e, resistivity, lparallel):
    """ 0.51 sqrt(mu_e)/Myra"""
    return 0.51*sqrt( abs( mu_e))/Myra( resistivity, lparallel)

def pref(R_0, beta):
    """ presssure reference p = pref * Psip """
    return 1/R_0/beta

def nu_perp(resistivity, mu_e):
    return 0.96/viscosity_i(resistivity, mu_e)

def GAM_inv (R_0) :
    """ inverse GAM frequency """
    return 1.0/R_0

# other scales
def paraview_scale(rho_s):
    return 1.0/rho_s

def GBS_timescale(R_0, omega_0_inv) :
    """ timescale in GBS units """
    return R_0*omega_0_inv

def GBS_resistivity(R_0, resistivity):
    """ resistivity in GBS units """
    return resistivity*R_0

def CFL_diff( R_0, epsilon_a_inv, scale_R, N_z, viscosity_e):
    """ CFL condition due to parallel diffusion """
    return (2.*np.pi/Nz*R_0*(1-scale_R*epsilon_a_inv))**2/viscosity_e


proton_mass   = cte.physical_constants["proton mass"][0]
deuteron_mass = cte.physical_constants["deuteron mass"][0]
triton_mass   = cte.physical_constants["triton mass"][0]

def physical2numerical( physical, numerical):
    """Compute numerical parameters from physical parameters"""
    numerical["mu"] = mu_e( physical["m_i"])
    numerical["tau_i"] = tau_i( physical["T_e"], physical["T_i"])
    numerical["beta"] = beta( physical["n_0"], physical["T_e"], physical["B_0"])
    numerical["resistivity"] = resistivity( physical["n_0"], physical["T_e"], physical["B_0"])
    numerical["R_0"] = R_0( physical["B_0"], physical["m_i"], physical["T_e"], physical["R_0"])
    numerical["epsilon_D"] = epsilon_D( physical["n_0"], physical["B_0"], physical["m_i"])

def numerical2physical( numerical, physical):
    """Invert numerical parameters to physical parameters"""

    physical["m_i"] = -cte.m_e/numerical["mu"]

    if ( ("epsilon_D" in numerical) and (numerical["epsilon_D"] != 0) ) :
        if (not("beta" in numerical) or (numerical["beta"] == 0)):
            raise ValueError("beta must be present if epsilon_D is")
        def to_invert0( x,
                       for_e_D,for_mu_e, for_tau_i, for_beta, for_eta, for_R_0) :
            T_e, T_i, n_0, B_0, R = x
            m_i = -cte.m_e/for_mu_e
            return (tau_i(T_e, T_i) - for_tau_i,
                (beta( n_0, T_e, B_0) - for_beta)/for_beta,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (epsilon_D( n_0, B_0, m_i) - for_e_D)/for_e_D,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        print( "Invert for given numerical parameters")
        physical["T_e"], physical["T_i"], physical["n_0"], physical["B_0"], physical["R_0"]\
        = opt.fsolve( to_invert0, [1,1,1,1,1],args=(
            numerical["epsilon_D"],numerical["mu"],numerical["tau_i"],
            numerical["beta"],numerical["resistivity"],numerical["R_0"]))

    elif (("beta" in numerical) and (numerical["beta"] != 0) ):
        print( "Invert for given R_0")
        def to_invert1( x, R,
              for_mu_e, for_tau_i, for_beta, for_eta, for_R_0) :
            T_e, T_i, n_0, B_0 = x
            m_i = -cte.m_e/for_mu_e
            return (tau_i(T_e, T_i) - for_tau_i,
                (beta( n_0, T_e, B_0) - for_beta)/for_beta,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        physical["T_e"], physical["T_i"], physical["n_0"], physical["B_0"]\
        = opt.fsolve( to_invert1, [1,1,1,1],args=(
            physical["R_0"],numerical["mu"],numerical["tau_i"],
            numerical["beta"],numerical["resistivity"],numerical["R_0"]))
        numerical["epsilon_D"] = epsilon_D(physical["n_0"],
                                        physical["B_0"],
                                        physical["m_i"])
    else:
        print( "Invert for given R_0 and B_0")
        if ((not "R_0" in physical) or (physical["R_0"] == 0) or
            (not "B_0" in physical) or (physical["B_0"] == 0)):
            raise ValueError ( "B_0 and R_0 must be present in physical and be different from zero")
        def to_invert2( x, B_0, R,
                  for_mu_e, for_tau_i, for_eta, for_R_0) :
            T_e, T_i, n_0 = x
            m_i = -cte.m_e/for_mu_e
            return (tau_i(T_e, T_i) - for_tau_i,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        physical["T_e"], physical["T_i"], physical["n_0"]\
        = opt.fsolve( to_invert2, [1,1,1],args=(
            physical["B_0"],physical["R_0"],numerical["mu"],
            numerical["tau_i"],numerical["resistivity"],numerical["R_0"]))
        numerical["beta"] = beta( physical["n_0"],
                                    physical["T_e"],
                                    physical["B_0"])
        numerical["epsilon_D"] = epsilon_D(physical["n_0"],
                                        physical["B_0"],
                                        physical["m_i"])

def load_calibration_default():
    """ generate default feltor input parameters for calibration"""
    inputfile= {
        "grid":
        {
            "n"  : 3,
            "Nx" : 48,
            "Ny" : 96,
            "Nz" : 8,
            "scaleR" : [1.25,1.2],
            "scaleZ" : [2.2, 2.0]
        },
        "timestepper":
        {
            "tableau": "TVB-3-3",
            "dt" : 1e-2
        },
        "regularization":
        {
            "order" : 2,
            "direction": "centered",
            "nu_perp_n" : 1e-5,
            "nu_perp_u" : 1e-5
        },
        "elliptic":
        {
            "stages" : 2,
            "eps_pol"    : [1e-6, 10],
            "eps_gamma" : 1e-7,
            "eps_ampere" : 1e-7,
            "direction" : "centered",
            "jumpfactor" : 1.0
        },
        "FCI":
        {
            "refine" : [1,1],
            "rk4eps" : 1e-6,
            "periodify": True,
            "bc" : "along_field"
        },
        "physical":
        {
            "mu"          : -0.000272121,
            "tau"         : 0.5,
            "beta"        : 1e-5,
            "resistivity" : 1e-4,
            "viscosity" : "braginskii"
        },
        "output":
        {
            "type" : "netcdf",
            "inner_loop" : 2,
            "itstp"  : 2,
            "maxout" : 1,
            "compression" : [1,1]
        },
        "source":
        {
            "type" : "influx",
            "rate" : 2e-3,
            "ntilde": {
                "type" : "turbulence",
                "amplitude"  : 0.1,
                "rk4eps" : 1e-6,
                "refine" : [5,5],
                "revolutions" : 1,
                "sigma_z" : 0.5
            },
            "profile":
            {
                "type" : "aligned",
                "amplitude" : 1.0
            },
            "damping":
            {
                "type": "aligned",
                "alpha" : 0.2,
                "boundary" : 0.55
            }
        },
        "init":
        {
            "type" : "ne",
            "potential" : "zero",
            "ntilde": {
                "type" : "turbulence",
                "amplitude"  : 0.1,
                "rk4eps" : 1e-6,
                "refine" : [5,5],
                "revolutions" : 1,
                "sigma_z" : 0.5
            },
            "profile":
            {
                "type" : "aligned",
                "amplitude" : 1.0
            },
            "damping":
            {
                "type": "aligned",
                "alpha" : 0.2,
                "boundary" : 1.0
            }
        },
        "boundary":
        {
            "wall":
            {
                "type": "sol_pfr",
                "boundary": [1.1,0.998],
                "alpha": [0.10,0.10],
                "modify-B" : True,
                "penalization" : 1e+0,
                "penalize-rhs" : False
            },
            "sheath":
            {
                "type": "bohm",
                "boundary": 0.30,
                "alpha": 0.2,
                "penalization" : 1e+0,
                "penalize-rhs" : False
            },
            "bc" :
            {
                "density" : ["DIR", "DIR"],
                "velocity": ["NEU", "NEU"],
                "potential":["DIR", "DIR"],
                "aparallel":["NEU", "NEU"]
            }
        },
        "flags" : [ "symmetric", "calibrate"],
        "magnetic_field":
        {
            "curvmode" : "toroidal",
            "input" : "params",
        }
    }
    return inputfile

