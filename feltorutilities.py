import numpy as np
import scipy.constants as cte
import scipy.optimize as opt

tasks = {}
# design copied from stackoverflow
task = lambda f: tasks.setdefault( f.__name__, f)

#WE PASS **kwargs TO ALL FUNCTIONS SO THAT WE CAN CALL ALL WITH DICTIONARIES

# identities (not sure we need them any more)
#@task
#def T_e( T_e, **kwargs) :
#    return T_e
#@task
#def T_i( T_i, **kwargs) :
#    return T_i
#@task
#def m_i( m_i, **kwargs) :
#    return m_i
#@task
#def R( R, **kwargs) :
#    return R
#@task
#def n_0( n_0, **kwargs) :
#    return n_0
#@task
#def B_0( B_0, **kwargs) :
#    return B_0


# scale variables
@task
def rho_s( B_0, m_i, T_e, **kwargs) :
    """ Spatial scale [m]: ion Larmor radius at electron temperature"""
    return  np.sqrt( T_e *cte.eV * m_i)  / cte.e / B_0

@task
def omega_0_inv( B_0, m_i, **kwargs):
    """ Time scale [s]: inverse ion Larmor frequency"""
    return m_i / cte.e / B_0

@task
def c_s( m_i, T_e, **kwargs) :
    """ Velocity scale [m/s]: Ion sound speed at electron temperature"""
    return np.sqrt( T_e*cte.eV/m_i)

@task
def Phi_0( T_e, **kwargs) :
    """Scale of potential [eV]"""
    return T_e*cte.eV/cte.e

@task
def Psip_0(B_0, m_i, T_e, R, **kwargs):
    """ Scale of magnetic flux []"""
    return B_0*rho_s( B_0, m_i, T_e)*R

@task
def lambda_D( T_e, n_0, **kwargs):
    """Debye length [m]"""
    return np.sqrt(cte.epsilon_0*T_e*cte.eV/n_0/1e19/cte.e**2)

# independent numerical variables
@task
def mu ( m_i, **kwargs) : # proton mass, deuteron mass, triton mass
    """ electron mass ratio"""
    return -cte.m_e/ m_i

@task
def tau ( T_e, T_i, **kwargs) :
    """ ion to electron temperature ratio"""
    return T_i / T_e

@task
def beta( n_0, T_e, B_0, **kwargs) :  # n_0 in 1e-19, T_e in eV, B_0 in T
    """ electron plasma beta"""
    return n_0*1e19*T_e*cte.eV / (B_0**2 / cte.mu_0)

@task
def resistivity( n_0, T_e, B_0, **kwargs):
    """plasma resistivity"""
    return 0.51*np.sqrt( 2*cte.m_e)*cte.e**3*10/ 12/np.pi**(3/2)/\
            cte.epsilon_0**2 * n_0*1e19 / B_0 / (T_e*cte.eV)**(3/2)

@task
def R_0 ( B_0, m_i, T_e, R, **kwargs) :
    """machine radius relative to Larmor radius"""
    return R / rho_s(B_0, m_i, T_e)

@task
def epsilon_D ( n_0, B_0, m_i, **kwargs) :
    """Square Debye length relative to Larmor radius"""
    return cte.epsilon_0/n_0/1e19/m_i*B_0**2

# dependent numerical variables

@task
def viscosity_e( n_0, T_e, B_0, **kwargs) :
    """ parallel electron viscosity"""
    return 0.37/resistivity(n_0, T_e, B_0)

@task
def viscosity_i( n_0, T_e, B_0, m_i, **kwargs) :
    """ parallel ion viscosity"""
    #return 0.96/np.sqrt( m_i)/cte.e**3/10* 12*np.pi**(3/2)*\
    #        cte.epsilon_0**2 * B_0 * (T_i*cte.eV)**(3/2) /(n_0*1e19)
    return 0.69/resistivity(n_0, T_e, B_0)*np.sqrt( np.abs(mu(m_i)))

# Miscellaneous ( inverseaspectratio, q)
@task
def inverseaspectratio( a, R, **kwargs) :
    return a/R

@task
def aspectratio( a, R, **kwargs) :
    return R/a
@task
def a_0( B_0, m_i, T_e, a, **kwargs) :
    return R_0(B_0, m_i, T_e, a)

@task
def lparallel( q, R, **kwags):
    return q*R
@task
def lparallel_rhos( q, B_0, m_i, T_e, R, **kwargs):
    return q*R_0(B_0, m_i, T_e, R)

#def alpha_MHD( q, inverseaspectratio, beta, DeltaN):
#    """From Eich2020 paper Eq. (3)"""
#    return q**2/inverseaspectratio*beta*DeltaN

@task
def Myra (n_0, T_e, B_0, m_i, q, R, **kwargs):
    return resistivity(n_0, T_e, B_0) * lparallel_rhos(q, B_0, m_i, T_e, R)

@task
def Knudsen (m_i, n_0, T_e, B_0, q, R, **kwargs):
    """ 0.51 sqrt(mu)/Myra"""
    return 0.51*np.sqrt( abs( mu(m_i)))/ Myra (n_0, T_e, B_0, m_i, q, R)

@task
def pref(m_i, n_0, T_e, B_0, R, **kwargs):
    """ presssure reference p = pref * Psip """
    return 1/R_0(B_0, m_i, T_e, R)/beta( n_0, T_e, B_0)

@task
def nu_perp(n_0, T_e, B_0, m_i, **kwargs) :
    return 0.96/viscosity_i(n_0, T_e, B_0, m_i)

@task
def GAM_inv (B_0, m_i,T_e, R, **kwargs) :
    """ inverse GAM frequency """
    return 1.0/R_0(B_0, m_i, T_e, R)

# other scales
@task
def paraview_scale(B_0, m_i, T_e, **kwargs):
    return 1.0/rho_s(B_0, m_i, T_e)

@task
def GBS_timescale(B_0, m_i, T_e, R, **kwargs) :
    """ timescale in GBS units """
    return R_0(B_0, m_i, T_e, R)*omega_0_inv(B_0,m_i)

@task
def GBS_resistivity(n_0, T_e, B_0, m_i, R, **kwargs):
    """ resistivity in GBS units """
    return resistivity(n_0, T_e, B_0)*R_0(B_0, m_i, T_e, R)

@task
def CFL_diff( n_0, B_0, m_i, T_e, R, a, scaleR, Nz, **kwargs):
    """ CFL condition due to parallel diffusion """
    return (2.*np.pi/Nz*R_0(B_0, m_i, T_e, R)*(
        1-scaleR*a/R))**2/viscosity_e( n_0, T_e, B_0)


proton_mass   = cte.physical_constants["proton mass"][0]
deuteron_mass = cte.physical_constants["deuteron mass"][0]
triton_mass   = cte.physical_constants["triton mass"][0]

def numerical2physical( numerical, physical):
    """Invert numerical parameters to physical parameters

    Execution path depends on which keys are present in the two dictionaries
    numerical (dict) : "epsilon_D", "resistivity", "mu", "tau", "R_0", "beta"
    physical  (dict) : "R", "B_0", "T_e", "T_i", "n_0", "m_i"

    """

    if "mu" in numerical :
        physical["m_i"] = -cte.m_e/numerical["mu"]
    else :
        numerical["mu"] = -cte.m_e/physical["m_i"]

    if ( ("epsilon_D" in numerical) and (numerical["epsilon_D"] != 0) ) :
        if (not("beta" in numerical) or (numerical["beta"] == 0)):
            raise ValueError("beta must be present if epsilon_D is")
        def to_invert0( x,
                       for_e_D,for_mu, for_tau, for_beta, for_eta, for_R_0) :
            T_e, T_i, n_0, B_0, R = x
            m_i = -cte.m_e/for_mu
            return (tau(T_e, T_i) - for_tau,
                (beta( n_0, T_e, B_0) - for_beta)/for_beta,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (epsilon_D( n_0, B_0, m_i) - for_e_D)/for_e_D,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        print( "Invert for given numerical parameters")
        physical["T_e"], physical["T_i"], physical["n_0"], physical["B_0"], physical["R"]\
        = opt.fsolve( to_invert0, [1,1,1,1,1],args=(
            numerical["epsilon_D"],numerical["mu"],numerical["tau"],
            numerical["beta"],numerical["resistivity"],numerical["R_0"]))

    elif (("beta" in numerical) and (numerical["beta"] != 0) ):
        print( "Invert for given R")
        def to_invert1( x, R,
              for_mu, for_tau, for_beta, for_eta, for_R_0) :
            T_e, T_i, n_0, B_0 = x
            m_i = -cte.m_e/for_mu
            return (tau(T_e, T_i) - for_tau,
                (beta( n_0, T_e, B_0) - for_beta)/for_beta,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        physical["T_e"], physical["T_i"], physical["n_0"], physical["B_0"]\
        = opt.fsolve( to_invert1, [1,1,1,1],args=(
            physical["R"],numerical["mu"],numerical["tau"],
            numerical["beta"],numerical["resistivity"],numerical["R_0"]))
        numerical["epsilon_D"] = epsilon_D(physical["n_0"],
                                        physical["B_0"],
                                        physical["m_i"])
    else:
        print( "Invert for given R and B_0")
        if ((not "R" in physical) or (physical["R"] == 0) or
            (not "B_0" in physical) or (physical["B_0"] == 0)):
            raise ValueError ( "B_0 and R must be present in physical and be different from zero")
        def to_invert2( x, B_0, R,
                  for_mu, for_tau, for_eta, for_R_0) :
            T_e, T_i, n_0 = x
            m_i = -cte.m_e/for_mu
            return (tau(T_e, T_i) - for_tau,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        physical["T_e"], physical["T_i"], physical["n_0"]\
        = opt.fsolve( to_invert2, [1,1,1],args=(
            physical["B_0"],physical["R"],numerical["mu"],
            numerical["tau"],numerical["resistivity"],numerical["R_0"]))
        numerical["beta"] = beta( physical["n_0"],
                                    physical["T_e"],
                                    physical["B_0"])
        numerical["epsilon_D"] = epsilon_D(physical["n_0"],
                                        physical["B_0"],
                                        physical["m_i"])

def quantities() :
    """ A List of all available quantities """
    return tasks.keys()

def parameters2quantity( parameters, quantity) :
    if quantity in parameters.keys() :
        return  parameters[quantity]
    else :
        return tasks[ quantity]( **parameters)


def parameters2quantities( parameters, quantities) :
    """ Return a list of values corresponding to quantities """
    myList = []
    for quantity in quantities:
        myList.append( parameters2quantity( parameters, quantity))
    return myList

def load_calibration_default():
    """ generate default feltor input parameters for calibration"""
    inputfile= {
        "grid":
        {
            "n"  : 3,
            "Nx" : 48,
            "Ny" : 96,
            "Nz" : 8,
            "scaleR" : [1.45,1.25],
            "scaleZ" : [2.6, 2.25]
        },
        "advection":
        {
            "slope-limiter" : "none"
        },
        "timestepper":
        {
            "tableau" : "Bogacki-Shampine-4-2-3",
            "type" : "adaptive",
            "rtol": 1e-5,
            "atol" : 1e-6,
            "output-mode" : "deltaT",
            "deltaT" : 100,
            "reject-limit" : 10
        },
        "regularization":
        {
            "order" : 2,
            "direction": "forward",
            "nu_perp_n" : 1e-5,
            "nu_perp_u" : 1e-5,
            "nu_parallel_n" : 5e2
        },
        "elliptic":
        {
            "stages" : 2,
            "eps_pol"    : [1e-6, 0.5],
            "eps_gamma" : 1e-7,
            "eps_ampere" : 1e-7,
            "direction" : "forward",
            "jumpfactor" : 1.0
        },
        "FCI":
        {
            "refine" : [5,5],
            "rk4eps" : 1e-6,
            "periodify": True,
            "bc" : "along_field",
            "interpolation-method" : "linear"
        },
        "physical":
        {
            "mu"          : -0.000272121,
            "epsilon_D" : 4.1458919332419e-05,
            "tau"         : 0.5,
            "beta"        : 1e-5,
            "resistivity" : 1e-4,
            "viscosity" : "braginskii"
        },
        "output":
        {
            "type" : "netcdf",
            "itstp"  : 2,
            "maxout" : 1,
            "compression" : [1,1],
            "equations":[
                "Basic",
                "Mass-conserv",
                "Energy-theorem",
                "Toroidal-momentum",
                "Parallel-momentum",
                "Zonal-Flow-Energy",
                "COCE"
            ]
        },
        "source":
        {
            "minne" : 0.2,
            "minrate" : 1.0,
            "minalpha" : 0.1,
            "type" : "influx",
            "rate" : 1e-4,
            "ntilde":
            {
                "type" : "zero",
            },
            "profile":
            {
                "type" : "aligned",
                "npeak" : 1.0,
                "nsep" : 0.0,
                "background" : 0.0
            },
            "damping":
            {
                "type": "alignedX",
                "alpha" : 0.2,
                "boundary" : 0.55
            }
        },
        "init":
        {
            "type" : "fields",
            "density" :
            {
                "type" : "ne",
                "ntilde":
                {
                    "type" : "turbulence",
                    "amplitude"  : 1e-4,
                    "revolutions" : 1,
                    "parallel" : "gaussian",
                    "sigma_z" : 0.5
                },
                "profile":
                {
                    "type" : "aligned",
                    "npeak" : 4.0,
                    "nsep" : 1.0,
                    "background" : 0.2
                },
                "damping":
                {
                    "type": "alignedPFR",
                    "alpha" : [0.1,0.04],
                    "boundary" : [1.14,0.96]
                }
            },
            "potential" : {"type" : "zero_pol"},
            "velocity" :
            {
                "type" : "ui",
                "profile" : "linear_cs"
            },
            "aparallel" : {"type" : "zero"}
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
                "penalize-rhs" : False,
                "nwall" : 0.2,
                "uwall" : 0.0
            },
            "sheath":
            {
                "type": "bohm",
                "boundary": 0.09375, # 3/32
                "alpha": 0.0625, # 2/32
                "penalization" : 1e+0,
                "penalize-rhs" : False,
                "coordinate" : "s",
                "max_angle" : 4
            },
            "bc" :
            {
                "density" : ["NEU", "NEU"],
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

def load_default_config ():
    configfile = {
        # for feltordiag
        "n": 3,
        "Npsi": 64,
        "Neta": 640,
        "Kphi": 10,
        "fx_0" : 0.125,
        "fsa" : "toroidal-average", # or "convoluted-toroidal-average"
        "diagnostics":
        [
            "fsa",
            "fsa2d",
            "cta2d",
            "cta2dX",
            "fluc2d",
            "ifs",
            "std_fsa",
            "ifs_lcfs",
            "ifs_norm"
        ],
        # for interpolate in 3d
        "fine-grid-factor" : 2,
        "time-reduction-factor" : 10
    }
    return configfile
