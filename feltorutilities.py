import numpy as np
import scipy.constants as cte
import scipy.optimize as opt

tasks = {}
task = lambda f: tasks.setdefault( f.__name__, f)

#WE PASS **kwargs TO ALL FUNCTIONS SO THAT WE CAN CALL ALL WITH DICTIONARIES

# identities
@task
def T_e( T_e, **kwargs) :
    return T_e
@task
def T_i( T_i, **kwargs) :
    return T_i
@task
def m_i( m_i, **kwargs) :
    return m_i
@task
def R( R, **kwargs) :
    return R
@task
def n_0( n_0, **kwargs) :
    return n_0
@task
def B_0( B_0, **kwargs) :
    return B_0


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
def tau_i ( T_e, T_i, **kwargs) :
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
def minor_radius( inverseaspectratio, R, **kwargs) :
    """as a function of inverse aspect ratio and R"""
    return inverseaspectratio*R
@task
def minor_radius_rhos( inverseaspectratio, B_0, m_i, T_e, R, **kwargs) :
    """as a function of inverse aspect ratio and R"""
    return inverseaspectratio*R_0(B_0, m_i, T_e, R)

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
    return 0.51*sqrt( abs( mu(m_i)))/Myra( resistivity(n_0, T_e, B_0),
            lparallel_rhos(q, B_0, m_i, T_e, R))

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
def CFL_diff( n_0, B_0, m_i, T_e, R, inverseaspectratio, scaleR, Nz, **kwargs):
    """ CFL condition due to parallel diffusion """
    return (2.*np.pi/Nz*R_0(B_0, m_i, T_e, R)*(
        1-scale_R*inverseaspectratio))**2/viscosity_e( n_0, T_e, B_0)


proton_mass   = cte.physical_constants["proton mass"][0]
deuteron_mass = cte.physical_constants["deuteron mass"][0]
triton_mass   = cte.physical_constants["triton mass"][0]

def physical2numerical( physical, numerical):
    """Compute numerical parameters from physical parameters"""
    numerical["mu"] = mu( physical["m_i"])
    numerical["tau_i"] = tau_i( physical["T_e"], physical["T_i"])
    numerical["beta"] = beta( physical["n_0"], physical["T_e"], physical["B_0"])
    numerical["resistivity"] = resistivity( physical["n_0"], physical["T_e"], physical["B_0"])
    numerical["R_0"] = R_0( physical["B_0"], physical["m_i"], physical["T_e"], physical["R"])
    numerical["epsilon_D"] = epsilon_D( physical["n_0"], physical["B_0"], physical["m_i"])

def numerical2physical( numerical, physical):
    """Invert numerical parameters to physical parameters

    Execution path depends on which keys are present in the two dictionaries

    """

    physical["m_i"] = -cte.m_e/numerical["mu"]

    if ( ("epsilon_D" in numerical) and (numerical["epsilon_D"] != 0) ) :
        if (not("beta" in numerical) or (numerical["beta"] == 0)):
            raise ValueError("beta must be present if epsilon_D is")
        def to_invert0( x,
                       for_e_D,for_mu, for_tau_i, for_beta, for_eta, for_R_0) :
            T_e, T_i, n_0, B_0, R = x
            m_i = -cte.m_e/for_mu
            return (tau_i(T_e, T_i) - for_tau_i,
                (beta( n_0, T_e, B_0) - for_beta)/for_beta,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (epsilon_D( n_0, B_0, m_i) - for_e_D)/for_e_D,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        print( "Invert for given numerical parameters")
        physical["T_e"], physical["T_i"], physical["n_0"], physical["B_0"], physical["R"]\
        = opt.fsolve( to_invert0, [1,1,1,1,1],args=(
            numerical["epsilon_D"],numerical["mu"],numerical["tau_i"],
            numerical["beta"],numerical["resistivity"],numerical["R_0"]))

    elif (("beta" in numerical) and (numerical["beta"] != 0) ):
        print( "Invert for given R_0")
        def to_invert1( x, R,
              for_mu, for_tau_i, for_beta, for_eta, for_R_0) :
            T_e, T_i, n_0, B_0 = x
            m_i = -cte.m_e/for_mu
            return (tau_i(T_e, T_i) - for_tau_i,
                (beta( n_0, T_e, B_0) - for_beta)/for_beta,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        physical["T_e"], physical["T_i"], physical["n_0"], physical["B_0"]\
        = opt.fsolve( to_invert1, [1,1,1,1],args=(
            physical["R"],numerical["mu"],numerical["tau_i"],
            numerical["beta"],numerical["resistivity"],numerical["R_0"]))
        numerical["epsilon_D"] = epsilon_D(physical["n_0"],
                                        physical["B_0"],
                                        physical["m_i"])
    else:
        print( "Invert for given R and B_0")
        if ((not "R" in physical) or (physical["R"] == 0) or
            (not "B_0" in physical) or (physical["B_0"] == 0)):
            raise ValueError ( "B_0 and R_0 must be present in physical and be different from zero")
        def to_invert2( x, B_0, R,
                  for_mu, for_tau_i, for_eta, for_R_0) :
            T_e, T_i, n_0 = x
            m_i = -cte.m_e/for_mu
            return (tau_i(T_e, T_i) - for_tau_i,
                (resistivity(n_0, T_e, B_0) - for_eta)/for_eta,
                (R_0(B_0, m_i,T_e, R) - for_R_0)/for_R_0
               )
        physical["T_e"], physical["T_i"], physical["n_0"]\
        = opt.fsolve( to_invert2, [1,1,1],args=(
            physical["B_0"],physical["R"],numerical["mu"],
            numerical["tau_i"],numerical["resistivity"],numerical["R_0"]))
        numerical["beta"] = beta( physical["n_0"],
                                    physical["T_e"],
                                    physical["B_0"])
        numerical["epsilon_D"] = epsilon_D(physical["n_0"],
                                        physical["B_0"],
                                        physical["m_i"])

def parameters2quantities( parameters, quantities) :
    """ Return a list of values corresponding to quantities """
    myList = []
    for quantity in quantities:
        myList.append(tasks[ quantity]( **parameters))
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