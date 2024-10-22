import numpy as np

import feltorutilities as fp

# run with pytest -s .

n_0 = 0.1  # 1e19
B_0 = 0.07  # T
m_i = fp.deuteron_mass
T_e = 10  # eV
R = 0.6


def test_elemental_functions():
    # Just print results
    print("TEST ELEMENTAL")
    print("Physical ", R, m_i, T_e, n_0, B_0)
    R_0 = fp.R_0(B_0, m_i, T_e, R)
    omega_0_inv = fp.omega_0_inv(B_0, m_i)
    collisionality = fp.collisionality(n_0, T_e, m_i, B_0)
    print(f" mue {fp.mue(fp.proton_mass):.2e}, {fp.mue(fp.deuteron_mass):.2e},\
            {fp.mue(fp.triton_mass):.2e}")
    print(f" beta\t\t{fp.beta( 0.1, 10, 0.07):.2e}")
    print(f" collisionality\t{collisionality:.2e}")
    print(f" viscosity_e\t{fp.viscosity_e( 0.1, 10, 0.07):.2e}")
    print(f" viscosity_i\t{fp.viscosity_i( 0.1, 10, 0.07, fp.deuteron_mass):.2e}")
    print(f" R_0\t\t{R_0:.2f}")
    print(f" omega^(-1)\t{omega_0_inv:.2e}")


def test_thermal_conversion():
    print("TEST THERMAL CONVERSION")
    physical = {"R": 0.6}
    numerical = {
        "mue": 2.72e-4,
        "beta": 4.11e-4,
        "collisionality": 4.52e-3,
        "R_0": 91.94,
    }
    fp.numerical2physical(numerical, physical, verbose=True)
    assert np.abs(physical["R"] - R) / R < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - T_e) / T_e < 1e-2
    assert np.abs(physical["n_0"] - n_0) / n_0 < 1e-2
    assert np.abs(physical["B_0"] - B_0) / B_0 < 1e-2
    # print( physical)
    # print( "Should be: ", R, m_i, T_e, n_0, B_0)

    numerical["epsilon_D"] = 0
    numerical["beta"] = 0
    fp.numerical2physical(numerical, physical, verbose=True)
    assert np.abs(physical["R"] - R) / R < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - T_e) / T_e < 1e-2
    assert np.abs(physical["n_0"] - n_0) / n_0 < 1e-2
    assert np.abs(physical["B_0"] - B_0) / B_0 < 1e-2
    # print( physical)
    # print( "Should be: ", R, m_i, T_e, n_0, B_0)

    print()
    fp.numerical2physical(numerical, physical, verbose=True)
    # print( physical)
    # print( "Should be: ", R, m_i, T_e, n_0, B_0)
    assert np.abs(physical["R"] - R) / R < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - T_e) / T_e < 1e-2
    assert np.abs(physical["n_0"] - n_0) / n_0 < 1e-2
    assert np.abs(physical["B_0"] - B_0) / B_0 < 1e-2


def test_feltor_conversion():
    print("TEST FELTOR CONVERSION")
    physical = {"R": 0.6}
    numerical = {
        "mu": -2.72e-4,
        "tau": 1,
        "beta": 4.11e-4,
        "resistivity": 3.81e-5,
        "R_0": 91.94,
    }
    fp.numerical2physical(numerical, physical, verbose=True)
    # print( physical)
    # print( "Should be: ", R, m_i, T_e, n_0, B_0)
    assert np.abs(physical["R"] - R) / R < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - T_e) / T_e < 1e-2
    assert np.abs(physical["n_0"] - n_0) / n_0 < 1e-2
    assert np.abs(physical["B_0"] - B_0) / B_0 < 1e-2
    assert np.abs(physical["T_i"] - T_e) / T_e < 1e-2

    print()
    numerical["epsilon_D"] = 0
    numerical["beta"] = 0
    fp.numerical2physical(numerical, physical, verbose=True)
    # print( physical)
    # print( "Should be: ", R, m_i, T_e, n_0, B_0)
    assert np.abs(physical["R"] - R) / R < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - T_e) / T_e < 1e-2
    assert np.abs(physical["n_0"] - n_0) / n_0 < 1e-2
    assert np.abs(physical["B_0"] - B_0) / B_0 < 1e-2
    assert np.abs(physical["T_i"] - T_e) / T_e < 1e-2

    print()
    fp.numerical2physical(numerical, physical, verbose=True)
    # print( physical)
    # print( "Should be: ", R, m_i, T_e, n_0, B_0)
    assert np.abs(physical["R"] - R) / R < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - T_e) / T_e < 1e-2
    assert np.abs(physical["n_0"] - n_0) / n_0 < 1e-2
    assert np.abs(physical["B_0"] - B_0) / B_0 < 1e-2
    assert np.abs(physical["T_i"] - T_e) / T_e < 1e-2


def test_tcv_conversion():
    print("TEST FELTOR CONVERSION 2")
    physical = {
        "beta": 9.637284991047368e-05,
        "epsilon_D": 0.00022903529747047053,
        "mu": -0.00027244371074816386,
        "nu_parallel": [10100.839709595633, 3359.5918043924403],
        "resistivity": 3.401796165299197e-06,
        "tau": 1.0,
        "viscosity": "value",
    }
    numerical = {"R_0": 900.4806396380083, **physical}
    fp.numerical2physical(numerical, physical, verbose=True)
    print( physical)
    assert np.abs(physical["R"] - 0.9) / 0.9 < 1e-2
    assert np.abs(physical["m_i"] - m_i) / m_i < 1e-2
    assert np.abs(physical["T_e"] - 41.4) / 41.4 < 1e-2
    assert np.abs(physical["n_0"] - 1.0) / 1.0 < 1e-2
    assert np.abs(physical["B_0"] - 0.93) / 0.93 < 1e-2
    assert np.abs(physical["T_i"] - 41.4) / 41.4 < 1e-2


def test_list_parameters():
    print()
    print("TETS LIST EVALUATION")
    physical = {
        "a": 0.15,
        "q": 2,
        "scaleR": 1.2,
        "Nz": 32,
        "n_0": 0.1,
        "B_0": 0.07,
        "m_i": fp.deuteron_mass,
        "T_e": 10,
        "T_i": 10,
        "R": 0.6,
    }
    table_list = fp.quantities()
    print(table_list)
    print(fp.parameters2quantities(physical, table_list))


def test_compass():
    # not sure what the test is here...

    print("COMPASS TEST")
    physicals = []
    for res in np.sort([3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6]):
        params = {
            "name": "COMPASS",
            "beta": 1e-4,
            "resistivity": res,  # change both to change n_0
            "tau": 1,
            "m_i": fp.deuteron_mass,
            "R_0": 545,
            "R": 0.545,
            "a": 0.175,
            "q": 2,
            "scaleR": 1.45,
            "Nz": 32,
        }
        physical = {
            "m_i": fp.deuteron_mass,
            "R": 0.545,
            "a": 0.175,
            "q": 2,
            "scaleR": 1.45,
            "Nz": 32,
        }
        fp.numerical2physical(params, physical)
        physicals.append(physical)

    # print(physicals[0])
    print()
    # numerical parameters that only differ by eps are the same
    test = []
    epsilon_D = fp.epsilon_D(**physicals[0])
    print("epsilon_D is ", epsilon_D)
    # for epsilon in [0.1*epsilon_D, epsilon_D]:
    params = {
        "name": "COMPASS",
        "beta": 1e-4,
        "resistivity": 1e-6,  # change both to change n_0
        "tau": 1,
        "R_0": 545,
        "epsilon_D": epsilon_D,
        "a": 0.175,
        "q": 2,
        "scaleR": 1.45,
        "Nz": 32,
    }
    physical = {"m_i": fp.deuteron_mass}
    fp.numerical2physical(params, physical, verbose=True)
    # print(physical)
    test.append(physical)
    print(fp.epsilon_D(**physical))
