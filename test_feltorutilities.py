import feltorutilities as fp
import numpy as np


R_0 = fp.R_0(0.07, fp.deuteron_mass, 10, 0.6)
omega_0_inv = fp.omega_0_inv( 0.07, fp.deuteron_mass)
resistivity = fp.resistivity( 0.1, 10, 0.07)
print(f"Hello mu {fp.mu(fp.proton_mass):.2e}, {fp.mu(fp.deuteron_mass):.2e},\
        {fp.mu(fp.triton_mass):.2e}")
print(f"Hello beta\t\t{fp.beta( 0.1, 10, 0.07):.2e}")
print(f"Hello resistivity\t{resistivity:.2e}")
print(f"Hello viscosity_e\t{fp.viscosity_e( 0.1, 10, 0.07):.2e}")
print(f"Hello viscosity_i\t{fp.viscosity_i( 0.1, 10, 0.07, fp.deuteron_mass):.2e}")
print(f"Hello R_0\t\t{R_0:.2f}")
print(f"Hello omega^(-1)\t{omega_0_inv:.2e}")

physical = {"R" : 0.6}
numerical = {"mu" : -2.72e-4,"tau":1,"beta":4.11e-4,"resistivity" : 3.81e-5,"R_0" : 91.94}
fp.numerical2physical( numerical, physical)

print( physical)
print( numerical)
numerical["epsilon_D"] = 0
numerical["beta"] = 0
fp.numerical2physical( numerical, physical)

print( physical, numerical)
fp.numerical2physical( numerical, physical)
print( physical, numerical)

table_list = fp.quantities()
physical["a"] = 0.15
physical["q"] = 2
physical["scaleR"] = 1.2
physical["Nz"] = 32
print( table_list)
print( fp.parameters2quantities( physical, table_list ))
