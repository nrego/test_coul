# Analytically estimate free energy of removing charge in vacuum
from __future__ import division

import numpy as np



#H3            H6                 x|
#  \          /                    |
#   \        /                     |
#    C1--- O2                      |-------z
#   / \                           / 
#  /   \                         /
# H4   H5                      y/


# Charges
fudge = 0.83333333 # All H-H coulombic interactions are 1-4's
ke = 138.9354859  # coul conversion factor to get kJ/mol from (e**2)/(nm)

T = 300
beta = 1 / (8.3144598e-3 * T)

# Dihedral torsional energy from angle (in degrees), in kJ/mol
def dihed_energy(phi, k_phi=0.69733334, mult=3, phase=0):
    phi_rad = phi*np.pi/180.0

    return k_phi*(1 + np.cos(mult*phi_rad+phase) )

## constants (equil bond lengths, angles - treated as constraints)
# in nm or radians, resp
r_co = 0.1426
r_ch = 0.1093
r_oh = 0.0974

thet_coh = theta2 = 108.16004647*np.pi / (180.)
thet_och = theta1 = 109.88004686*np.pi / (180.)
thet_hch = theta3 = 109.55004724*np.pi / (180.)

# c1
r1 = np.array([0.,0.,0.])
# o2
r2 = np.array([0., 0., r_co])
# h6
r6 = np.array([r_oh*np.sin(theta2), 0., -r_oh*np.cos(theta2)+r_co])

# normal to plane C1-O2-H6
# (r6-r2) x (r1-r2) / |n2|
n2 = np.array([0., 1., 0.])

# Distance between two methyl hydrogens - fixed because angle H-C-H and dist r_ch are fixed
r_hh = np.sqrt(2*r_ch**2 * (1-np.cos(theta3)))

# z coord for r3,r4,r5 (methyl H's). All are fixed by theta 1 and C-H
r_z = r_ch * np.cos(theta1)

def get_r(phi):
    phi_rad = (phi * np.pi) / 180.

    r3_x = np.sqrt(r_ch**2 - r_z**2) * np.cos(phi_rad)

    r3_y = np.sqrt(r_ch**2 - r3_x**2 - r_z**2)
    if phi < 180:
        r3_y = -r3_y

    return np.array([r3_x, r3_y, r_z])


phivals = np.linspace(0, 360, 1000)

# first col is at lambda=0 (coul on), second at lambda=1 (coul off)
energies = np.zeros((phivals.size, 2))

h1_charge = 0.02770
ho_charge = 0.39710

for i, phi_val in enumerate(phivals):

    phi2 = (phi_val + 120) % 360
    phi3 = (phi2 + 120) % 360

    r3 = get_r(phi_val)
    r4 = get_r(phi2)
    r5 = get_r(phi3)

    r36 = np.linalg.norm(r3-r6)
    r46 = np.linalg.norm(r4-r6)
    r56 = np.linalg.norm(r5-r6)

    e_coul = fudge*ke*(h1_charge*ho_charge) * ( (1/r36) + (1/r46) + (1/r56))
    e_dihed = dihed_energy(phi_val) + dihed_energy(phi2) + dihed_energy(phi3)


    energies[i, 0] = beta*(e_coul+e_dihed)
    energies[i, 1] = beta*(e_dihed)





