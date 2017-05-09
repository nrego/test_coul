## Do vdw calculation from md trajectory##

from __future__ import division, print_function

import MDAnalysis
import numpy as np

import argparse
try:
    u_for_prev = u_for
except:
    pass
univ = MDAnalysis.Universe('run.tpr', 'traj.trr')
from mdtools import dr



atm_indices = range(2)
alc_indices = range(2)

atmtypes = [('DUM_HC', 'DUM_HC'),
            ('DUM_OW_spc', 'DUM_OW_spc')]


charges = [(0.0, 0.5),
           (0.0, -0.5)]


ke = 138.9354859 # conv factor to get kJ/mol

fudge = 1.0



n_frames = univ.trajectory.n_frames
my_diffs = np.zeros((n_frames, 2))

for i_frame in range(n_frames):
    univ.trajectory[i_frame]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    coul_a = 0.0
    coul_b = 0.0
    
    my_diffs[i_frame, 0] = univ.trajectory.time

    for i in alc_indices:
        atm_i = univ.atoms[i]

        i_type_a, i_type_b = atmtypes[i]
        assert atm_i.type == i_type_a

        i_charge_a, i_charge_b = charges[i]
        assert atm_i.charge == i_charge_a

        for j in alc_indices:
            if j <= i:
                continue

            atm_j = univ.atoms[j]

            j_type_a, j_type_b = atmtypes[j]
            assert atm_j.type == j_type_a

            j_charge_a, j_charge_b = charges[j]
            assert atm_j.charge == i_charge_a

            r_sq = np.sum((atm_i.position - atm_j.position)**2)

            coul_a += ke*fudge*((i_charge_a*j_charge_b)/(r_sq))
            coul_b += ke*fudge*((i_charge_b*j_charge_b)/(r_sq))


    my_diffs[i_frame, 1] = coul_b - coul_a

