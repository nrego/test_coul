## Do vdw calculation from md trajectory##

from __future__ import division, print_function

import MDAnalysis
import numpy as np
import itertools
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


charges = [0.1,-0.1]


ke = 138.9354859 # conv factor to get kJ/mol

fudge = 1.0

box = univ.dimensions[:3] / 10.0

n_shifts = 1
shifts = []
print("generating shift vectors for n_shift: {}!".format(n_shifts))
for n_shift in range(n_shifts+1):
    shifts.append(-n_shift)
    shifts.append(n_shift)
shifts = np.unique(shifts)
shift_vectors = np.array([p for p in itertools.product(shifts, repeat=3)])
print("  ...{} shift vectors generated!".format(shift_vectors.shape[0]))

n_frames = univ.trajectory.n_frames
my_diffs = np.zeros((n_frames, 2))

for i_frame in range(n_frames):
    if i_frame % 100 == 0:
        print("frame {} of {}".format(i_frame, n_frames))
    univ.trajectory[i_frame]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    coul = 0.0

    atm_i = univ.atoms[0]
    atm_j = univ.atoms[1]

    pos_i = atm_i.position
    pos_j = atm_j.position

    my_diffs[i_frame, 0] = univ.trajectory.time

    for this_shift in shift_vectors:
        shift = this_shift * box
        r = np.linalg.norm(pos_i - pos_j+shift)
        coul += ke*(0.01)*(1/r)
        if not np.array_equal(shift, np.zeros(3)):
            coul += 2*ke*(-0.01)*(1/np.linalg.norm(shift))

    my_diffs[i_frame, 1] = coul

