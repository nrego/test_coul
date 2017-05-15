## Do vdw calculation from md trajectory##

from __future__ import division, print_function

import MDAnalysis
import numpy as np

import argparse
try:
    u_for_prev = u_for
except:
    pass

from mdtools import dr
import itertools
from IPython import embed


parser = argparse.ArgumentParser('calculate coulombic interactions for methanol')
parser.add_argument('--lmbda', required=True, type=float,
                    help='This lambda value')
parser.add_argument('--for-lmbda', required=True, type=str,
                    help='comma separated array of foreign lambdas')
parser.add_argument('-s', '--tprfile', required=True, type=str,
                    help='tpr file')
parser.add_argument('-f', '--trajfile', required=True, type=str,
                    help='traj file (trr or xtc)')
parser.add_argument('-o', '--outfile', default='my_diffs.dat',
                    help='output file for diffs')
parser.add_argument('--n-periodic-images', type=int, default=0,
                    help='Number of periodic multiples of each box vector to consider when' \
                    'calculating coulombic interactions; e.g if 1 will consider all 26 additional' \
                    'adjacent box images (by shifting each component of the box vector -1, 0, or 1).' \
                    'default is 0; i.e. do not consider any periodic images')


args = parser.parse_args()

univ = MDAnalysis.Universe(args.tprfile, args.trajfile)

# box vector
box = univ.dimensions[:3] / 10.0

#shifts = [0]
n_shifts = args.n_periodic_images
#assert n_shifts >= 0, "--n-periodic-shifts arg must be an integer >= 0"
shifts = []
print("generating shift vectors for n_shift: {}!".format(n_shifts))
for n_shift in range(n_shifts+1):
    shifts.append(-n_shift)
    shifts.append(n_shift)
shifts = np.unique(shifts)
shift_vectors = np.array([p for p in itertools.product(shifts, repeat=3)])
print("  ...{} shift vectors generated!".format(shift_vectors.shape[0]))
# Sanity
#embed()

n_images = shift_vectors.shape[0]

outfile = args.outfile

lmbda = args.lmbda

for_lmbda_str = args.for_lmbda
for_lmbdas = sorted([float(for_lmbda) for for_lmbda in for_lmbda_str.split(',')])

n_for_lmbdas = len(for_lmbdas)

atm_indices = range(univ.atoms.n_atoms)
alc_indices = range(6)

atmtypes = [('C3', 'C3'),
            ('OH', 'OH'),
            ('H1', 'H1'),
            ('H1', 'H1'),
            ('H1', 'H1'),
            ('HO', 'HO')]


charges = [(0.12010, 0.0),
           (-0.60030, 0.0),
           (0.02770, 0.0),
           (0.02770, 0.0),
           (0.02770, 0.0),
           (0.39710, 0.0)]

excls = {0: (0,1,2,3,4,5),
         1: (0,1,2,3,4,5),
         2: (0,1,2,3,4,5),
         3: (0,1,2,3,4,5),
         4: (0,1,2,3,4,5),
         5: (0,1,2,3,4,5)}

# 1-4 pairs
pairs = {0: (),
         1: (),
         2: [5],
         3: [5],
         4: [5],
         5: [2,3,4]}


ke = 138.9354859 # conv factor to get kJ/mol

fudge = 0.83333333

n_frames = univ.trajectory.n_frames
n_frames = 10
my_diffs = np.zeros((n_frames, n_for_lmbdas+1))

for i_frame in range(n_frames):
    if i_frame % 10 == 0:
        print("doing frame {} of {}".format(i_frame, n_frames))
    univ.trajectory[i_frame]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    my_diffs[i_frame, 0] = univ.trajectory.time

    for for_lmbda_idx, for_lmbda in enumerate(for_lmbdas):
        this_coul = 0.0 # at this lambda
        for_coul = 0.0 # at foreign lambda
        
        for i in alc_indices:
            # Everything but this atom's exclusions
            incl_indices = np.setdiff1d(atm_indices, excls[i])
            pair_indices = pairs[i]
            atm_i = univ.atoms[i]

            i_type_a, i_type_b = atmtypes[i]
            assert atm_i.type == i_type_a

            i_charge_a, i_charge_b = charges[i]
            np.testing.assert_almost_equal(i_charge_a, atm_i.charge)

            # any regular (not 14) included atoms
            # optionally include periodic images
            
            for j in incl_indices:
                if j <= i:
                    continue
                if j in pair_indices:
                    continue

                atm_j = univ.atoms[j]

                #j_type_a, j_type_b = atmtypes[j]
                #assert atm_j.type == j_type_a

                #j_charge_a, j_charge_b = charges[j]
                #np.testing.assert_almost_equal(atm_j.charge, j_charge_a)
                j_charge_a = j_charge_b = atm_j.charge
                for this_shift in shift_vectors:
                    j_pos = atm_j.position + this_shift*box
                    r = np.sqrt(np.sum((atm_i.position - j_pos)**2))

                    a_coul = ke*((i_charge_a * j_charge_a)/ r)
                    b_coul = ke*((i_charge_b * j_charge_b)/ r)

                    this_coul += (1-lmbda)*a_coul + (lmbda)*b_coul
                    for_coul += (1-for_lmbda)*a_coul + (for_lmbda)*b_coul

            for j in pair_indices:
                if j <= i:
                    continue

                atm_j = univ.atoms[j]

                j_type_a, j_type_b = atmtypes[j]
                assert atm_j.type == j_type_a

                j_charge_a, j_charge_b = charges[j]
                np.testing.assert_almost_equal(atm_j.charge, j_charge_a)

                r = np.sqrt(np.sum((atm_i.position - atm_j.position)**2))

                a_coul = ke*fudge*((i_charge_a * j_charge_a)/ r)
                b_coul = ke*fudge*((i_charge_b * j_charge_b)/ r)

                this_coul += (1-lmbda)*a_coul + (lmbda)*b_coul
                for_coul += (1-for_lmbda)*a_coul + (for_lmbda)*b_coul

        #embed()
        my_diffs[i_frame, for_lmbda_idx+1] = for_coul - this_coul

headerstr = 'time (ps)'
headerstr += '    '.join(str(for_lmbda) for for_lmbda in for_lmbdas)
np.savetxt(outfile, my_diffs, header=headerstr, fmt='%2.4f  %2.8f')    


