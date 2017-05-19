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
univ = MDAnalysis.Universe('run.tpr', 'traj_pbc.trr')
from mdtools import dr

from math import erf, erfc


atm_indices = range(2)
alc_indices = range(2)

atmtypes = [('DUM_HC', 'DUM_HC'),
            ('DUM_OW_spc', 'DUM_OW_spc')]


charges = [0.1,-0.1]

rtol = 1e-6
rc = 0.9

# determine ewald stuff
xvals = np.linspace(0,5,1000000)
vals = np.zeros_like(xvals)
for i, x in enumerate(xvals):
    vals[i] = erfc(x)

etol_ind = np.argmax(vals < rtol)
beta_smooth = xvals[etol_ind] / rc
print("smoothing width: {}".format(beta_smooth))

ke = 138.9354859 # conv factor to get kJ/mol

fudge = 1.0

ewald_self = -beta_smooth*ke/np.sqrt(np.pi) * 2 * (0.01)
print("self-self Ewald correction: {} (kJ/mol)".format(ewald_self))

box = univ.dimensions[:3] / 10.0

n_shifts = 0
shifts = []
print("generating shift vectors for n_shift: {}!".format(n_shifts))
for n_shift in range(n_shifts+1):
    shifts.append(-n_shift)
    shifts.append(n_shift)
shifts = np.unique(shifts)
shift_vectors = np.array([p for p in itertools.product(shifts, repeat=3)])
print("  ...{} shift vectors generated!".format(shift_vectors.shape[0]))

n_waves = 5
waves = []
print("generating wave vectors for n_waves: {}!".format(n_waves))
for n_wave in range(1,n_waves+1):
    waves.append(-n_wave)
    waves.append(n_wave)
waves = np.unique(waves)
wave_vectors = np.array([p for p in itertools.product(waves, repeat=3)])
print("  ...{} wave vectors generated!".format(wave_vectors.shape[0]))

n_frames = univ.trajectory.n_frames
n_frames = 1
my_diffs = np.zeros((n_frames, 6))

for i_frame in range(n_frames):
    if i_frame % 100 == 0:
        print("frame {} of {}".format(i_frame, n_frames))
    #univ.trajectory[i_frame]
    univ.trajectory[50]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    coul = 0.0
    ewald_sr = 0.0
    ewald_lr = 0.0
    ewald_lr2 = 0.0

    atm_i = univ.atoms[0]
    atm_j = univ.atoms[1]

    pos_i = atm_i.position
    pos_j = atm_j.position

    my_diffs[i_frame, 0] = univ.trajectory.time

    for this_shift in shift_vectors:
        shift = this_shift * box
        r = np.linalg.norm(pos_i - pos_j+shift)
        coul += ke*(-0.01)*(1/r)
        if not np.array_equal(shift, np.zeros(3)):
            coul += ke*(0.01)*(1/np.linalg.norm(shift))
        else:
            ewald_sr += ke*(-0.01)*(1/r)*erfc(beta_smooth*r)

    for this_wave in wave_vectors:
        k = this_wave * 2 * np.pi / box

        s_k = (0.1)*np.exp(1j*np.dot(k,pos_i)) - (0.1)*np.exp(1j*np.dot(k,pos_j))
        k_sq = np.dot(k,k)
        r_ij = pos_i - pos_j
        exp1 = np.exp(1j * np.dot(k, r_ij))
        exp2 = np.exp(-k_sq/(2*beta_smooth)**2)
        ewald_lr += (4*np.pi*ke / box.prod())*(-0.01)*(1/k_sq)*exp1*exp2
        ewald_lr += (4*np.pi*ke / box.prod())*(0.01)*(1/k_sq)*exp2
        ewald_lr2 += (2*np.pi*ke / box.prod()) * (1/k_sq) * np.abs(s_k)**2 * exp2

    r = np.linalg.norm(pos_i - pos_j)
    my_diffs[i_frame, 1] = r
    my_diffs[i_frame, 2] = ke * (-0.01) * (1/r)
    #beyond_cut = my_diffs[:,1] > rc
    #my_diffs[beyond_cut, 2] = 0.0
    my_diffs[i_frame, 3] = coul
    my_diffs[i_frame, 4] = ewald_sr
    my_diffs[i_frame, 5] = ewald_lr2 #+ ewald_self

header = 'time(ps)        r(nm)        Coul_simple(kJ/mol)        Coul_with_periodic(kJ/mol)         U_dir        U_rec'
np.savetxt('my_diffs.dat', my_diffs)

#plt.plot(my_diffs[:,0], -my_diffs[:,5], label='calc')
#plt.plot(my_diffs[:,0], -gmx_wave, label='gmx')
#plt.legend()