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

for_lmbdas = [0.5,1.0]
n_for_lmbdas = len(for_lmbdas)
this_lmbda = 0.0

atm_indices = range(2)
alc_indices = range(2)

excls = {0: (),
         1: ()}

atmtypes = [('DUM_HC', 'DUM_HC'),
            ('DUM_OW_spc', 'DUM_OW_spc')]


charges = {0:(0.5, 0.0),
           1:(-0.5, 0.0)}

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

fudge = 0.8333333

ewald_self_this = 0.0
ewald_self_for = [0.0 for for_lmbda in for_lmbdas]

for alc_idx in alc_indices:
    charge_a, charge_b = charges[alc_idx]
    ewald_self_this += -ke*beta_smooth/np.sqrt(np.pi) * ((1-this_lmbda)*charge_a**2 + (this_lmbda)*charge_b**2)
    for i_for, for_lmbda in enumerate(for_lmbdas):
        ewald_self_for[i_for] += -ke*beta_smooth/np.sqrt(np.pi) * ((1-for_lmbda)*charge_a**2 + (for_lmbda)*charge_b**2)
print("self-self Ewald correction: {} (kJ/mol)".format(ewald_self_this))

box = univ.dimensions[:3] / 10.0

n_waves = 7
waves = []
print("generating wave vectors for n_waves: {}!".format(n_waves))
for n_wave in range(0,n_waves+1):
    waves.append(-n_wave)
    waves.append(n_wave)
waves = np.unique(waves)
wave_vectors = np.array([p for p in itertools.product(waves, repeat=3)])
print("  ...{} wave vectors generated!".format(wave_vectors.shape[0]))

n_frames = univ.trajectory.n_frames
#n_frames = 1
my_diffs = np.zeros((n_frames, len(for_lmbdas), 6))


for i_frame in range(n_frames):
    if i_frame % 100 == 0:
        print("frame {} of {}".format(i_frame, n_frames))
    univ.trajectory[i_frame]
    #univ.trajectory[250]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    coul_this = 0.0
    coul_for = [0.0 for for_lmbda in for_lmbdas]
    ewald_sr_this = 0.0
    ewald_sr_for = [0.0 for for_lmbda in for_lmbdas]
    ewald_lr_this = 0.0
    ewald_lr_for = [0.0 for for_lmbda in for_lmbdas]

    my_diffs[i_frame, :, 0] = univ.trajectory.time

    for i in alc_indices:
        incl_indices = np.setdiff1d(atm_indices, excls[i])
        atm_i = univ.atoms[i]
        pos_i = atm_i.position

        charge_i_a, charge_i_b = charges[i]

        for j in atm_indices:
            atm_j = univ.atoms[j]
            if j in alc_indices:
                if j <= i:
                    continue
                charge_j_a, charge_j_b = charges[j]
            else:
                charge_j_a = charge_j_b = atm_j.charge
                
            pos_j = atm_j.position

            r = np.linalg.norm(pos_i - pos_j)

            coul_this += ke*(1/r)*( charge_i_a*charge_j_a*(1-this_lmbda) + charge_i_b*charge_j_b*(this_lmbda) )
            ewald_sr_this += ke*(1/r)*erfc(beta_smooth*r) * \
                            ( charge_i_a*charge_j_a*(1-this_lmbda) + charge_i_b*charge_j_b*(this_lmbda) )
            for i_for, for_lmbda in enumerate(for_lmbdas):
                coul_for[i_for] += ke*(1/r)*( charge_i_a*charge_j_a*(1-for_lmbda) + charge_i_b*charge_j_b*(for_lmbda))
                ewald_sr_for[i_for] += ke*(1/r)*erfc(beta_smooth*r) * \
                                ( charge_i_a*charge_j_a*(1-for_lmbda) + charge_i_b*charge_j_b*(for_lmbda) )   

            for this_wave in wave_vectors:
                if np.array_equal(this_wave, np.zeros(3)):
                    continue
                k = this_wave * 2 * np.pi / box

                s_k_this = (1-this_lmbda)*(charge_i_a)*np.exp(1j*np.dot(k,pos_i)) + \
                           (this_lmbda)*(charge_i_b)*np.exp(1j*np.dot(k,pos_i)) + \
                           (1-this_lmbda)*(charge_j_a)*np.exp(1j*np.dot(k,pos_j)) + \
                           (this_lmbda)*(charge_j_b)*np.exp(1j*np.dot(k,pos_j))

                k_sq = np.dot(k,k)
                r_ij = pos_i - pos_j
                exp1 = np.exp(1j * np.dot(k, r_ij))
                exp2 = np.exp(-k_sq/(2*beta_smooth)**2)
                ewald_lr_this += (2*np.pi*ke / box.prod()) * (1/k_sq) * np.abs(s_k_this)**2 * exp2

                for i_for, for_lmbda in enumerate(for_lmbdas):
                    s_k_for_a = (charge_i_a)*np.exp(1j*np.dot(k,pos_i)) + \
                                (charge_j_a)*np.exp(1j*np.dot(k,pos_j))
                    
                    s_k_for_b = (charge_i_b)*np.exp(1j*np.dot(k,pos_i)) + \
                                (charge_j_b)*np.exp(1j*np.dot(k,pos_j))
                    s_k_for_b = 0.0+0j
                    ewald_lr_for[i_for] += (2*np.pi*ke / box.prod()) * (1/k_sq) * exp2 * ((1-for_lmbda)*np.abs(s_k_for_a)**2 + (for_lmbda)*np.abs(s_k_for_b)**2)


    my_diffs[i_frame, :, 1] = r
    for i_for in range(n_for_lmbdas):
        my_diffs[i_frame, i_for, 2] = coul_for[i_for] - coul_this
        my_diffs[i_frame, i_for, 3] = ewald_sr_for[i_for] - ewald_sr_this
        my_diffs[i_frame, i_for, 4] = ewald_lr_for[i_for] - ewald_lr_this + ewald_self_for[i_for] - ewald_self_this
        my_diffs[i_frame, i_for, 5] = my_diffs[i_frame, i_for, 3:5].sum() 

