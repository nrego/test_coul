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

from math import erf, erfc
from scipy.special import erfinv

for_lmbdas = [1.0]
n_for_lmbdas = len(for_lmbdas)
this_lmbda = 0.0

atm_indices = np.arange(univ.atoms.n_atoms)
alc_indices = np.arange(872, 888)
#alc_indices = [877]

excls = {i: None for i in alc_indices}

excls[872]=[853, 868, 869, 870, 871, 872, 873, 874, 
875, 876, 877, 878, 879, 880, 884, 888, 889, 890, 891, 892]
excls[873]=[868, 870, 871, 872, 873, 874, 875, 876, 
877, 878, 888, 889, 890]
excls[874]=[868, 870, 871, 872, 873, 874, 875, 876, 
877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 
889, 890]
excls[875]=[870, 872, 873, 874, 875, 876, 877, 878, 
879, 880, 884, 888]
excls[876]=[870, 872, 873, 874, 875, 876, 877, 878, 
879, 880, 884, 888]
excls[877]=[870, 872, 873, 874, 875, 876, 877, 878, 
879, 880, 884, 888]
excls[878]=[870, 872, 873, 874, 875, 876, 877, 878, 
879, 880, 881, 882, 883, 884, 885, 886, 887, 888]
excls[879]=[872, 874, 875, 876, 877, 878, 879, 880, 
881, 882, 883, 884, 885, 886, 887]
excls[880]=[872, 874, 875, 876, 877, 878, 879, 880, 
881, 882, 883, 884, 885, 886, 887]
excls[881]=[874, 878, 879, 880, 881, 882, 883, 884]
excls[882]=[874, 878, 879, 880, 881, 882, 883, 884]
excls[883]=[874, 878, 879, 880, 881, 882, 883, 884]
excls[884]=[872, 874, 875, 876, 877, 878, 879, 880, 
881, 882, 883, 884, 885, 886, 887]
excls[885]=[874, 878, 879, 880, 884, 885, 886, 887]
excls[886]=[874, 878, 879, 880, 884, 885, 886, 887]
excls[887]=[874, 878, 879, 880, 884, 885, 886, 887]

pairs = {i: None for i in alc_indices}
pairs[872] = [853, 869, 879, 880, 884, 891, 892]                # CA
pairs[873] = [868, 871, 875, 876, 877, 878, 889, 890]           # HA
pairs[874] = [868, 871, 881, 882, 883, 885, 886, 887, 889, 890] # CB
pairs[875] = [870, 879, 880, 884, 888]                          # HB1
pairs[876] = [870, 879, 880, 884, 888]                          # HB2
pairs[877] = [870, 888]                                         # HB3
pairs[878] = [870, 888]                                         # CG
pairs[879] = [881, 882, 883, 885, 886, 887]                     # HG
pairs[880] = [885, 886, 887]                                    # CD1
pairs[881] = [884]                                              # HD11
pairs[882] = [884]                                              # HD12
pairs[883] = [884]                                              # HD13
pairs[884] = []                                                 # CD2
pairs[885] = []                                                 # HD21
pairs[886] = []                                                 # HD22
pairs[887] = []                                                 # HD23


alc_info = {i: None for i in alc_indices}

alc_info[872]={'typeA':  2, 'typeB':  2,  'mA': 1.20100e+01, 'qA':-5.18000e-02, 'mB': 1.20100e+01, 'qB': 3.37000e-02, 'resind':   62, 'atomnumber':  6}
alc_info[873]={'typeA':  8, 'typeB':  8,  'mA': 1.00800e+00, 'qA': 9.22000e-02, 'mB': 1.00800e+00, 'qB': 8.23000e-02, 'resind':   62, 'atomnumber':  1}
alc_info[874]={'typeA':  2, 'typeB':  2,  'mA': 1.20100e+01, 'qA':-1.10200e-01, 'mB': 1.20100e+01, 'qB':-1.82500e-01, 'resind':   62, 'atomnumber':  6}
alc_info[875]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 4.57000e-02, 'mB': 1.00800e+00, 'qB': 6.03000e-02, 'resind':   62, 'atomnumber':  1}
alc_info[876]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 4.57000e-02, 'mB': 1.00800e+00, 'qB': 6.03000e-02, 'resind':   62, 'atomnumber':  1}
alc_info[877]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 0.00000e+00, 'mB': 1.00800e+00, 'qB': 6.03000e-02, 'resind':   62, 'atomnumber':  1}
alc_info[878]={'typeA':  2, 'typeB':  2,  'mA': 1.20100e+01, 'qA': 3.53100e-01, 'mB': 1.20100e+01, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  6}
alc_info[879]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA':-3.61000e-02, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}
alc_info[880]={'typeA':  2, 'typeB':  2,  'mA': 1.20100e+01, 'qA':-4.12100e-01, 'mB': 1.20100e+01, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  6}
alc_info[881]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 1.00000e-01, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}
alc_info[882]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 1.00000e-01, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}
alc_info[883]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 1.00000e-01, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}
alc_info[884]={'typeA':  2, 'typeB':  2,  'mA': 1.20100e+01, 'qA':-4.12100e-01, 'mB': 1.20100e+01, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  6}
alc_info[885]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 1.00000e-01, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}
alc_info[886]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 1.00000e-01, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}
alc_info[887]={'typeA':  4, 'typeB':  4,  'mA': 1.00800e+00, 'qA': 1.00000e-01, 'mB': 1.00800e+00, 'qB': 0.00000e+00, 'resind':   62, 'atomnumber':  1}

other_charges = {
    853: -0.0518,  # CA(-)
    868:  0.5973,  # C(-)
    869: -0.5679,  # O(-)
    870: -0.4157,  # N
    871:  0.2719,  # H
    888:  0.5973,  # C
    889: -0.5679,  # O
    890: -0.4157,  # N(+),
    891:  0.2719,  # H(+)
    892:  0.0213   # CA(+)
}
rtol = 1e-6
rc = 1.0

# determine ewald stuff
beta_smooth = erfinv(1-rtol) / rc
print("smoothing width: {}".format(beta_smooth))

ke = 138.9354859 # coul conv factor to get kJ/mol

fudge = 0.8333

ewald_self_this = 0.0
ewald_self_for = [0.0 for for_lmbda in for_lmbdas]


## Ewald LR self-correction
for alc_idx in alc_indices:
    atm_info = alc_info[alc_idx]
    charge_a, charge_b = atm_info['qA'], atm_info['qB']
    ewald_self_this += -ke*beta_smooth/np.sqrt(np.pi) * ((1-this_lmbda)*charge_a**2 + (this_lmbda)*charge_b**2)
    for i_for, for_lmbda in enumerate(for_lmbdas):
        ewald_self_for[i_for] += -ke*beta_smooth/np.sqrt(np.pi) * ((1-for_lmbda)*charge_a**2 + (for_lmbda)*charge_b**2)
print("self-self Ewald correction: {} (kJ/mol)".format(ewald_self_this))

box = univ.dimensions[:3] / 10.0

n_waves = 8
waves = []
print("generating wave vectors for n_waves: {}!".format(n_waves))
for n_wave in range(0,n_waves+1):
    waves.append(-n_wave)
    waves.append(n_wave)
waves = np.unique(waves)
wave_vectors = np.array([p for p in itertools.product(waves, repeat=3)])
print("  ...{} wave vectors generated!".format(wave_vectors.shape[0]))

n_frames = univ.trajectory.n_frames
n_frames = 1
my_diffs = np.zeros((n_frames, len(for_lmbdas), 6))

for i_frame in range(n_frames):
    if i_frame % 100 == 0:
        print("frame {} of {}".format(i_frame, n_frames))
    #univ.trajectory[i_frame]
    univ.trajectory[-1]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    ## Regular 1/r coul interactions with cut-off
    coul_this = 0.0
    coul_for = [0.0 for for_lmbda in for_lmbdas]
    ## 1-4 pair interactions (treated as regular 1/r coul interactions)
    coul_14_this = 0.0
    coul_14_for = [0.0 for for_lmbda in for_lmbdas]
    ## Short-ranged part of ewald sum (depends on rc, the cutoff, and rtol)
    ewald_sr_this = 0.0
    ewald_sr_for = [0.0 for for_lmbda in for_lmbdas]
    ## Long-ranged (fourier) part of ewald sum (also depends on rc and rtol, which together determine beta_smooth)
    ewald_lr_this = 0.0
    ewald_lr_for = [0.0 for for_lmbda in for_lmbdas]

    my_diffs[i_frame, :, 0] = univ.trajectory.time

    for i in alc_indices:
        incl_indices = np.setdiff1d(atm_indices, excls[i])
        pair_indices = pairs[i]
        atm_i = univ.atoms[i]
        pos_i = atm_i.position
        atm_i_info = alc_info[i]

        charge_i_a, charge_i_b = atm_i_info['qA'], atm_i_info['qB']

        for j in incl_indices:
            continue
            assert j not in pair_indices, "double-counting a pair!"
            atm_j = univ.atoms[j]
            if j in alc_indices:
                if j <= i:
                    continue
                atm_j_info = alc_info[j]
                charge_j_a, charge_j_b = atm_j_info['qA'], atm_j_info['qB']
            else:
                #charge_j_a = charge_j_b = np.round(atm_j.charge, decimals=4)
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

        # 1-4 pairs; regular coulomb only
        this_contrib = 0.0
        for_contrib = 0.0
        for j in pair_indices:
            atm_j = univ.atoms[j]
            if j in alc_indices:
                if j <= i:
                    continue
                atm_j_info = alc_info[j]
                charge_j_a, charge_j_b = atm_j_info['qA'], atm_j_info['qB'] 
            else:
                #charge_j_a = charge_j_b = np.round(atm_j.charge, decimals=4)
                #charge_j_a = charge_j_b = other_charges[j]
                charge_j_a = charge_j_b = atm_j.charge
            pos_j = atm_j.position

            r = np.linalg.norm(pos_i - pos_j)

            coul_14_this += fudge*ke*(1/r)*( charge_i_a*charge_j_a*(1-this_lmbda) + charge_i_b*charge_j_b*(this_lmbda) )
            this_contrib += fudge*ke*(1/r)*( charge_i_a*charge_j_a*(1-this_lmbda) + charge_i_b*charge_j_b*(this_lmbda) )
            for i_for, for_lmbda in enumerate(for_lmbdas):
                coul_14_for[i_for] += fudge*ke*(1/r)*( charge_i_a*charge_j_a*(1-for_lmbda) + charge_i_b*charge_j_b*(for_lmbda) )
                for_contrib += fudge*ke*(1/r)*( charge_i_a*charge_j_a*(1-for_lmbda) + charge_i_b*charge_j_b*(for_lmbda) )

        print("14 interactions for atm {}".format(atm_i.name))
        print("  lambda=0:  {}".format(this_contrib))
        print("  lambda=1:  {}".format(for_contrib))

    my_diffs[i_frame, :, 1] = r
    for i_for in range(n_for_lmbdas):
        my_diffs[i_frame, i_for, 2] = coul_for[i_for] - coul_this
        my_diffs[i_frame, i_for, 3] = ewald_sr_for[i_for] - ewald_sr_this
        my_diffs[i_frame, i_for, 4] = ewald_lr_for[i_for] - ewald_lr_this + ewald_self_for[i_for] - ewald_self_this
        my_diffs[i_frame, i_for, 5] = my_diffs[i_frame, i_for, 3:5].sum() 

