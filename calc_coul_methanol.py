## Do vdw calculation from md trajectory##

from __future__ import division, print_function

import MDAnalysis
import numpy as np

import argparse
try:
    u_for_prev = u_for
except:
    pass

from math import erf, erfc

from mdtools import dr
import itertools
#from IPython import embed


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
parser.add_argument('--n-wave-images', type=int, default=10,
                    help='Number of wave numbers of each box vector to consider when' \
                    'calculating long-range ewald coulombic interactions; e.g if 1 will consider all 26 additional' \
                    'adjacent box reciprocal wavelengths images (by shifting each component of the box vector -1, 0, or 1).')


#args = parser.parse_args()
args = parser.parse_args(['-s', 'run.tpr', '-f', 'traj.trr', '--lmbda', '0.0', '--for-lmbda', '1.0'])
univ = MDAnalysis.Universe(args.tprfile, args.trajfile)

# box vector
box = univ.dimensions[:3] / 10.0

outfile = args.outfile

this_lmbda = args.lmbda

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

'''
charges = {0:(0.12010, 0.0),
           1:(-0.60030, 0.0),
           2:(0.02770, 0.0),
           3:(0.02770, 0.0),
           4:(0.02770, 0.0),
           5:(0.39710, 0.0)}


HW_charge = 0.417
OW_charge = -0.834

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
'''
charges = {0:(0.0, 0.0),
           1:(0.0, 0.0),
           2:(0.25, 0.0),
           3:(0.25, 0.0),
           4:(0.0, 0.0),
           5:(-0.5, 0.0)}


HW_charge = 0.417
OW_charge = -0.834

excls = {0: (0,1,2,3,4,5),
         1: (0,1,2,3,4,5),
         2: (0,1,2,3,4),
         3: (0,1,2,3,4,5),
         4: (0,1,2,3,4,5),
         5: (0,1,3,4,5)}


# 1-4 pairs

pairs = {0: (),
         1: (),
         2: (),
         3: [],
         4: [],
         5: []}

ke = 138.9354859 # conv factor to get kJ/mol

fudge = 0.83333333
## Do vdw calculation from md trajectory##

rtol = 1e-6
rc = 1.0

# determine ewald stuff
xvals = np.linspace(0,5,1000000)
vals = np.zeros_like(xvals)
for i, x in enumerate(xvals):
    vals[i] = erfc(x)

etol_ind = np.argmax(vals < rtol)
beta_smooth = xvals[etol_ind] / rc
print("smoothing width: {}".format(beta_smooth))

ewald_self_this = 0.0
ewald_self_for = [0.0 for for_lmbda in for_lmbdas]

for alc_idx in alc_indices:
    charge_a, charge_b = charges[alc_idx]
    ewald_self_this += -ke*beta_smooth/np.sqrt(np.pi) * ((1-this_lmbda)*charge_a**2 + (this_lmbda)*charge_b**2)
    for i_for, for_lmbda in enumerate(for_lmbdas):
        ewald_self_for[i_for] += -ke*beta_smooth/np.sqrt(np.pi) * ((1-for_lmbda)*charge_a**2 + (for_lmbda)*charge_b**2)
print("self-self Ewald correction: {} (kJ/mol)".format(ewald_self_this))

box = univ.dimensions[:3] / 10.0

n_waves = args.n_wave_images
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
my_diffs = np.zeros((n_frames, n_for_lmbdas, 6))


for i_frame in range(n_frames):
    if i_frame % 100 == 0:
        print("frame {} of {}".format(i_frame, n_frames))
    univ.trajectory[i_frame]
    # Get positions in nm
    univ.atoms.positions = univ.atoms.positions / 10.0
    coul_this = 0.0
    coul_for = [0.0 for for_lmbda in for_lmbdas]
    ewald_sr_this = 0.0
    ewald_sr_for = [0.0 for for_lmbda in for_lmbdas]
    ewald_lr_this = 0.0
    ewald_lr2 = 0.0
    ewald_lr_for = [0.0 for for_lmbda in for_lmbdas]

    coul_14_this = 0.0
    coul_14_for = [0.0 for for_lmbda in for_lmbdas]

    my_diffs[i_frame, :, 0] = univ.trajectory.time
    
    for this_wave in wave_vectors:
        if np.array_equal(this_wave, np.zeros(3)):
            continue

        k = this_wave * 2 * np.pi / box
        k_sq = np.dot(k,k)
        exp2 = np.exp(-k_sq/(2*beta_smooth)**2)
        s_k = 0.0
        for i in atm_indices:
            atm_i = univ.atoms[i]
            pos_i = atm_i.position
            charge_i, blah = charges[i]
            s_k += charge_i * np.exp(1j*np.dot(k, pos_i))
            #print("s_k: {}".format(s_k))

        ewald_lr2 += (2*np.pi*ke / box.prod()) * (1/k_sq) * np.abs(s_k)**2 * exp2
    
    for i in alc_indices:
        print("atom {} of {}".format(i, alc_indices[-1]))
        incl_indices = np.setdiff1d(atm_indices, excls[i])
        pair_indices = pairs[i]
        atm_i = univ.atoms[i]
        pos_i = atm_i.position

        charge_i_a, charge_i_b = charges[i]

        # all atoms must be considered for LR
        for j in atm_indices:
            atm_j = univ.atoms[j]
            if j in alc_indices:
                if j <= i:
                    continue
                charge_j_a, charge_j_b = charges[j]
            else:
                charge_j_a = charge_j_b = atm_j.charge
            pos_j = atm_j.position
            #print("doing other atom {}".format(j))

            #if j in pair_indices:
            #    continue
            ## EWALD SR
            r = np.linalg.norm(pos_i - pos_j)
            if j in incl_indices:
                coul_this += ke*(1/r)*( charge_i_a*charge_j_a*(1-this_lmbda) + charge_i_b*charge_j_b*(this_lmbda) )
                ewald_sr_this += ke*(1/r)*erfc(beta_smooth*r) * \
                                ( charge_i_a*charge_j_a*(1-this_lmbda) + charge_i_b*charge_j_b*(this_lmbda) )
                for i_for, for_lmbda in enumerate(for_lmbdas):
                    coul_for[i_for] += ke*(1/r)*( charge_i_a*charge_j_a*(1-for_lmbda) + charge_i_b*charge_j_b*(for_lmbda))
                    ewald_sr_for[i_for] += ke*(1/r)*erfc(beta_smooth*r) * \
                                    ( charge_i_a*charge_j_a*(1-for_lmbda) + charge_i_b*charge_j_b*(for_lmbda) )  
        
            #Ewald LR
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
                
                fac1 = np.pi * this_wave / beta_smooth
                exp_1 = np.exp(-(np.dot(fac1, fac1) ))
                fac2 = np.dot(2*np.pi*1j*this_wave, r_ij)
                exp_2 = np.exp(fac2)
                m_sq = np.dot(this_wave, this_wave)
                #ewald_lr2 += ke/(2*np.pi*box.prod()) * charge_i_a * charge_j_b * exp_1 * exp_2 / m_sq

                for i_for, for_lmbda in enumerate(for_lmbdas):
                    s_k_for_a = (charge_i_a)*np.exp(1j*np.dot(k,pos_i)) + \
                                (charge_j_a)*np.exp(1j*np.dot(k,pos_j))
                    
                    s_k_for_b = (charge_i_b)*np.exp(1j*np.dot(k,pos_i)) + \
                                (charge_j_b)*np.exp(1j*np.dot(k,pos_j))
                    s_k_for_b = 0.0+0j
                    ewald_lr_for[i_for] += (2*np.pi*ke / box.prod()) * (1/k_sq) * exp2 * ((1-for_lmbda)*np.abs(s_k_for_a)**2 + (for_lmbda)*np.abs(s_k_for_b)**2)

        # Treat 1-4 pair interactions with normal coulombic 1/r potential
        for j in pair_indices:
            if j <= i:
                continue
            # sanity - for methanol
            assert i in alc_indices
            atm_j = univ.atoms[j]
            pos_j = atm_j.position
            charge_j_a, charge_j_b = charges[j]

            r = np.linalg.norm(pos_i - pos_j)

            coul_14_this += fudge*ke*(1/r)*( (1-this_lmbda)*charge_i_a*charge_j_a + (this_lmbda)*charge_i_b*charge_j_b )

            for i_for, for_lmbda in enumerate(for_lmbdas):
                coul_14_for[i_for] += fudge*ke*(1/r)*( (1-for_lmbda)*charge_i_a*charge_j_a + (for_lmbda)*charge_i_b*charge_j_b )


    for i_for in range(n_for_lmbdas):
        my_diffs[i_frame, i_for, 1] = coul_for[i_for] - coul_this
        my_diffs[i_frame, i_for, 2] = ewald_sr_for[i_for] - ewald_sr_this
        my_diffs[i_frame, i_for, 3] = ewald_lr_for[i_for] - ewald_lr_this
        my_diffs[i_frame, i_for, 4] = coul_14_for[i_for] - coul_14_this
        my_diffs[i_frame, i_for, 5] = my_diffs[i_frame, i_for, 3:6].sum() + ewald_self_for[i_for] - ewald_self_this


#header = 'time(ps)        r(nm)        Coul_simple(kJ/mol)        Coul_with_periodic(kJ/mol)         U_dir        U_rec      U_pairs    U_ewald_tot'
#np.savetxt('my_diffs.dat', my_diffs)

#plt.plot(my_diffs[:,0], my_diffs[:,5], label='calc')
#plt.plot(my_diffs[:,0], gmx_wave, 'o', label='gmx')
#plt.legend()


