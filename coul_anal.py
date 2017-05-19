from mdtools import dr
import numpy as np


plt.clf()
pme = dr.loadXVG('pme.xvg').data[1.0]
no_pme = dr.loadXVG('no_pme.xvg').data[1.0]
no_pme_energy = np.loadtxt('no_pme_energy.xvg', comments=['@','#'])

ewald_lr = np.loadtxt('ewald_lr.xvg', comments=['@', '#'])
ewald_sr = np.loadtxt('ewald_sr.xvg', comments=['@', '#'])

my_diffs = np.loadtxt('my_diffs.dat')

plt.plot(my_diffs[:,0], -my_diffs[:,2], label='calculated 1/r coulombic potential')
#plt.plot(my_diffs[:,0], -my_diffs[:,3], label='calculated 1/r coul potential (with periodic images)')
#plt.plot(my_diffs[:,0], -my_diffs[:,4], label='calculated ewald SR (direct)')
plt.plot(my_diffs[:,0], -my_diffs[:,5], label='calculated ewald LR (recip)')


#plt.plot(my_diffs[:,0], -ewald_sr[:,1], 'o', markersize=10, label='Ewald SR (gmx output)')
plt.plot(my_diffs[:,0], (ewald_lr[:,1]-ewald_lr[0,1]), 'o', markersize=10, label='Ewald LR (gmx output)')
#plt.plot(pme_potential[:,0::2], -pme_potential[:,1::2], 'o', markersize=10, label='ewald LR (recip) (gmx energy output)')

plt.plot(no_pme[::2], 'o', markersize=10, label='cut-off (gmx fe output)')
plt.plot(no_pme_energy[:,0], -no_pme_energy[:,1], 'o', label='cut-off (GMX potential energy output)')
#plt.plot(pme[::2], 'o', label='pme (gmx fe output)')

plt.legend()
plt.xlim(0,2)
plt.ylim(0,10)
plt.savefig('blah.png')