from mdtools import dr
import numpy as np


plt.clf()
pme = dr.loadXVG('pme.xvg').data[1.0]
no_pme = dr.loadXVG('no_pme.xvg').data[1.0]
#ewald_potential = np.loadtxt('ewald_energy.xvg', comments=['@', '#'])
pme_potential = np.loadtxt('ewald_energy.xvg', comments=['@', '#'])

my_diffs = np.loadtxt('my_diffs.dat')

#plt.plot(my_diffs[:,0], -my_diffs[:,2], label='calculated coulombic potential')
plt.plot(my_diffs[:,0], -my_diffs[:,3], label='calculated 1/r coul potential (with periodic images)')
plt.plot(my_diffs[:,0], -my_diffs[:,4], label='calculated ewald SR (direct)')
#plt.plot(my_diffs[:,0], -my_diffs[:,5], label='calculated ewald LR (recip)')

#plt.plot(ewald_potential[:,0], -ewald_potential[:,1::2], 'o', markersize=10, label='total potential energy (gmx output w/ ewald)')
plt.plot(pme_potential[:,0::2], -pme_potential[:,1::2], 'o', label='total ewald energy (gmx output w/ pme)')

plt.plot(no_pme[::2], 'o', label='cut-off (gmx fe output)')
plt.plot(pme[::2], 'o', label='pme (gmx fe output)')

plt.legend()
plt.xlim(0,2)
plt.ylim(0,10)
plt.savefig('blah.png')