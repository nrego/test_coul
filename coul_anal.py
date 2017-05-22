from mdtools import dr
import numpy as np


plt.clf()
fig, ax = plt.subplots()
pme = dr.loadXVG('pme.xvg').data[1.0]
pme_energy = np.loadtxt('pme_energy.xvg', comments=['@','#'])
no_pme = dr.loadXVG('no_pme.xvg').data[1.0]
no_pme_energy = np.loadtxt('no_pme_energy.xvg', comments=['@','#'])

recip = np.loadtxt('pme_energy_lr.xvg', comments=['@', '#'])
gmx_ewald_lr = recip[:,1] - recip[0,1]
gmx_wave = gmx_ewald_lr
ewald_sr = np.loadtxt('ewald_sr.xvg', comments=['@', '#'])

my_diffs = np.loadtxt('my_diffs.dat')

#plt.plot(my_diffs[:,0], my_diffs[:,3], label='calculated 1/r coul potential (with periodic images)')

ax.plot(my_diffs[:,0], my_diffs[:,2], label='calculated 1/r coulombic potential', linewidth=4)
#ax.plot(my_diffs[:,0], my_diffs[:,4], label='calculated ewald SR potential', linewidth=4)
#ax.plot(-no_pme, 'o', markersize=10, label=r'cut-off (gmx fe output, $\Delta U=U_0-U_1$)')
#ax.plot(no_pme_energy[:,0], no_pme_energy[:,1], 'o', markersize=6, label=r'cut-off, potential energy (gmx energy output)')

#ax.plot(my_diffs[:,0], ewald_sr[:,1], 'o', markersize=4, label='ewald SR potential (gmx energy output)')

ax.plot(my_diffs[:,0], my_diffs[:,5], label='calculated ewald LR (recip)')
ax.plot(my_diffs[:,0], ewald_lr, 'o', markersize=10, label='ewald LR (gmx output)')
#ax.plot(my_diffs[:,0], ewald_lr-my_diffs[:,5], label='difference (GMX output minus calc)')
#ax.plot(my_diffs[:,0], ewald_lr/my_diffs[:,5], label='quotient (GMX output by calc)')

#ax.plot(-pme, 'o', markersize=10, label='PME (gmx fe output, $\Delta U=U_0-U_1$)')
#ax.plot(pme_energy[:,0], pme_energy[:,1], 'o', markersize=6, label=r'PME, potential energy (gmx energy output)')

ax.set_xticklabels([0.0, 0.5, 1.0, 0.5, 0.0])
plt.legend(loc=0)
plt.xlim(0,2)
plt.ylim(0,5)
plt.savefig('blah.png')