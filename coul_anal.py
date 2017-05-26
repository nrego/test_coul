from mdtools import dr
import numpy as np


plt.clf()
fig, ax = plt.subplots()
pme = dr.loadXVG('pme.xvg').data[1.0]
pme_energy = np.loadtxt('pme_energy.xvg', comments=['@','#'])
#no_pme = dr.loadXVG('no_pme.xvg').data[1.0]
#no_pme_energy = np.loadtxt('no_pme_energy.xvg', comments=['@','#'])

pme_lr = np.loadtxt('pme_lr.xvg', comments=['@', '#'])
#pme_lr = recip[:,1] - recip[0,1]
pme_sr = np.loadtxt('pme_sr.xvg', comments=['@', '#'])

#ax.plot(my_diffs[:,1,0], my_diffs[:,1,2], label='calculated 1/r coulombic potential', linewidth=4)
ax.plot(my_diffs[:,1,0], -my_diffs[:,1,3], label='calculated ewald SR potential', linewidth=6)
ax.plot(pme_sr[:,0::2], pme_sr[:,1::2], 'o', markersize=8, label='ewald SR (gmx energy output)')

ax.plot(my_diffs[:,1,0], -my_diffs[:,1,4], label='calculated ewald LR (recip)', linewidth=6)
ax.plot(pme_lr[:,0::2], pme_lr[:,1::2], 'o', markersize=8, label='ewald LR (gmx energy output)')

#ax.plot(my_diffs[:,1,0], my_diffs[:,1,5], label='Calculated Ewald $\Delta U$', linewidth=4)
#ax.plot(pme, 'o', markersize=10, label='PME (gmx fe output, $\Delta U=U_1-U_0$)')
#ax.plot(pme_energy[:,0], pme_energy[:,1], 'o', markersize=6, label=r'PME, potential energy (gmx energy output)')

ax.set_xticklabels([0.0, 0.5, 1.0, 0.5, 0.0])
plt.legend(loc=0)
#plt.xlim(0,2)
plt.ylim(-10,0)
plt.xlabel('r (nm)')
plt.ylabel('U (kJ/mol)')

plt.savefig('blah.png')