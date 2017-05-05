
from mdtools import dr

lmbda_for = 1.0
fe_ds = dr.loadXVG('fe.xvg')
fe_dat = fe_ds.data[lmbda_for]
plt.plot(-fe_dat, '-o', label=r'(Gromacs FE)')
plt.plot(my_diffs[:,0], -my_diffs[:,1], '-o', label=r'direct calc')

plt.legend()
plt.show()
