
from mdtools import dr

gmx_ds = dr.loadXVG('dhdl.xvg')
my_ds = np.loadtxt('my_diffs.dat')

lmbdas = [ 0.0,0.25, 0.50, 0.75, 1.00]
this_lmbda = 1.0

for i, for_lmbda in enumerate(lmbdas):
    gmx_data = gmx_ds.data[for_lmbda][::10]
    my_data = my_ds[:, i+1]
    fig, ax1 = plt.subplots()
    ax1.plot(gmx_data, 'o', label=r'GMX output', markersize=10)
    ax1.plot(my_ds[:,0], my_data, '--o', label=r'Calculated', linewidth=2)
    ax1.set_ylabel(r'$\Delta U$ (kJ/mol)')
    ax1.set_title(r'$\lambda_0={}; \; \lambda_1={}$'.format(this_lmbda, for_lmbda))
    ax1.set_xlim(0,100)
    ax1.set_xlabel('time (ps)')
    ax1.legend()

    #plt.show()
    left, bottom, width, height = [0.25, 0.6, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.plot(gmx_data-my_data, '-o')
    ax2.set_ylabel(r'$\Delta \Delta U$ (kJ/mol)')
    ax2.set_xlim(0,100)

    plt.show()