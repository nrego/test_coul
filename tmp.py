fig, ax = plt.subplots()

rects = ax.bar([0,2,4], [-dgs_bulk_coul[0]+dgs_vac_coul[0],-dgs_bulk_vdw[0], -dgs_bulk_coul[0]+dgs_vac_coul[0]-dgs_bulk_vdw[0]], width=1,yerr=[0.026209213908471179,  0.059392227117696116, 0.064918098675947042], color='b', label='my results', error_kw={'ecolor':'k', 'elinewidth':10, 'capsize':10})

rects = ax.bar([1,3,5], [-8.638604854916057, 2.801256331594139, -5.837348523321918], yerr=[0.016774,0.016774,0.016774], width=1,color='r', label='published', error_kw={'ecolor':'k', 'elinewidth':10, 'capsize':10})

ax.set_xticks([1,3,5])
ax.set_xticklabels((r'$\Delta\Delta G_{\rm{coul}}$', r'$\Delta G_{\rm{VdW}}$',r'$\Delta G_{\rm{hyd}}$'))
#plt.title(r'$\Delta G_{\rm{vdw}}^{\rm{solv}}$')
plt.legend()
plt.ylabel(r'$\Delta G \; (k_B T)$')