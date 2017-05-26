

charge_i = -0.5
charge_j = 0.25
i = 5
j = 3
pos_i = univ.atoms[i].position
pos_j  =univ.atoms[j].position
for this_wave in wave_vectors:
    if np.array_equal(this_wave, np.zeros(3)):
        continue
    k = this_wave * 2 * np.pi / box
    sk = charge_i*np.exp(1j*np.dot(k, pos_i)) + charge_j*np.exp(1j*np.dot(k,pos_j))
    k_sq = np.dot(k,k)
    exp2 = np.exp(-k_sq/(2*beta_smooth)**2)
    ewald_lr += (2*np.pi*ke / box.prod()) * (1/k_sq) * exp2 * np.abs(sk)**2
