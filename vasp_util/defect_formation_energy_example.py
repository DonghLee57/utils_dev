# python3.X
# charge defect: Co interstitial in Si
#
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
font = {'size': 15}
mpl.rc('font', **font)
mpl.rcParams['xtick.major.size'] = font['size']
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.size'] = font['size']
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['axes.linewidth'] = 1

ref_Si = -1171.67153566
vbm_Si = 5.6167
cbm_Si = 6.2318
u_Co = -14.08523189/2

Co_int_T_0 = -1177.41232600 #  0 neutral
Co_int_T_1 = -1183.28660572 # +1 charged
Co_int_T_2 = -1188.85461681 # +2 charged
Co_int_Tm1 = -1171.33722155 # -1 charged

x = [0,cbm_Si-vbm_Si]
Co0  = [Co_int_T_0 - (ref_Si + u_Co),              Co_int_T_0 - (ref_Si + u_Co)]
Co_1 = [Co_int_T_1 - (ref_Si + u_Co) + (vbm_Si),   Co_int_T_1 - (ref_Si + u_Co) + (cbm_Si)]
Co_2 = [Co_int_T_2 - (ref_Si + u_Co) + 2*(vbm_Si), Co_int_T_2 - (ref_Si + u_Co) + 2*(cbm_Si)]
Com1 = [Co_int_Tm1 - (ref_Si + u_Co) - (vbm_Si),   Co_int_Tm1 - (ref_Si + u_Co) - (cbm_Si)]

f_0 = np.poly1d(np.polyfit(x,Co0,1))
f_1 = np.poly1d(np.polyfit(x,Co_1,1))
f_2 = np.poly1d(np.polyfit(x,Co_2,1))
fm1 = np.poly1d(np.polyfit(x,Com1,1))
xx = np.linspace(0, cbm_Si-vbm_Si)
tot = np.array([f_0(xx),f_1(xx),f_2(xx),fm1(xx)])
fig, ax = plt.subplots()

ax.plot(x, Co_2,':',label='+2')
ax.plot(x, Co_1,':',label='+1')
ax.plot(x, Co0,':',label='0')
ax.plot(x, Com1,':',label=r'$\mathrm{-}$1')
ax.plot(xx, np.min(tot,axis=0),'r',label='min')
ax.set_xlim(x)
ax.set_ylim([1,2])
ax.set_xlabel('Fermi level (eV)')
ax.set_ylabel('Formation energy (eV)')
plt.legend()
plt.tight_layout()
plt.show()
