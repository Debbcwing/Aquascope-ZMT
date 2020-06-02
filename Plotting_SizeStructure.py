import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import netCDF4
from scipy.interpolate import make_interp_spline


################
### Figure 3 ###
################
"""
Fig 3 shows results of phytoplantkon size distribution of overall long-term and steady state solutions
for detailed set-up please refer to UptkGraz_10yr.py
"""
Pnum = 150  # specify the nos of phytoplankton size groups
Ynum = 10  # no. of modeling year
nut_rgm = [1, 15, 50]  # µmol NL-1

# Import NetCDF data (These are the only dataset for the results)
os.chdir('wd')
out1 = netCDF4.Dataset('StdRun1.nc')
out15 = netCDF4.Dataset('StdRun15.nc')
out50 = netCDF4.Dataset('StdRun50.nc')

# set up arrays to store extracted data
P_sc = np.zeros((3650, Pnum, 3, 3))  # time; Psize class; nutrient regime; mixing regime

# Extract (in loop)
outlist = [out1, out15, out50]
for a in range(len(nut_rgm)):
    print(outlist[a].variables['Nut'][:].shape,
          outlist[a].variables['Zoo'][:].shape,
          outlist[a].variables['Phy'][:].shape)  # (3650, 3) (3650, 2, 3) (3650, 150, 3)
    P_sc[:, :, a, 0] = outlist[a].variables['Phy'][:, :, 0]
    P_sc[:, :, a, 1] = outlist[a].variables['Phy'][:, :, 2]
    P_sc[:, :, a, 2] = outlist[a].variables['Phy'][:, :, 1]  # simlist str = (time; concentration; MixRgm)

"""
Plotting
"""
Nutlab = ['Oligotrophic', 'Eutrophic', 'Hypertrophic']
Mixlist = ['Constant', 'Medium mixing', 'High mixing']
lw, fs = 1.5, 11
Psize = np.logspace(0, 2, 150)  # define phytoplankton size array
Zsize = [5, 200]

# Gross phytoplankton population, gross growth and gross grazing (sum 10 yrs)
fig, axs = plt.subplots(3, 3, figsize=(8, 6), sharex='all', sharey='all')
plt.subplots_adjust(hspace=0.1, wspace=0.15)
for i in range(3):
    for j in range(3):  # mixing
        axs[j, i].plot(Psize, np.sum(P_sc, axis=0)[:, i, j], color='#454242')
        axs[j, i].bar(Psize, np.sum(P_sc, axis=0)[:, i, j], color='#454242')
        axs[j, i].loglog()
        axs[j, i].axis([1e0, 1e2, 1e-2, 1e5])
        axs[0, i].set_title(Nutlab[i], fontweight='bold', fontsize=10)
        axs[j, i].set_xticks([1e0, 1e1, 1e2])
        axs[j, i].set_xticklabels([1, 10, 100])
        axs[j, i].set_yticks([1e-2, 1e0, 1e2, 1e5])
        axs[j, 2].yaxis.set_label_position('right')
        axs[j, 2].yaxis.set_ticks_position('right')
        axs[j, i].yaxis.set_ticks_position('both')
        axs[j, 2].set_ylabel(Mixlist[j], fontweight='bold', fontsize=10, rotation=270, labelpad=13)
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        axins2 = inset_axes(axs[j, i], width=0.65, height=0.75, borderpad=0.65)
        axins2.bar(Psize, np.sum(P_sc[3285:], axis=0)[:, i, j], color='#454242')
        axins2.loglog()
        axins2.axis([1e0, 1e2, 1e-2, 1e4])
        axins2.set_xticks([1e0, 1e1, 1e2])
        axins2.set_xticklabels([1, 10, 100])
        axins2.set_yticks([1e-2, 1e0, 1e2, 1e4])
        axins2.tick_params(labelsize=7)
axs[2, 1].set_xlabel('Phytoplankton size, $S_{i}$ [µm ESD]', fontsize=fs)
axs[1, 0].set_ylabel('Phytoplankton biomass, $P_{i}$ [µmol NL$^{-1}$]', fontsize=fs)
#fig.savefig(fname='/Users/szewing/Desktop/PhD_work/SizeMod/Plotting_report/Fig3.png')
