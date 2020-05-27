import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4


################
### Figure 4 ###
################
"""
Fig 4 shows the productivity-diversity relationships
for detailed set-up please refer to SizeMatrix.py
"""
Pnum = 150  # specify the nos of phytoplankton size groups
Ynum = 10  # no. of modeling year
nut_rgm = [1, 15, 50]  # µmol NL-1

# Import NetCDF data (These are the only dataset for the results)
os.chdir('/Users/szewing/Desktop/PhD_work/SizeMod/Results/StdRun/')
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

def WMeanSize(b, s):
    """
    Weighted mean size (Acevedo-Trejos et al. 2018)
    It calculates the mean size weighted with total biomass
    ----------
    b = biomass of a specific size class;
    s = the size of that size class (i.e. size array)
    """
    return np.sum(b * s, axis=1) / np.sum(b, axis=1)


def WSizeVar(b, s):
    """
    Weighted size variance (Acevedo-Trejos et al. 2018)
    It calculates size variance weighted with total biomass
    ----------
    b = biomass of a specific size class;
    s = the size of that size class (i.e. size array)
    """
    return (np.sum(b * (s ** 2), axis=1) / np.sum(b, axis=1)) - (WMeanSize(b, s) ** 2)


# Calculate indices
indi = np.zeros((365, 3, 3, 3))  # time array; nutrient regimes; mixing regimes; indices
numSC = np.zeros((3650, 3, 3))

for j in range(3):  # nutrient regimes
    for k in range(3):  # mixing regimes
        indi[:, j, k, 0] = np.sum(P_sc[3285:, :, j, k], axis=1)     # total production
        indi[:, j, k, 1] = WMeanSize(P_sc[3285:, :, j, k], np.logspace(0, 2, 150))
        indi[:, j, k, 2] = WSizeVar(P_sc[3285:, :, j, k], np.logspace(0, 2, 150))
        for t in range(3650):
            numSC[t, j, k] = sum(P_sc[t, 4:, j, k] >= 1e-2)
        #indi[:, j, k, 5] = numSC[3285:, j, k]

# annual average
AnnMean = np.zeros((3, 3, 3))  # nutrient regimes; mixing regimes; indices
AnnMean[:, :, 0] = np.sum(indi[:, :, :, 0], axis=0) / 1000  # annual total biomass
AnnMean[:, :, 1] = np.mean(indi[:, :, :, 1], axis=0)  # annual mean indices
AnnMean[:, :, 2] = np.mean(indi[:, :, :, 2], axis=0)  # annual mean indices

"""
Plotting
"""
Nutclr = ['steelblue', 'gray', 'palegreen']
fig, axs = plt.subplots(1, 2, figsize=(9, 4.3), sharey='all')
plt.tight_layout(pad=6)
plt.subplots_adjust(wspace=0.09, hspace=0.1)
fs, mks = 12, 7
for i in range(2):
    for j in range(3):  # nutrient regimes
        for k in range(3):  # mixing regimes
            axs[i].scatter((AnnMean[j, k, i+1]), AnnMean[j, k, 0], color=Nutclr[j], s=[30, 70, 120][k], marker='o')
            axs[0].set_ylabel('Annual total biomass\n[mmol NL$^{-1}$]', fontsize=fs)
            axs[0].set_xlabel('Annual mean Size\n[µm ESD]', fontsize=fs)
            axs[1].set_xlabel('Annual size Variance\n[µm ESD$^{2}$]', fontsize=fs)
            axs[i].tick_params(labelsize=fs)
            axs[0].axis([0.8, 4.2, -0.2, 5.2])
            axs[1].axis([-0.1, 2, -0.2, 5.2])
            axs[1].set_xticks([0, 0.5, 1, 1.5, 2])
            axs[1].set_xticklabels([0, 0.5, 1, 1.5, 2])
from matplotlib.lines import Line2D
legend = [Line2D([0], [0], marker='o', color=Nutclr[0], ls='none', ms=mks, label='Oligotrophic'),
          Line2D([0], [0], marker='o', color=Nutclr[1], ls='none', ms=mks, label='Eutrophic'),
          Line2D([0], [0], marker='o', color=Nutclr[2], ls='none', ms=mks, label='Hypertrophic'),
          Line2D([0], [0], marker='o', color='black', ls='none', ms=5, label='Constant mixing'),
          Line2D([0], [0], marker='o', color='black', ls='none', ms=8, label='Medium mixing'),
          Line2D([0], [0], marker='o', color='black', ls='none', ms=12, label='High mixing')]
axs[0].legend(handles=legend, fontsize=10, loc=2)
#fig.savefig(fname='/Users/szewing/Desktop/PhD_work/SizeMod/Plotting_report/Fig4.png')
