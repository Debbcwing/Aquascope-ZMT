import os
import matplotlib.pyplot as plt
import numpy as np
import netCDF4

################
### Figure 5 ###
################
"""
Fig 5 shows the uptake and grazing forces
for detailed set-up please refer to UptkGraz_10yr.py
"""
Pnum = 150  # specify the nos of phytoplankton size groups
Ynum = 10  # no. of modeling year
nut_rgm = [1, 15, 50]  # µmol NL-1

# Import NetCDF data
os.chdir('wd')
out1 = netCDF4.Dataset('StdRun1.nc')
out15 = netCDF4.Dataset('StdRun15.nc')
out50 = netCDF4.Dataset('StdRun50.nc')

# set up arrays to store extracted data
N = np.zeros((3650, 3, 3))  # model result: nutrient
P_sc = np.zeros((3650, Pnum, 3, 3))  # time; Psize class; nutrient regime; mixing regime
Z_sc = np.zeros((3650, 2, 3, 3))  # time; Zsize class; nutrient regime; mixing regime

# Extract (in loop)
outlist = [out1, out15, out50]
for a in range(len(nut_rgm)):
    print(outlist[a].variables['Nut'][:].shape,
          outlist[a].variables['Zoo'][:].shape,
          outlist[a].variables['Phy'][:].shape)  # (3650, 3) (3650, 2, 3) (3650, 150, 3)
    N[:, a, 0] = outlist[a].variables['Nut'][:, 0]
    N[:, a, 1] = outlist[a].variables['Nut'][:, 2]
    N[:, a, 2] = outlist[a].variables['Nut'][:, 1]
    Z_sc[:, :, a, 0] = outlist[a].variables['Zoo'][:, :, 0]
    Z_sc[:, :, a, 1] = outlist[a].variables['Zoo'][:, :, 2]
    Z_sc[:, :, a, 2] = outlist[a].variables['Zoo'][:, :, 1]
    P_sc[:, :, a, 0] = outlist[a].variables['Phy'][:, :, 0]
    P_sc[:, :, a, 1] = outlist[a].variables['Phy'][:, :, 2]
    P_sc[:, :, a, 2] = outlist[a].variables['Phy'][:, :, 1]  # simlist str = (time; concentration; MixRgm)


def Vol2ESD(v):
    """
    Vol2ESD(1) = 1.2407
    -----------
    This function is used to convert the parameter values taken from literatures that
    uses biovolume as their plankton size unit
    """
    esd = ((v * 6) / np.pi) ** (1 / 3)
    return esd


def Allo_mumax(PhySize):
    """
    Allometric max. uptake rate function
    mumax constraints the max. nutrient uptake of phytoplankton
    """
    mu_alpha = -0.36
    mu_beta = 10 ** 0.69
    mumax = mu_beta * ((PhySize / Vol2ESD(1)) ** mu_alpha)  # allometric mumax
    return mumax

def AlloUptk_normalized(nut, PhySize, Pi):
    """
    To normalize uptake against phytoplantkon biomass on each size class in order to visualize if the model
    run according to hypothesized theories
    """
    Kn_alpha = 0.52
    Kn_beta = 10 ** -0.71
    Kn = Kn_beta * ((PhySize / Vol2ESD(1)) ** Kn_alpha)
    Kn_MM = nut / (nut + Kn)    # Michelis-Menton formulation on nutrient saturation
    uptake = (Allo_mumax(PhySize) * Kn_MM) * Pi
    return uptake / (Pi*Allo_mumax(PhySize))

def AlloGraz_normalized(PhySize, ZooSize, Pi, SizeTolerance, Zi):
    """
    PhySize = Size class/array of phytoplankton population;
    ZooSize = Size class/array of zooplankton population;
    Pi = Phytoplankton biomass
    """
    Kp = 3  # Phytoplankton half-saturation constant for zooplankton  [µmol N/L]
    imax_alpha = -0.4  # Hansen et al. 1994, Banas, 2011
    imax_beta = 26  # Hansen et al. 1994, Banas, 2011 [day^-1]
    Imax = imax_beta * ((ZooSize / 1) ** imax_alpha)
    OPS = 0.65 * ((ZooSize / 1) ** 0.56)  # µm ESD
    gpP = np.exp(-((np.log10(PhySize) - np.log10(OPS)) / SizeTolerance) ** 2)
    denom = np.sum(gpP * Pi)
    graz = Imax * ((gpP * Pi) / (Kp + denom)) * Zi
    return graz / (Pi*Allo_mumax(PhySize))

# set up arrays to store calculated data
NU_norm = np.zeros((3650, Pnum, 3, 3))  # NU metric (i.e. 0-1)
GrazS_norm = np.zeros((3650, Pnum, 3, 3))   # grazing metric (i.e. 0-1)
GrazL_norm = np.zeros((3650, Pnum, 3, 3))

    # Size arrays
Psize = np.logspace(0, 2, 150)  # define phytoplankton size array
Zsize = [5, 200]

#### nutrient uptake ----
for i in range(3):  # nutrient regime
    for j in range(3):  # mixing regime
        for k in range(150):  # Pnum
            NU_norm[:, k, i, j] = AlloUptk_normalized(N[:, i, j], Psize[k], P_sc[:, k, i, j])
            GrazS_norm[:, k, i, j] = AlloGraz_normalized(Psize[k], 5, P_sc[:, k, i, j], 0.2, Z_sc[:, 0, i, j])
            GrazL_norm[:, k, i, j] = AlloGraz_normalized(Psize[k], 200, P_sc[:, k, i, j], 0.5, Z_sc[:, 1, i, j])

# Total and average nutrient uptake of Psize in 10 yrs
NU_norm_mean10yr = np.mean(NU_norm, axis=0)
NU_norm_meanlastyr = np.mean(NU_norm[3285:, :, :, :], axis=0)   # steady-state solution

# Total and average grazing of Psize in 10 yrs
GrazS_norm_sum10yr = np.sum(GrazS_norm, axis=0)
GrazL_norm_sum10yr = np.sum(GrazL_norm, axis=0)
GrazS_norm_mean10yr = np.mean(GrazS_norm, axis=0)
GrazL_norm_mean10yr = np.mean(GrazL_norm, axis=0)
GrazS_norm_meanlastyr = np.mean(GrazS_norm[3285:, :, :, :], axis=0)
GrazL_norm_meanlastyr = np.mean(GrazL_norm[3285:, :, :, :], axis=0)

"""
Plotting
"""
# Plot settings
Nutclr = ['steelblue', 'gray', 'palegreen']
Nutlab = ['Oligotrophic', 'Eutrophic', 'Hypertrophic']
Mixlist = ['Constant', 'Medium mixing', 'High mixing']
lw, fs = 1.5, 11

## Figure 6
fig, axs = plt.subplots(3, 3, figsize=(8, 6), sharex='all', sharey='all')
plt.subplots_adjust(hspace=0.1, wspace=0.15)
for i in range(3):
    for j in range(3):  # mixing
        #axs[j, i].plot(Psize, NU_norm_mean10yr[:, i, j], color='#52b9c4', label='Nutrient uptake')
        axs[j, i].bar(Psize, NU_norm_mean10yr[:, i, j], color='#52b9c4', label='Nutrient uptake')
        #axs[j, i].plot(Psize, GrazS_norm_mean10yr[:, i, j], color='gray', label='Small grazer')
        axs[j, i].bar(Psize, GrazS_norm_mean10yr[:, i, j], color='gray', label='Small grazer')
        #axs[j, i].plot(Psize, GrazL_norm_mean10yr[:, i, j], color='orange', label='Large grazer')
        axs[j, i].bar(Psize, GrazL_norm_mean10yr[:, i, j], color='orange', label='Large grazer')
        axs[j, i].set_xscale('log')
        axs[j, i].axis([1e0, 1e2, 0.001, 0.5])
        axs[0, i].set_title(Nutlab[i], fontweight='bold', fontsize=10)
        axs[j, i].set_xticks([1e0, 1e1, 1e2])
        axs[j, i].set_xticklabels([1, 10, 100])
        axs[j, i].set_yticks([0.01, 0.25, 0.5])
        axs[j, 2].yaxis.set_label_position('right')
        axs[j, 2].yaxis.set_ticks_position('right')
        axs[j, i].yaxis.set_ticks_position('both')
        axs[j, 2].set_ylabel(Mixlist[j], fontweight='bold', fontsize=10, rotation=270, labelpad=13)
axs[2, 1].set_xlabel('Phytoplankton size, $S_{i}$ [µm ESD]', fontsize=fs+2)
axs[1, 0].set_ylabel('Normalized nutrient uptake and grazing', fontsize=fs)
from matplotlib.lines import Line2D
legend = [Line2D([0], [0], color='#52b9c4', lw=5, label='Nutrient uptake'),
          Line2D([0], [0], color='gray', lw=5, label='Grazing - Z1'),
          Line2D([0], [0], color='orange', lw=5, label='Grazing - Z2')]
axs[0, 0].legend(handles=legend, loc=1, handlelength=1.0, fontsize=10, ncol=3, bbox_to_anchor=(3.0, 1.5))
#fig.savefig(fname='/Users/szewing/Desktop/PhD_work/SizeMod/Plotting_report/Fig5.png')

