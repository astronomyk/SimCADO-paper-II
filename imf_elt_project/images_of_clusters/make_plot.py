import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import scipy.optimize as opt
from scipy.stats import linregress
from scipy.stats.mstats import theilslopes
from scipy.signal import fftconvolve

from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import vstack, hstack, Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky

# from photutils import DAOStarFinder

import simcado as sim
from imf_elt_project.simcado_sandbox import imf


############################


dname = "postage_stamps/"

# psf = fits.getdata("PSF_MAORY_SCAO_Ks_2.fits")
# psf = spi.zoom(psf, 0.5)
# psf = spi.shift(psf, (-0.5,-0.5))

# psf2 = fits.PrimaryHDU(psf)
# psf2.header["CDELT1"] = 0.004
# psf2.header["WAVELENG"] = 2.2
# psf2.writeto("PSF_MAORY_SCAO_Ks_2.fits", clobber=True)


bins = 30
# mass = 6000
lims = [[0.01, 0.08], [0.08, 0.5], [0.5, 300]]
alphas = [0.3, 1.3, 2.3]

exptime = 3600
ndit = 1

# m_range = [1E2, 3E2, 1E3, 3E3, 1E4, 3E4, 1E5, 3E5]
m_range = [1E2, 1E3, 1E4, 1E5]
# m_range = np.logspace(2, 6, 31)
# d = np.array([8.5E3, 20E3, 50E3, 200E3, 800E3, 2E6, 5E6])
d = np.array([8.5E3, 20E3, 50E3, 200E3])
dist_mods = 5 * np.log10(d) - 5

hdus = []
#############################################################################
for mass in m_range:
    for dist_mod in dist_mods:  # GC=14.5, LMC=18.5, Leo I=21.5, M31=24.5

        px_fov = 500
        radius = 4

        masses = imf.imf_population(mass=mass * 8 / 5, lims=lims, alphas=alphas)
        lums = imf.luminosity_from_mass(masses)
        n = len(lums)
        a = px_fov

        x, y = np.random.randint(radius, a - radius, (2, n))

        mags = imf.abs_mag_from_mass(masses) + dist_mod
        dist = 10 ** (1 + dist_mod / 5)

        xs = 0.004 * (x - a // 2)  # / (dist / 8.5E3)
        ys = 0.004 * (y - a // 2)  # / (dist / 8.5E3)

        star_density = mass
        print("star_density", star_density, "dist", dist, "mass", mass)

        chip_layout = """
        #  id    x_cen    y_cen   x_len   y_len   angle    gain
        #        arcsec   arcsec   pixel  pixel     deg  e-/ADU
            0        0        0    """ + str(a) + "  " + str(a) + """   0.0     1.0
        """

        chip_layout = """
        id    x_cen   y_cen   xhw   yhw  x_len   y_len  pixsize  angle    gain
        #        mm      mm    mm    mm  pixel   pixel       mm    deg  e-/ADU
        1      0.00    0.00  7.68  7.68    512     512    0.015    0.0     1.0
        """

        cmd = sim.UserCommands()
        cmd["OBS_EXPTIME"] = exptime
        cmd["SCOPE_PSF_FILE"] = "PSF_MAORY_SCAO_Ks_2.fits"
        cmd["INST_FILTER_TC"] = "Ks"
        cmd["FPA_LINEARITY_CURVE"] = "none"
        cmd["FPA_USE_NOISE"] = "no"
        cmd["FPA_CHIP_LAYOUT"] = "small"

        opt = sim.OpticalTrain(cmd)
        fpa = sim.Detector(cmd, small_fov=False)

        src = sim.source.stars(mags=mags, x=xs, y=ys)
        src.apply_optical_train(opt, fpa)

        hdus += fpa.read_out(OBS_EXPTIME=3600)


##############


plt.figure(figsize=(12, 16))
axs = [plt.axes([x, y, 0.28, 0.22]) for y in [0.1, 0.32, 0.54, 0.76] for x in
       [0.1, 0.38, 0.66]]

j = 0
# for i in [0,1,2, 4,5,6, 8,9,10, 12,13,14]:
for i in [0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11]:
    # for i in [0,1,2,3]:
    plt.sca(axs[j])

    if j == 0:
        fs = 14
        plt.text(-50, 250, "$8.5\,\mathrm{kpc}$", fontsize=fs,
                 verticalalignment="center", rotation=90)
        plt.text(-50, 750, "$20\,\mathrm{kpc}$", fontsize=fs,
                 verticalalignment="center", rotation=90)
        plt.text(-50, 1250, "$50\,\mathrm{kpc}$", fontsize=fs,
                 verticalalignment="center", rotation=90)
        plt.text(-50, 1750, "$200\,\mathrm{kpc}$", fontsize=fs,
                 verticalalignment="center", rotation=90)
        plt.text(-100, 1000, "Distance", fontsize=18,
                 verticalalignment="center", rotation=90)

        plt.text(250, -50, "$100$ $\mathrm{stars}$ $\mathrm{arcsec}^{-2}$",
                 fontsize=fs, horizontalalignment="center")
        plt.text(750, -50, "$1000$ $\mathrm{stars}$ $\mathrm{arcsec}^{-2}$",
                 fontsize=fs, horizontalalignment="center")
        plt.text(1250, -50, "$10\,000$ $\mathrm{stars}$ $\mathrm{arcsec}^{-2}$",
                 fontsize=fs, horizontalalignment="center")
        plt.text(750, -100, "True stellar density", fontsize=18,
                 horizontalalignment="center")

    if j == 8:
        plt.text(620, 0, 'Raw (stacked) pixel counts for 1 hour exposure',
                 fontsize=fs, verticalalignment="center", rotation=90)

    if j == 9:
        plt.plot([50, 300], [460, 460], "w", linewidth=3)
        plt.text(175, 430, '1"', color="w", fontsize=14,
                 horizontalalignment="center")

        plt.plot()

    im = plt.imshow(hdus[i].data[6:-6, 6:-6], interpolation="none",
                    norm=LogNorm(), cmap="hot", vmin=8e4, vmax=1e7,
                    aspect="auto", origin="lower")
    plt.gca().set_xticklabels([])
    plt.gca().set_yticklabels([])

    # plt.text(250,250,str(j)+" "+str(i), color="w")

    j += 1

plt.subplots_adjust(wspace=0, hspace=0)

cax = plt.axes([0.94, 0.1, 0.03, 0.88])
cbar = plt.colorbar(mappable=im, cax=cax)

plt.savefig("clusters_post_stamps_test.png", format="png")