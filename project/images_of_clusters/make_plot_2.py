from copy import deepcopy
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import simcado as sim


def make_cluster_hdus(mass=2E4, half_light_radius=1, exptime=3600, show=False):

    cmd = sim.UserCommands()
    cmd["OBS_EXPTIME"] = exptime
    cmd["SCOPE_PSF_FILE"] = "PSF_MAORY_SCAO_Ks_2.fits"
    cmd["INST_FILTER_TC"] = "Ks"
    cmd["FPA_LINEARITY_CURVE"] = "none"
    cmd["FPA_USE_NOISE"] = "no"
    cmd["FPA_CHIP_LAYOUT"] = "FPA_chip_layout_centre.dat"

    opt = sim.OpticalTrain(cmd)
    fpa = sim.Detector(cmd, small_fov=False)

    dist_orig = 1E3
    cluster = sim.source.cluster(mass=mass, distance=dist_orig, half_light_radius=half_light_radius)

    x_orig = deepcopy(cluster.x)
    y_orig = deepcopy(cluster.y)
    f_orig = deepcopy(cluster.weight)

    d = np.array([8.5E3, 20E3, 50E3, 200E3, 800E3])
    dist_mods = 5 * np.log10(d) - 5

    plt.figure(figsize=(15, 5))

    for ii, dist_mod in enumerate(dist_mods):  # GC=14.5, LMC=18.5, Leo I=21.5, M31=24.5

        dist = 10 ** (1 + dist_mod / 5)
        cluster.x = x_orig * dist_orig / dist
        cluster.y = y_orig * dist_orig / dist
        cluster.weight = f_orig * (dist_orig / dist)**2

        cluster.apply_optical_train(opt, fpa)
        hdu = fpa.read_out(OBS_EXPTIME=exptime, filename=f"{int(dist)}kpc_{int(mass)}Msun.fits")

        if show:
            plt.subplot(1, 5, ii+1)
            plt.imshow(hdu[0].data, norm=LogNorm(), vmin=1.61E5, vmax=1.8E5)

    if show:
        plt.show()


def plot_cluster_hdus():
    from glob import glob
    from astropy.io import fits

    dists = [8.5, 20, 50, 200, 800]

    fnames = glob("*Msun.fits")
    jj = np.argsort([int(f.split("kpc")[0]) for f in fnames])

    plt.figure(figsize=(15, 6.03))
    for ii, fname in enumerate(np.array(fnames)[jj]):
        cmap = "hot"
        im = fits.getdata(fname)
        vmin = np.min(im)
        vmax = 100 * np.median(im)

        cax = plt.axes([0 + 0.2*ii, 0.5, 0.2, 0.5])
        # plt.subplot(2, 5, 1 + ii)
        cax = plt.gca()
        cax.imshow(im, norm=LogNorm(), vmin=vmin, vmax=vmax,
                   cmap=cmap, origin="lower")

        cax.text(200, 3700, f"{dists[ii]} kpc", color="w",
                 horizontalalignment="left", fontsize=14)
        cax.axis("off")

        cax.plot([1, 1798, 1798, 2298, 2298, 1798, 2298, 4095],
                 [1, 2298, 1798, 1798, 2298, 2298, 2298, 1], c="w", lw=1)

        if ii == 4:
            cax.plot([3840, 3840], [2048-1250, 2048+1250], c="w", lw=3)
            cax.text(3760, 2048, '10"', fontsize=14, color="w",
                     horizontalalignment="right", verticalalignment="center")

        # plt.subplot(2, 5, 6 + ii)
        # cax = plt.gca()
        cax = plt.axes([0 + 0.2*ii, 0.0, 0.2, 0.5])

        m, d = 2048, 250
        cax.imshow(im[m - d:m + d, m - d:m + d], norm=LogNorm(), vmin=vmin,
                   vmax=vmax, cmap=cmap, origin="lower")

        cax.axis("off")
        if ii == 4:
            cax.plot([480, 480], [125, 375], c="w", lw=3)
            cax.text(470, 250, '1"', fontsize=14, color="w",
                     horizontalalignment="right", verticalalignment="center")

    # plt.subplots_adjust(wspace=0., hspace=0.)

make_cluster_hdus(exptime=60)
plot_cluster_hdus()

plt.savefig("5_clusters.png", format="png")
plt.savefig("5_clusters.pdf", format="pdf")

plt.show()