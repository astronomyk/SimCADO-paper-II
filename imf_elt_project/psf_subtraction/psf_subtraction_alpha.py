import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scipy.stats.mstats import theilslopes

from astropy import units as u
from astropy.io import fits
from astropy.table import hstack, Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky

import simcado as sim
from imf_elt_project.simcado_sandbox import imf

dname = ""

# psf_name = "PSF_MAORY_SCAO_Ks_15deg.fits"
# psf_name = "../images_of_clusters/Default_PSF_SCAO.fits"
psf_name = "../images_of_clusters/PSF_MAORY_SCAO_Ks_2.fits"
psf = fits.getdata(psf_name)

bins = 30
# mass = 6000
lims = [[0.01, 0.08], [0.08, 0.5], [0.5, 300]]
alphas = [0.3, 1.3, 2.3]

exptime = 3600
ndit = 1

m_range = [1E3, 1E4]
# m_range = np.logspace(2, 6, 17)
# d = np.array([8.5E3, 20E3, 50E3, 200E3, 800E3, 2E6])
# d = np.logspace(4, 6, 9)
d = np.array([50E3])
dist_mods = 5 * np.log10(d) - 5

alphas_list = [[0.3, 1.3, 2.3]]
scale_factors = [1, 1, 1]

#############################################################################

for dist_mod in dist_mods[:1]:  # GC=14.5, LMC=18.5, Leo I=21.5, M31=24.5
    for mass in m_range:
        for alphas in alphas_list[:1]:
            # try:
            px_fov = 512

            masses = imf.imf_population(mass=mass * 8 / 5, lims=lims,
                                        alphas=alphas,
                                        scale_factors=scale_factors)
            lums = imf.luminosity_from_mass(masses)
            mags = imf.abs_mag_from_mass(masses) + dist_mod

            radius = 4
            n = len(lums)
            a = px_fov
            im = np.zeros((a, a))
            x, y = np.random.randint(radius, a - radius, (2, n))
            xs = 0.004 * (x - a // 2)
            ys = 0.004 * (y - a // 2)

            dist = 10 ** (1 + dist_mod / 5)
            dist = round(dist)
            star_density = round(mass)
            print("star_density", star_density, "dist", dist, "mass", mass)

            cmd = sim.UserCommands()
            cmd["OBS_EXPTIME"] = exptime
            cmd["SCOPE_PSF_FILE"] = psf_name
            cmd["INST_FILTER_TC"] = "Ks"
            cmd["FPA_LINEARITY_CURVE"] = "none"
            cmd["FPA_USE_NOISE"] = "no"
            cmd["FPA_CHIP_LAYOUT"] = "FPA_chip_layout_2as.dat"

            opt = sim.OpticalTrain(cmd)
            fpa = sim.Detector(cmd, small_fov=False)

            src = sim.source.stars(mags=mags, x=xs, y=ys)
            src.apply_optical_train(opt, fpa)

            hdu = fpa.read_out()
            im = hdu[0].data.T
            for kk in range(1, ndit):
                hdu = fpa.read_out()
                im += hdu[0].data.T

            ###################################################################

            results, new_im = imf.iter_psf_photometry(im, psf, radius,
                                                      n_steps=7, e_limit=0.5)

            yr, xr = results["xcentroid"], results["ycentroid"]

            found_srcs = SkyCoord(ra=xr * u.arcsec, dec=yr * u.arcsec)
            init_srcs = SkyCoord(ra=x * u.arcsec, dec=y * u.arcsec)
            idx, d2d, d3d = match_coordinates_sky(found_srcs, init_srcs)
            matched = np.array(
                [x[idx], xr, y[idx], yr, lums[idx], results["m"]]).T
            names = ["x_orig", "x_match", "y_orig", "y_match", "lum_orig",
                     "lum_match"]

            tbl_matched = Table(data=matched, names=names)

            ###################################################################

            n = 26
            bins = np.logspace(-3, 2, n)

            mask = tbl_matched["lum_orig"] > 0.5
            f = theilslopes(tbl_matched["lum_match"][mask],
                            tbl_matched["lum_orig"][mask])[0]

            mass_orig = imf.mass_from_luminosity(tbl_matched["lum_orig"])
            mass_found = imf.mass_from_luminosity(results["m"] / f)
            mass_match = imf.mass_from_luminosity(tbl_matched["lum_match"] / f)
            mass_ratio = mass_match / mass_orig

            results.add_column(Column(data=mass_found, name="mass_found"))

            tbl_mass = Table(data=[mass_orig, mass_match],
                             names=["mass_orig", "mass_match"])
            tbl_matched = hstack(tbl_matched, tbl_mass)

            tbl_stats = imf.binned_clipped_stats(mass_orig, mass_ratio, bins)

            # mask = (tbl_stats["std"] > 0) * (tbl_stats["std"] < 99) * (bins > 1E-2)[1:]
            # if sum(mask) > 0:
            #     mass_trusted = np.interp(0.02, tbl_stats["std"][mask],
            #                              bins[1:][mask])
            # else:
            #     mass_trusted = 1
            #
            for i_row in range(len(tbl_stats[:-4])):
                stds = tbl_stats["std"][i_row:i_row+4]
                if all(stds < 1E-1) and all(stds > 1E-7):
                    mass_trusted = tbl_stats["x_max"][i_row]
                    break

            h_orig, be = np.histogram(mass_orig, bins=bins)
            h_match, be = np.histogram(mass_match, bins=bins)
            h_found, be = np.histogram(mass_found, bins=bins)
            tbl_hist = Table(data=[h_orig, h_match, h_found],
                             names=["n_orig", "n_match", "n_found"])

            tbl_stats = hstack([tbl_stats, tbl_hist])

            fname = dname + "tbl_stats_dist=" + str(dist) + "_rho=" + str(
                star_density) + ".dat"
            # tbl_stats.write(fname, format="ascii")
            # results.write(fname.replace("stats", "found"), format="ascii")

            ###################################################################

            q = radius
            plt.figure(figsize=(12, 12))

            ###################################################################

            plt.axes([0.05, 0.55, 0.4, 0.35])

            vmin_sf = 0.8
            vmax_sf = 10

            plt.title("Original star field")
            vmin = np.median(im) * vmin_sf
            vmax = vmin * vmax_sf
            plt.imshow(im, norm=LogNorm(), interpolation="none",
                       cmap="afmhot", vmin=vmin, vmax=vmax)
            plt.colorbar()

            ###################################################################

            vmin_sf = 0.8
            vmax_sf = 3

            vmin = np.median(new_im) * vmin_sf
            vmax = vmin * vmax_sf
            plt.axes([0.55, 0.55, 0.4, 0.35])
            plt.title("Residuals after fitting and subtraction")
            plt.imshow(new_im, norm=LogNorm(), interpolation="none",
                       cmap="afmhot", vmin=vmin, vmax=vmax)
            plt.colorbar()

            ###################################################################

            plt.axes([0.05, 0.05, 0.4, 0.4])

            plt.hist(masses, bins=bins, alpha=0.5, label="initial", rwidth=1,
                     color="g", edgecolor='black', linewidth=1.2)
            plt.hist(mass_found, bins=bins, alpha=0.5, label="found",
                     rwidth=0.3, color="r", edgecolor='black', linewidth=1.2)
            plt.axvline(0.08, ls=":", c="k")
            plt.axvline(0.5, ls=":", c="k")

            plt.semilogx()
            plt.yscale('log', nonposy='clip')
            plt.xlabel("Mass")
            plt.ylabel("Number of stars in Mass bin")
            plt.xlim(0.01, 100)
            plt.ylim(bottom=1E-1)
            plt.legend(loc=1)

            ###################################################################

            plt.axes([0.55, 0.2, 0.4, 0.25])

            plt.plot(mass_orig, mass_match / mass_orig, ".", alpha=0.2)
            plt.axhline(1, linewidth=3, c="r", zorder=0)
            plt.axvline(mass_trusted, ls=":", c="k")
            plt.loglog()
            plt.ylim(1E-2, 1E2)

            # plt.xlabel("Original Mass")
            plt.ylabel("Ratio of Recovered vs Original Mass ")
            plt.xticks([])
            plt.xlim(0.01, 100)

            ###################################################################

            plt.axes([0.55, 0.05, 0.4, 0.15])

            plt.xlabel("Original Mass")
            plt.yticks([])
            plt.gca().twinx()
            plt.ylabel("Standard Deviation of Recovery")

            plt.plot(bins[1:], tbl_stats["std"], "ko-", alpha=0.3)
            plt.axhline(0.1, ls=":", c="k")
            plt.axvline(mass_trusted, ls=":", c="k")

            plt.loglog()
            plt.xlim(0.01, 100)

            plt.savefig(fname.replace(".dat", ".png"), format="png")
            plt.savefig(fname.replace(".dat", ".pdf"), format="pdf")
