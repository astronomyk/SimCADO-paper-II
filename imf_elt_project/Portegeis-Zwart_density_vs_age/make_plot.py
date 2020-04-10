import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import ascii
from astropy import units as u

from imf_elt_project.simcado_sandbox import imf


def king_profile_core_fraction(r_core, r_eff):
    """
    Returns the percentage of stars inside r_core for a King profile
    """

    if r_core > r_eff:
        # raise ValueError("re must be greater than rc")
        print("r_core", r_core, ">", r_eff, "r_eff")
        return 0.3

    Io = 1
    r = np.linspace(0, r_eff, int(20 * r_eff / r_core))
    a = np.pi * (r[1:] ** 2 - r[:-1] ** 2)

    I = Io / (1 + (r / r_core) ** 2)
    n = I[1:] * a
    N = np.cumsum(n) / np.sum(n)

    return N[20]


def angle_in_arcseconds(distance, width):
    """
    Returns the angular distance of an object in arcseconds. Units must be consistent
    """

    return np.arctan2(width, distance) * u.rad.to(u.arcsec)


def telescope_diffraction_limit(aperture_size, wavelength, distance=None):
    """
    Returns the diffraction limit of a telescope

    Parameters
    ----------
    aperture_size : float
        [m] The diameter of the primary mirror

    wavelength : float
        [um] The wavelength for diffarction

    distance : float, optional
        Default is None. If ``distance`` is given, the transverse distance for
        the diffraction limit is returned in the same units as ``distance``


    Returns
    -------
    diff_limit : float
        [arcsec] The angular diffraction limit.
        If distance is not None, diff_limit is in the same units as distance

    """

    diff_limit = (((wavelength * u.um) / (aperture_size * u.m)) * u.rad).to(
        u.arcsec).value

    if distance is not None:
        diff_limit *= distance / u.pc.to(u.AU)

    return diff_limit



lmc_oc = ascii.read("LMC_clusters.csv", data_start=1, header_start=0, format="csv")

lmc_oc = lmc_oc[np.invert(lmc_oc["r_core"].mask)]

masses = 10**lmc_oc["logM_phot"]
r_core = lmc_oc["r_core"]
r_eff = lmc_oc["r_eff"]
names = lmc_oc["Name"]
dist = lmc_oc["Distance"] * 1E3
age = lmc_oc["Age"]
mass_lims = lmc_oc["limit_mass"]
gal = lmc_oc["Galaxy"]


n_stars_core = [np.sum(imf.imf_population(mass) > lim) for mass, lim in zip(masses, mass_lims)]
frac = np.array([king_profile_core_fraction(rc, re) for rc, re in zip(r_core, r_eff)])
area = np.pi*angle_in_arcseconds(dist, r_core)**2

rho = 2 * frac * np.array(n_stars_core) / area
ang_eff = angle_in_arcseconds(dist, r_core)




plt.figure(figsize=(9, 6))

clr_dict = {"MW": "#C3610F", "LMC": "#FF8C79", "SMC": "#FF8C79",
            "NGC6822": "#EDB279", "M31": "#FEDB50", "M33": "#EECB40",
            "NGC1569": "#E9ECB8"}

ud = np.unique(gal)
k = 0

for i in [4, 0, 11, 1]:
    mask = (gal == ud[i])
    j = np.where(mask)[0][0]
    if dist[j] > 2.7E6: continue

    k += 1
    axs = []
    sctr = plt.scatter(age[mask], rho[mask], c=clr_dict[ud[i]], label=gal[j],
                       s=25 + 50 * ang_eff[mask], edgecolors='w', zorder=100)  #

    for i in range(len(age)):
        if names[i] in ["Arches", "Westerlund 1", "R136", "NGC330"]:
            x, y = age[i], rho[i]
            if names[i] == "R136":
                x *= 0.95
                y *= 0.85
            plt.text(x, y, "    " + names[i] + "   ", fontsize=8,
                     rotation=45, verticalalignment="bottom",
                     horizontalalignment="left")
n_fwhm = 2.5
wavelength = 1.6
for t, ap in zip(["ELT / MICADO", "JWST / NIRCam", "HST / WFC3-IR"],
                 [38, 6.5, 2.5]):
    y = (n_fwhm * telescope_diffraction_limit(ap, wavelength)) ** -2
    plt.axhline(y, ls=":", c="grey")
    plt.text(1.1E0, 1.1 * y, t, color="grey", verticalalignment="bottom",
             fontsize=14)

# plt.text(1.1E0, 0.9E6, """Limiting magnitude is K=28""", verticalalignment="top")
# Star densities of clusters respect the magnitude limit
# Resolvability limits assume an are FWHM for $\lambda$=1.2$\mu$m
# Data taken from Portegies Zwart (2010)""")

plt.text(1E2, 2E4, "(M31)")
# plt.colorbar()

plt.xlabel("Cluster Age [Myr]", fontsize=10)
plt.ylabel("Observable stellar density [arcsec$^{-2}$]", fontsize=10)

plt.loglog()
plt.xlim(1E0, 2E2)
plt.ylim(1E0, 1E5)

# legend1 = plt.gca().legend(*sctr.legend_elements(),
#                            loc="lower left", title="Classes")
# plt.gca().add_artist(legend1)

#plt.legend(scatterpoints=1, loc=4)

xs, ys = [120]*3, [12, 5, 2],
c = [clr_dict["MW"], clr_dict["LMC"], clr_dict["M31"]]
names = ["MW", "LMC", "M31"]
plt.scatter(xs, ys, c=c, s=100)
for x, y, name in zip(xs, ys, names):
    plt.text(1.1 * x, y, name, verticalalignment="center")


plt.savefig("star_density_vs_age.png", format="png")
plt.savefig("star_density_vs_age.pdf", format="pdf")
