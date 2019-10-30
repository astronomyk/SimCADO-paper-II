"""
Makes a plot of cumulative star formation vs distance. This SFR is converted
into an "average number of clusters", assuming the 140 M31 number of SF regions
can be equivalent to a SFR of 0.5 Msun / yr.

The number of galaxies that will be seen by MICADO are limited by the latitude
of Armazones: -25 deg +/- 60 deg --> [-85, +35]

"""

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker
import numpy as np
from astropy.io import ascii

TBL = ascii.read("../combined_karachentsev.dat", format="fixed_width")
dec = np.array([float(s+str(d)) for s, d in zip(TBL["DE-"].data, TBL["DEd"].data)])
mask = (TBL["SFRa"] != 0.0) * (dec < 35) * (dec > -85) * (TBL["Dist"] > 0.01)
TBL = TBL[mask]

LIT = ascii.read("../literature_review_cluster_numbers.dat", format="fixed_width")
# mask = LIT["Dist"] > 0.01
# LIT = LIT[mask]

XMIN = 0.02
XMAX = 3
YMIN = 0
YMAX = 3000

log_scale = plt.semilogx


def parsec2distmod(x):
    return 5 * np.log10(x) - 5


def cum_sfr(x):
    return np.sum(10**TBL["SFRa"][TBL["Dist"] < x])


def cum_nclust(x):
    return np.sum(LIT["Num_Clust"][LIT["Dist"] < x])


def plot_literature_number_of_star_clusters(dist, offset=-590):
    cclust = np.array([cum_nclust(d) for d in dist]) + offset
    # plt.scatter(LIT["Dist"], LIT["Num_Clust"])

    plt.plot(dist, cclust, "b", label="Literature")
    log_scale()
    plt.xlabel("Distance [Mpc]")
    plt.ylabel("Cumulative Number of Young Clusters")
    plt.xlim(XMIN, XMAX)
    plt.ylim(YMIN, YMAX)
    plt.text(1.05 * XMIN, 0.95 * cclust[0], "Milky Way clusters", color="b", verticalalignment="top", alpha=0.5)


def plot_cum_sfr_vs_dist(dist):

    csfr = np.array([cum_sfr(d) for d in dist])

    plt.plot(dist, csfr, "g")
    log_scale()
    plt.xlim(XMIN, XMAX)

    plt.xlabel("Distance [Mpc]")
    plt.ylabel("Cumulative Extragalactic Star Formation [$M_\odot/yr$]")


def plot_nclusters_vs_dist(dist, scale_factor=1520, offset=0):
    csfr = np.array([cum_sfr(d) for d in dist])
    nclust = csfr * scale_factor + offset     # Caldwall 2009, 140 young clusters in M31 @ 0.5 Msun/yr

    plt.plot(dist, nclust, "r--", label="$H_{\\alpha}$ Projection")
    log_scale()
    plt.xlim(XMIN, XMAX)
    plt.ylim(YMIN, YMAX)


# def plot_scatter_sfr_vs_dist():
#     plt.scatter(TBL["Dist"], 10**TBL["SFRa"])
#     plt.loglog()
#     plt.xlim(1e-2, 1e2)
#     plt.ylim(1e-7, 1e1)


def plot_spec_type_dropouts():
    top = 2900
    abs_mags = [-6.4, -4.5, -3., -0.7, 1.1, 1.76, 3., 3.8, 5.4, 10.4, 12.7]
    names =  ["B0I",    "O5V",      "B0V",      "B3V",      "A0V",      "F0V",      "G0V",      "K0V",      "M0V",      "M9V",      "L9V"]
    masses = ["",       "",         "18",       "5.5",      "2.3",      "1.6",      "1.1",      "0.85",     "0.6",      "0.065",    ""]
    clrs =  ["#B2EEEC", "#B2EEEC",  "#DAE0E0",  "#DAE0E0",  "#E9ECB8",  "#E4EE40",  "#EECB40",  "#EDB279",  "#ED8C79",  "#DD7C69",  "#C3610F", "#C3610F"]
    mag_lim = 28

    start, end = 4, -1
    for name, mag, clr, mass in zip(names[start:end], abs_mags[start:end], clrs[start:end], masses[start:end]):
        plt.axvline(mag_lim-mag, c=clr, lw=3, alpha=0.5)
        s = name +" ("+mass+"$M_{\odot}$)" if mass != "" else name
        plt.text(mag_lim-mag, top, s, rotation=270, horizontalalignment="left", verticalalignment="top", alpha=0.7)
    plt.xlabel("Distance Modulus [mag]")
    plt.xlim(parsec2distmod(XMIN*1e6), parsec2distmod(XMAX*1e6))

    text = """Observational horizons 
for spectral types 
assuming an apparent 
magnitude limit of 
$Ks=28^m$ """
    # plt.text(mag_lim - 7., top, text, horizontalalignment="center", verticalalignment="top", alpha=0.7)


def plot_notable_galaxy_distances():
    start, end = 0, -1
    for d, name in zip([0.05, 0.847, 2.148, 3.53][start:end], ["LMC", "M33", "NGC300", "M82"][start:end]):
        plt.axvline(d, c="k", ls=":")
        plt.text(1.03*d, 100, name, rotation=90, horizontalalignment="left", verticalalignment="bottom", alpha=0.9)


def plot_imf_curves(scale_factor=0.5, plot_kroupa=False, plot_chabrier=False):
    import imf
    x = np.logspace(np.log10(XMIN), np.log10(XMAX)) / scale_factor
    if plot_chabrier:
        yc = imf.chabrier(x) * x
        plt.loglog(x * scale_factor, yc, c="k", ls="--", alpha=0.3)
        text = "Chabrier" if plot_chabrier else "IMF"
        plt.text(1.03 * XMIN, 1.35 * yc[0], text, rotation=17, alpha=0.5)
    if plot_kroupa:
        yk = imf.kroupa(x) * x
        plt.loglog(x * scale_factor, yk, c="k", ls=":", alpha=0.3)
        text = "Kroupa" if plot_chabrier else "IMF"
        plt.text(1.03 * XMIN, 1.2 * yk[0], text, rotation=14, alpha=0.5)
    plt.ylim(ymax=1)
    plt.gca().get_yaxis().set_visible(False)


################################################################################


def make_plot_num_clusters_estimated_from_halpha_sfr():
    dist = np.logspace(np.log10(XMIN), np.log10(XMAX), 1001)
    plt.figure(figsize=(12, 5))
    ax1 = plt.gca()
    plot_cum_sfr_vs_dist(dist)
    plot_notable_galaxy_distances()
    ax2 = ax1.twinx()
    plot_nclusters_vs_dist(dist)
    ax3 = ax1.twiny()
    plot_spec_type_dropouts()

    plt.show()


def make_plot_catalogued_sfr_vs_halpha_sfr():
    """
    The red and blue lines now only show the total number of clusters that we
    expect to be able to see with MICADO. The blue line is the cumulative sum
    of clusters from the literature. The red line is the projected number of
    clusters that we can infer to exist based on the total Halpha flux of the
    galaxies. The conversion is

    """

    dist = np.logspace(np.log10(XMIN), np.log10(XMAX), 1001)
    plt.figure(figsize=(12, 5))
    ax1 = plt.gca()

    ax3 = ax1.twiny()
    plot_spec_type_dropouts()

    ax2 = ax1.twinx()
    plot_imf_curves(plot_kroupa=True, plot_chabrier=False)

    plt.sca(ax1)
    plot_literature_number_of_star_clusters(dist, offset=0)
    plot_nclusters_vs_dist(dist, offset=590)
    plot_notable_galaxy_distances()
    # plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter())

    plt.legend(loc="best", bbox_to_anchor=(0.3, 0.2))
    plt.tight_layout()

    plt.savefig("young_clusters_within_2Mpc.pdf", format="pdf")
    plt.show()


make_plot_catalogued_sfr_vs_halpha_sfr()
