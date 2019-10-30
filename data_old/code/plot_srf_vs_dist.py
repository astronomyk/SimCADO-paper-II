from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from astropy.io import ascii

TBL = ascii.read("../combined_karachentsev.dat", format="fixed_width")
TBL = TBL[TBL["SFRa"] != 0.0]


def parsec2distmod(x):
    return 5 * np.log10(x) - 5


def cum_sfr(x):
    return np.sum(10**TBL["SFRa"][TBL["Dist"] < x])


def plot_cum_sfr_vs_dist(dist):

    csfr = np.array([cum_sfr(d) for d in dist])

    plt.plot(dist, csfr)
    plt.loglog()
    plt.xlim(0.02, 5)

    plt.xlabel("Distance [Mpc]")
    plt.ylabel("Cumulative Star Formation [$M_\odot/yr$]")


def plot_nclusters_vs_dist(dist):
    csfr = np.array([cum_sfr(d) for d in dist])
    nclust = csfr * (140 / 0.5)     # Caldwall 2009, 140 young clusters in M31 @ 0.5 Msun/yr

    plt.plot(dist, nclust, "b")
    plt.loglog()
    plt.ylabel("Number of Young Clusters")


def plot_scatter_sfr_vs_dist():
    plt.scatter(TBL["Dist"], 10**TBL["SFRa"])
    plt.loglog()
    plt.xlim(1e-2, 1e2)
    plt.ylim(1e-7, 1e1)


def plot_spec_type_dropouts():
    abs_mags = [-6.4, -4.5, -3., -0.7, 1.1, 1.76, 3., 3.8, 5.4, 10.4, 12.7]
    names =  ["B0I",    "O5V",      "B0V",      "B3V",      "A0V",      "F0V",      "G0V",      "K0V",      "M0V",      "M9V",      "L9V"]
    masses = ["",       "",         "18",       "5.5",      "2.3",      "1.6",      "1.1",      "0.85",     "0.6",      "0.065",    ""]
    clrs =  ["#B2EEEC", "#B2EEEC",  "#DAE0E0",  "#DAE0E0",  "#E9ECB8",  "#E4EE40",  "#EECB40",  "#EDB279",  "#ED8C79",  "#DD7C69",  "#C3610F", "#C3610F"]
    mag_lim = 28

    start = 4
    for name, mag, clr, mass in zip(names[start:-1], abs_mags[start:-1], clrs[start:-1], masses[start:-1]):
        plt.axvline(mag_lim-mag, c=clr, lw=3, zorder=1, alpha=0.5)
        s = name +" ("+mass+"$M_{\odot}$)" if mass != "" else name
        plt.text(mag_lim-mag, 60, s, rotation=270, horizontalalignment="left", verticalalignment="top", alpha=0.7)
    plt.xlabel("Distance Modulus [$M_K$]")
    plt.xlim(parsec2distmod(20e3), parsec2distmod(5e6))

    plt.text(mag_lim - 7.75, 70, """    
    Observational horizon 
    for spectral types 
    assuming an apparent 
    magnitude limit of 
    $Ks=28^m$
    """,
    horizontalalignment="center", verticalalignment="top", alpha=0.7)




def plot_notable_galaxy_distances():
    for d, name in zip([0.05, 0.81, 2, 3.5], ["LMC", "M33", "NGC300", "M82"]):
        plt.axvline(d, c="k", ls=":")
        plt.text(1.03*d, 5e-1, name, rotation=90, horizontalalignment="left", verticalalignment="bottom", alpha=0.7)


dist = np.logspace(-2, 1, 31)
plt.figure(figsize=(12, 5))
ax1 = plt.gca()
plot_cum_sfr_vs_dist(dist)
plot_notable_galaxy_distances()
ax2 = ax1.twinx()
plot_nclusters_vs_dist(dist)
ax3 = ax1.twiny()
plot_spec_type_dropouts()


plt.show()