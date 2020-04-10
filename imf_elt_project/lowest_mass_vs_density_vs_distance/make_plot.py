# taken from "plotting_results.ipynb"

from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import units as u
from astropy.table import Table, Column


def get_sky_density(mass, diam, dist):

    density = (2.5 / u.Msun) * mass / diam**2
    arcsec_per_pc = np.arctan(1*u.pc / dist).to(u.arcsec) / u.pc
    sky_density = density / arcsec_per_pc**2

    return sky_density


####################


fnames = glob("./results4/*.dat")


####################


d = np.array(
    [np.array(f[:-4].replace("=", "_").split("_"))[[3, 5]].astype(float) for f
     in fnames])

n = ["dist", "density"]
tbl = Table(data=d, names=n)
# tbl["dist"] = 10**(tbl["dist"]/5+1)

masses = np.arange(len(fnames), dtype=float)

for i in range(len(fnames)):
    data = ascii.read(fnames[i])
    data["std"] = np.nan_to_num(data["std"])

    mask = data["std"] > 0
    x = data["x_center"][mask]
    y = data["std"][mask]

    try:
        j = np.where(y < 0.1)[0][0]
        mass_01 = np.interp(0.1, y[j - 1:j + 1], x[j - 1:j + 1])
    except:
        mass_01 = x[0]
    masses[i] = mass_01
    # print(i, j)

tbl.add_column(Column(data=masses, name="masses"))


####################


clrs = "mrygcbk"[::-1]
xlim = [5E1, 5E5]

plt.figure(figsize=(10, 6))

for i, c, po, in zip(np.unique(tbl["dist"]), clrs, [4, 4, 4, 3, 3, 2, 2]):

    m = (tbl["dist"] == i) * (tbl["density"] < 1E6)

    density = tbl["density"][m]
    masses = tbl["masses"][m]
    j = np.argsort(density)

    xi, yi = density[j], masses[j]

    xf = xi
    if po > 0:
        p = np.polyfit(xi, yi, po)
        yf = np.polyval(p, xf)
    else:
        p = np.polyfit(xi, np.log10(yi), 1)
        yf = 10 ** np.polyval(p, xf)

    plt.errorbar(xi, yi, yerr=0.22 * yi, elinewidth=1, alpha=0.7, fmt=c + ".",
                 label=str(int(np.round(i / 1E3))) + " kpc")
    # plt.plot(xi, yi, c+"o", linewidth=3, alpha=0.3)
    plt.plot(xf, yf, c + "-", linewidth=3, alpha=0.3)
    plt.text(1.1 * xi[0], yf[0] + 0.005, str(int(np.round(i / 1E3))) + " kpc",
             horizontalalignment="left", verticalalignment="bottom",
             fontsize=10, color=c)

    f = np.polyfit(density[j], masses[j], 0)
    x = np.logspace(2, 6, 100)
    y = np.polyval(f, x)
    # plt.plot(x,y, c=c)

plt.gca().set_xscale("log", nonposx='clip')
plt.gca().set_yscale("log", nonposy='clip')

# plt.legend(loc=2, numpoints=1)

plt.xlabel("Star Density [stars/arcsec$^2$]")
plt.ylabel("Limiting observable mass [M$_\odot$]")

abs_mags = [-6.4, -4.5, -3., -0.7, 1.1, 1.76, 3., 3.8, 5.4, 10.4, 12.7]
names = ["B0I", "O5V", "B0V", "B3V", "A0V", "F0V", "G0V", "K0V", "M0V", "M9V",
         "L9V"]
masses = ["", "", "18", "5.5", "2.3", "1.6", "1.1", "0.85", "0.6", "0.065", ""]
clrs = ["#B2EEEC", "#B2EEEC", "#DAE0E0", "#DAE0E0", "#E9ECB8", "#E4EE40",
        "#EECB40", "#EDB279", "#ED8C79", "#DD7C69", "#C3610F", "#C3610F"]
#              O9I        O6          B0         B3      A0         F0        G0       K0         M0        L0        T0        T8

for mass, clr, name in zip(masses[2:-1], clrs[2:-1], names[2:-1]):
    m = float(mass)
    plt.plot([5E5, 1E6], [m, m], c=clr, alpha=0.9, linewidth=3)
    plt.text(5E5, m * 1.1, name)
# plt.grid("on")

plt.axhline(0.08, ls="-.", c="k", alpha=0.5)
plt.axhline(0.5, ls="-.", c="k", alpha=0.5)

plt.text(2E5, 0.5 * 0.9, "0.5 M$_\odot$", verticalalignment="top")
plt.text(2E5, 0.08 * 0.9, "0.08 M$_\odot$", verticalalignment="top")

#################################################################################################################

dist = np.unique(tbl["dist"]) * u.pc

#### Westerlund 1
mass = 63E3 * u.Msun
diam = 1 * u.pc
sds = get_sky_density(mass, diam, dist)

clrs = "mrygcbk"[::-1]
for sd, m, c in zip(sds[:-3], [0.028, 0.11, 0.45, 10.1], clrs):
    # plt.text(sd.value, 1.1*m, "Wd1 =", color=c, rotation=90, horizontalalignment="center", verticalalignment="top")
    plt.plot(sd.value, m, c + "s")

#### ONC
mass = 5E3 * u.Msun
diam = 0.45 * u.pc
sds = get_sky_density(mass, diam, dist)

clrs = "mrygcbk"[::-1]
for sd, m, c in zip(sds[:-3], [0.019, 0.065, 0.24, 3.2], clrs):
    # plt.text(sd.value, m, "= ONC", color=c, rotation=90, horizontalalignment="center", verticalalignment="bottom")
    plt.plot(sd.value, m, c + "o")

#### Ori 1c
mass = 6E3 * u.Msun
diam = 25 * u.pc
sds = get_sky_density(mass, diam, dist)

clrs = "mrygcbk"[:-5][::-1]
for sd, m, c in zip(sds[4:], [1.9, 6.3], clrs):
    # plt.text(sd.value, m, "= Ori 1c", color=c, rotation=90, horizontalalignment="center", verticalalignment="bottom") # weight="bold"
    plt.plot(sd.value, m, c + "^")

# plt.text(6E4, 1.2E-2, "Portegies Zwart et al (2010)")

plt.text(7E5, 3E-2, "Wd1 (YMC): ~60000 M$_\odot$ in ~1pc$^2$     ",
         horizontalalignment="right")
plt.text(7E5, 2E-2, "ONC (OC):   ~5000 M$_\odot$ in ~1pc$^2$     ",
         horizontalalignment="right")
plt.text(7E5, 1.3E-2, "Ori 1c (OB Assoc):   ~6000 M$_\odot$ in ~600pc$^2$",
         horizontalalignment="right")

plt.plot(1.6E4, 3.5E-2, "ks")
plt.plot(1.6E4, 2.2E-2, "ko")
plt.plot(1.6E4, 1.4E-2, "k^")

# plt.text(6E4, 1.2E-2, "Portegies Zwart et al (2010)")

plt.xlim(0.9E2, 1E6)
plt.ylim(0.8E-2, 0.5E2)
plt.gca().set_yticklabels([0, 0, 0.01, 0.1, 1, 10])

plt.savefig("old_trusted_mass.png", format="png")
plt.savefig("old_trusted_mass.pdf", format="pdf")