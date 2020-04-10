from astropy.io import ascii
from astropy import units as u

import numpy as np
import matplotlib.pyplot as plt
from imf_elt_project.simcado_sandbox import imf


def mu_to_d(mu):
    return 10**((mu+5)/5)


def d_to_mu(d):
    return 5*np.log10(d)-5


sensitivity_limit = 28

lims = [[0.01, 0.08], [0.08, 0.5], [0.5, 300]]
alphas = [0.3, 1.3, 2.3]
scale_factors = [1, 1, 1]
mass = 1E6

masses = imf.imf_population(mass=mass * 8 / 5, lims=lims,
                            alphas=alphas,
                            scale_factors=scale_factors)


# Colours for spectral types
clr = ["#B2EEEC", "#DAE0E0", "#E9ECB8", "#E4EE40", "#EECB40", "#EDB279", "#ED8C79", "#C3610F", "#C3610F"]
#          O          B        A         F         G         K         M        L          T

# Distance of well known objects in kpc
# obj = ["130pc (Pleiades)", "450pc (Orion)", "8.5kpc (Galactic Centre)",
#        "50kpc (LMC)", "770kpc (M31)", "2Mpc (NGC300)", "4Mpc (Cen A)",
#        "5Mpc (M83)", "20Mpc (Virgo)"]
obj = ["450pc (Orion)", "8.5kpc (Galactic Centre)", "50kpc (LMC)",
       "850kpc (M33)", "2Mpc (NGC300)", "5Mpc (M83)", "20Mpc (Virgo)"]
dist = [0.45, 8.5, 50, 850, 2150, 4610, 20000]
text_dist = np.array([0.8, 0.8, 1.12, 0.80, 0.80, 1.12, 1.12])
mu = d_to_mu(dist)

# Generate the IMF boxes
types = ["O", "B", "A", "F", "G", "K", "M", "L", "T"]
mass = [200, 20, 3.4, 1.7, 1.1, 0.8, 0.5, 0.08, 0.03, 0.01]
Mk = [-5.0, -3.2, 0.82, 1.66, 2.81, 3.69, 5.23, 10.41, 12.7, 18.5]
num = np.array([np.sum((masses < mass[i]) * (masses > mass[i + 1]))
                for i in range(len(mass) - 1)])
dm = sensitivity_limit - np.array(Mk)
d_stars = mu_to_d(dm) / 1E3

num = num.astype(float) / np.sum(num)

plt.figure(figsize=(15, 9))

for i in range(len(Mk) - 1):
    x = [d_stars[i], d_stars[i], d_stars[i + 1], d_stars[i + 1]]
    y = [1E-4, num[i], num[i], 1E-4]
    plt.plot(x, y, c=clr[i], linewidth=20)

# Text in graph for stellar type
text_dist_st = mu_to_d(0.5 * (dm[:-1] + dm[1:])) / 1E3
for i in range(len(types)):
    plt.text(text_dist_st[i], 4E-3, types[i], fontsize=24, rotation=0,
             horizontalalignment="center")

# Mass and M_K above the graph
# for i in range(len(types)): plt.text(d_stars[i], 2, str(Mk[i])  ,fontsize=18, rotation=0, horizontalalignment="center")
for i in range(len(mass)):
    plt.text(d_stars[i], 0.95E-3, str(mass[i]) + "M$_\odot$", fontsize=20,
             rotation=-90, verticalalignment="top", horizontalalignment="center")

plt.loglog()
plt.xlim(1E-1, 1E5)
plt.ylim(1E-3, 1.1)

# from matplotlib.ticker import ScalarFormatter
# plt.gca().xaxis.set_major_formatter(ScalarFormatter())
# plt.gca().yaxis.set_major_formatter(ScalarFormatter())
# plt.ticklabel_format(style="plain", axis="both")

plt.xlabel("Distance [kpc]", fontsize=18)
plt.ylabel("Fraction of stars", fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')

for i in range(-1, 6):
    plt.text(10 ** i, 3, f"{5 * (i + 2)} mag",
             horizontalalignment="center", fontsize=18)
    plt.text(100, 5, "Distance Modulus [K band]",
             horizontalalignment="center", fontsize=18)

# Well known objects - vertical dashed line and dot point
for j in range(len(dist)):
    plt.plot((dist[j], dist[j]), (1E-5, 1), "k--", alpha=0.5)
for j in range(len(obj)):
    plt.text(text_dist[j] * dist[j], 0.99, obj[j], verticalalignment="top",
             horizontalalignment="left", fontsize=18, rotation=90)

plt.tight_layout()

plt.savefig("imf_educational.pdf", format="pdf")
plt.savefig("imf_educational.png", format="png")