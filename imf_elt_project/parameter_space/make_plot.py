from operator import attrgetter
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u


class cluster():
    def __init__(self, half_light_radius, distance, mass,
                 mean_mass=0.4 * u.solMass, plate_scale=0.004 * u.arcsec):

        if np.isscalar(half_light_radius):
            half_light_radius *= u.pc
        if np.isscalar(distance):
            distance *= u.pc
        if np.isscalar(mass):
            mass *= u.solMass
        if np.isscalar(mean_mass):
            mean_mass *= u.solMass
        if np.isscalar(plate_scale):
            plate_scale *= u.arcsec

        self.radius = half_light_radius.to(u.pc)
        self.distance = distance.to(u.pc)
        self.mass = mass.to(u.solMass)
        self.mean_mass = mean_mass.to(u.solMass)
        self.plate_scale = plate_scale.to(u.arcsec)

        self.num_stars = self.mass / self.mean_mass

        self.radius_arcsec = (
            np.rad2deg(self.radius / self.distance * u.rad)).to(u.arcsec)
        self.radius_pixel = self.radius_arcsec / self.plate_scale

        self.area = np.pi * self.radius ** 2
        self.area_arcsec = np.pi * self.radius_arcsec ** 2
        self.area_pixel = np.pi * self.radius_pixel ** 2

        self.mean_density = self.num_stars / self.area
        self.mean_density_arcsec = self.num_stars / self.area_arcsec
        self.mean_density_pixel = self.num_stars / self.area_pixel

        self.area_per_star = 1 / self.mean_density
        self.area_per_star_arcsec = 1 / self.mean_density_arcsec
        self.area_per_star_pixel = 1 / self.mean_density_pixel

        self.mean_seperation = 2 * self.radius / np.sqrt(self.num_stars)
        self.mean_seperation_arcsec = 2 * self.radius_arcsec / np.sqrt(
            self.num_stars)
        self.mean_seperation_pixel = 2 * self.radius_pixel / np.sqrt(
            self.num_stars)


class clusters():
    def __init__(self, half_light_radius, distance, mass,
                 mean_mass=0.4 * u.solMass, plate_scale=0.004 * u.arcsec):

        if not isinstance(half_light_radius, (list, tuple, np.ndarray)):
            half_light_radius = [half_light_radius]

        if not isinstance(distance, (list, tuple, np.ndarray)):
            distance = [distance]

        if not isinstance(mass, (list, tuple, np.ndarray)):
            mass = [mass]

        self.radius = half_light_radius
        self.distance = distance
        self.mass = mass

        self.cluster_list = [None] * len(self.radius) * len(
            self.distance) * len(self.mass)

        i = 0
        for rad in self.radius:
            for dist in self.distance:
                for m in self.mass:
                    self.cluster_list[i] = cluster(rad, dist, m,
                                                   mean_mass=mean_mass,
                                                   plate_scale=plate_scale)
                    i += 1

    def __getitem__(self, i):

        if isinstance(i, int):
            print(i)
            return self.cluster_list[i]
        elif isinstance(i, str):
            f = attrgetter(i)
            vals = [f(cluster).value for cluster in self.cluster_list]
            unit = f(self.cluster_list[0]).unit
        return vals * unit


radii = np.array([0.1, 0.14, 0.3, 0.4, 1, 2, 3, 10, 30, 100])
masses = np.array([1E2, 2E2, 1E3, 1.6E3, 5E3, 1E4, 5E4])
distances = np.logspace(3,8, 2)
ylims = np.array([1E-3, 1E3])

c = clusters(radii, distances, masses)


#################################


plt.figure(figsize=(13, 8))
angle = -27

ax1 = plt.gca()

# ------------------ axis 3

ax3 = ax1.twinx()

# -------------- telescopes diffraction limits

rhos = [1E2, 3E2, 1E3, 3E3, 1E4, 3E4, 1E5, 3E5]
ds = [8E3, 20E3, 50E3, 200E3, 800E3, 2E6, 5E6]
for d in ds:
    for rho in rhos:
        plt.scatter(d, rho, c="k", marker="+", s=40)

ys = ylims ** -2
plt.ylim(ys[0], ys[-1])
plt.ylabel("Mean stellar density in core [stars / arcsec$^2$]")
plt.loglog()

# ------------------ axis 2

ax2 = ax1.twiny()
m_M = 5 * np.log10(distances) - 5
plt.xlim(m_M[0], m_M[-1])
plt.xlabel("Distance Modulus [mag]")

mag_lim = 28

abs_mags = [-6.4, -4.5, -3., -0.7, 1.1, 1.76, 3., 3.8, 5.4, 10.4, 12.7]
names = ["B0I", "O5V", "B0V", "B3V", "A0V", "F0V", "G0V", "K0V", "M0V", "M9V",
         "L9V"]
masses = ["", "", "18", "5.5", "2.3", "1.6", "1.1", "0.85", "0.6", "0.065", ""]
clrs = ["#B2EEEC", "#B2EEEC", "#DAE0E0", "#DAE0E0", "#E9ECB8", "#E4EE40",
        "#EECB40", "#EDB279", "#ED8C79", "#DD7C69", "#C3610F", "#C3610F"]
#              O9I        O6          B0         B3      A0         F0        G0       K0         M0        L0        T0        T8
plt.text(mag_lim - 7.75, 950, """
Observational horizon 
for spectral types 
assuming an apparent 
magnitude limit of 
$Ks=28^m$
""",
         horizontalalignment="center", verticalalignment="top", alpha=0.7)

for name, mag, clr, mass in zip(names, abs_mags, clrs, masses):
    plt.plot([mag_lim - mag, mag_lim - mag], [1E-3, 1E3], clr, linewidth=3,
             zorder=1, alpha=0.5)
    s = name + " (" + mass + "$M_{\odot}$)" if mass != "" else name
    plt.text(mag_lim - mag, 0.9E3, s, rotation=270, horizontalalignment="left",
             verticalalignment="top", alpha=0.7)

# ------------------ axis 1
plt.sca(ax1)

# ----------------- YMCs
mass = 1E3
rad_lims = [0.1, 0.4] * u.pc

mask1 = (c["mass"] == mass * u.solMass) * (c["radius"] == rad_lims[0])
mask2 = (c["mass"] == mass * u.solMass) * (c["radius"] == rad_lims[-1])
plt.fill_between(np.array(c["distance"][mask1]),
                 np.array(c["mean_seperation_arcsec"][mask1]),
                 np.array(c["mean_seperation_arcsec"][mask2]),
                 alpha=0.3, color="r", edgecolor="k", zorder=89)
plt.text(5E3,
         0.5 * c["mean_seperation_arcsec"][mask1][0].value,
         "YMC Cores", rotation=angle, verticalalignment="top", zorder=203)

# R136
mask1 = (c["mass"] == 5E3 * u.solMass) * (c["radius"] == rad_lims[0])
plt.plot(np.array(c["distance"][mask1]),
         np.array(c["mean_seperation_arcsec"][mask1]), "k--")
plt.text(3E3,
         0.3 * c["mean_seperation_arcsec"][mask1][0].value,
         "R136 Core", rotation=angle, verticalalignment="top", zorder=210)

# ----------------- Open Cluster
mass = 1E3
rad_lims = [0.1, 3] * u.pc

mask1 = (c["mass"] == mass * u.solMass) * (c["radius"] == rad_lims[0])
mask2 = (c["mass"] == mass * u.solMass) * (c["radius"] == rad_lims[-1])
plt.fill_between(np.array(c["distance"][mask1]),
                 np.array(c["mean_seperation_arcsec"][mask1]),
                 np.array(c["mean_seperation_arcsec"][mask2]),
                 alpha=0.3, color="y", edgecolor="k", zorder=101)
plt.text(5E3,
         1.9 * c["mean_seperation_arcsec"][mask1][0].value,
         str(int(mass)) + " $M_{\odot}$ Open Cluster", rotation=angle,
         verticalalignment="top", zorder=200)
plt.text(c["distance"][mask1][0].value,
         1.1 * c["mean_seperation_arcsec"][mask1][0].value,
         "Core Radius: " + str(rad_lims[0]), rotation=angle, zorder=201)
plt.text(c["distance"][mask2][0].value,
         0.9 * c["mean_seperation_arcsec"][mask2][0].value,
         "Core Radius: " + str(rad_lims[-1]), rotation=angle,
         verticalalignment="top", zorder=202)

# ----------------- OB associations
mass = 5E3
rad_lims = [10, 100] * u.pc

mask1 = (c["mass"] == mass * u.solMass) * (c["radius"] == rad_lims[0])
mask2 = (c["mass"] == mass * u.solMass) * (c["radius"] == rad_lims[-1])
plt.fill_between(np.array(c["distance"][mask1]),
                 np.array(c["mean_seperation_arcsec"][mask1]),
                 np.array(c["mean_seperation_arcsec"][mask2]),
                 alpha=0.3, color="b", edgecolor="k", zorder=100)
plt.text(5E3,
         0.8 * c["mean_seperation_arcsec"][mask1][0].value,
         str(int(mass)) + " $M_{\odot}$ OB Association", rotation=angle,
         verticalalignment="top")
plt.text(c["distance"][mask1][0].value,
         1.1 * c["mean_seperation_arcsec"][mask1][0].value,
         "Core Radius: " + str(rad_lims[0]), rotation=angle)
plt.text(c["distance"][mask2][0].value,
         0.9 * c["mean_seperation_arcsec"][mask2][0].value,
         "Core Radius: " + str(rad_lims[-1]), rotation=angle,
         verticalalignment="top")

# -------------- telescopes diffraction limits

n_fwhm = 2
pix_reses = [0.007, 0.012, 0.045, 0.066, 0.12]
scopes = ["E-ELT/MICADO J", "E-ELT/MICADO Ks", "JWST/NIRCam J", "VLT/NACO Ks",
          "HST/WFC3 F125W"]
textys = [0.004, 0.014, 0.028, 0.066, 0.13]
for pix_res, scope, y in zip(pix_reses, scopes, textys):
    plt.plot(distances, [pix_res * n_fwhm, pix_res * n_fwhm], "k--", alpha=0.5)
    plt.text(0.8 * distances[-1], y * n_fwhm * 1.4,
             scope, horizontalalignment="right", verticalalignment="top",
             alpha=0.7)
plt.text(0.8 * distances[-1], 1E-3 * n_fwhm, """
For a detection, mean 
separations should be 
~""" + str(n_fwhm) + """ PSF FWHMs""",
         horizontalalignment="right", verticalalignment="bottom", alpha=0.7,
         zorder=206)

galaxies = [8.5E3, 50E3, 220E3, 850E3, 2E6, 5E6, 20E6]
names = ["Sgr A* - 8.5 kpc", "LMC - 50 kpc", "Leo I Dwarf - 220 kpc",
         "M33 - 850 kpc", "NGC 300 - 2Mpc", "M83 - 5Mpc",
         "M87 (Virgo) 17Mpc)"]
for gal, name in zip(galaxies[:-1], names[:-1]):
    plt.plot([gal], [1.2E-3], "ko", alpha=0.5)
    plt.text(1.2 * gal, 1.2E-3, name, rotation=90, horizontalalignment="center",
             verticalalignment="bottom", alpha=0.7, zorder=1000)

plt.loglog()
plt.xlim(distances[0], distances[-1])
plt.ylim(ylims[0], ylims[-1])

plt.xlabel("Distance [pc]");
plt.ylabel("Mean separation between stars [arcsec]")

plt.savefig("resolved_stellar_densities.png", format="png")
plt.savefig("resolved_stellar_densities.pdf", format="pdf")
