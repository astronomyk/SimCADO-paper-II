import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import simcado as sim
from project.simcado_sandbox import imf

exptime = 3600
ndit = 1

cmd = sim.UserCommands()
cmd["OBS_EXPTIME"] = exptime
cmd["SCOPE_PSF_FILE"] = "PSF_MAORY_SCAO_Ks_2.fits"
cmd["INST_FILTER_TC"] = "Ks"
cmd["FPA_LINEARITY_CURVE"] = "none"
cmd["FPA_USE_NOISE"] = "no"
cmd["FPA_CHIP_LAYOUT"] = "small"

opt = sim.OpticalTrain(cmd)
fpa = sim.Detector(cmd, small_fov=False)

mass = 1E3
lims = [[0.01, 0.08], [0.08, 0.5], [0.5, 300]]
alphas = [0.3, 1.3, 2.3]
masses = imf.imf_population(mass=mass * 8 / 5, lims=lims, alphas=alphas)
lums = imf.luminosity_from_mass(masses)

border = 4
width = 500
x, y = np.random.randint(border, width - border, (2, len(lums)))
xs = 0.004 * (x - width // 2)  # / (dist / 8.5E3)
ys = 0.004 * (y - width // 2)  # / (dist / 8.5E3)

d = np.array([8.5E3, 20E3, 50E3, 200E3, 800E3])
dist_mods = 5 * np.log10(d) - 5


for dist_mod in dist_mods:  # GC=14.5, LMC=18.5, Leo I=21.5, M31=24.5

    mags = imf.abs_mag_from_mass(masses) + dist_mod
    dist = 10 ** (1 + dist_mod / 5)

    src = sim.source.stars(mags=mags, x=xs, y=ys)
    src.apply_optical_train(opt, fpa)

    fpa.read_out(OBS_EXPTIME=exptime, filename=f"{dist}pc_{mass}Msun.fits")



