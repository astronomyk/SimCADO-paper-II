import numpy as np
from matplotlib import pyplot as plt
from astropy.io import ascii


plt.figure(figsize=(10,5))

cat_OB = ascii.read("Melnik_1995-OB_assocs.txt")
rad = np.sqrt(cat_OB["rad_l"] * cat_OB["rad_b"])
log_rad = np.log10(rad)

mu = np.average(log_rad)
sig = np.std(log_rad)
print(10**(mu-sig), 10**(mu), 10**(mu+sig))

plt.hist(rad, bins=np.logspace(0.6,2.2,12), alpha=0.5, label="OB Associations")
plt.plot([10**(mu-sig), 10**(mu-sig)], [0,100], "k--", alpha=0.5)
plt.plot([10**(mu+sig), 10**(mu+sig)], [0,100], "k--", alpha=0.5)
plt.text(10**(mu-sig)*1.05, 20, "-1$\sigma$ = "+str(10**(mu-sig))[:4]+" pc", rotation=90, horizontalalignment="left" )
plt.text(10**(mu+sig)*1.05, 20, "1$\sigma$ = "+str(10**(mu+sig))[:4]+" pc", rotation=90, horizontalalignment="left" )


#------------

cat = ascii.read("Piskunov_2007-OC_masses.txt")
rad = 2*cat["rc"]
log_rad = np.log10(rad)

mu = np.average(log_rad)
sig = np.std(log_rad)
print(10**(mu-sig), 10**(mu), 10**(mu+sig))

plt.hist(rad, bins=np.logspace(-1,2,30), alpha=0.5, color="r", label="Open Clusters", zorder=0)
plt.plot([10**(mu-sig), 10**(mu-sig)], [0,100], "k--", alpha=0.5)
plt.plot([10**(mu+sig), 10**(mu+sig)], [0,100], "k--", alpha=0.5)
plt.text(10**(mu-sig), 33, "-1$\sigma$ = "+str(10**(mu-sig))[:4]+" pc", rotation=90, horizontalalignment="right" )
plt.text(10**(mu+sig), 33, "1$\sigma$ = "+str(10**(mu+sig))[:4]+" pc", rotation=90, horizontalalignment="right" )

#------------

plt.semilogx()
plt.xlim(4E-1, 2E2)
plt.ylim(0, 35)
plt.legend()
plt.ylabel("Number of clusters")
plt.xlabel("Core readius of cluster [pc]")

plt.show()