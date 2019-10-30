from os import path as pth
import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column

PKG_DIR = pth.abspath(pth.dirname(__file__))


def combine_karachantsev_tables(write_to_disk=False):

    dname = pth.join(PKG_DIR, "../", "karachantsev_sfr_catalogue/")
    sfr = ascii.read(dname+"table3.dat", readme=dname+"ReadMe.txt")
    dname = pth.join(PKG_DIR, "../", "karachantsev_list_of_galaxies/")
    phot = ascii.read(dname+"table2.dat", readme=dname+"readme.txt")
    dist = ascii.read(dname+"table6.dat", readme=dname+"readme.txt")

    wanted_col_names = ["Name", "RAh", "RAm", "RAs", "DE-", "DEd", "DEm", "DEs", "TT", "SFRa"]
    dtypes =           ["<U17", "<i4", "<i4", "<f8", "<U1", "<i4", "<i4", "<i4", "<i4", "<f8"]
    tbl = Table(names=wanted_col_names + ["KLum"] + ["DM"], dtype=dtypes+["<f8", "<f8"])

    for row in dist:
        if row["Name"] in sfr["Name"] and row["Name"] in phot["Name"]:
            sfr_ii = np.where(sfr["Name"] == row["Name"])[0][0]
            phot_ii = np.where(phot["Name"] == row["Name"])[0][0]
            vals = [sfr[sfr_ii][col] for col in wanted_col_names] + [phot[phot_ii]["KLum"]] + [row["DM"]]
            tbl.add_row(vals)

    dist_col = Column(name="Dist", data=np.round(10**(0.2*tbl["DM"]+1) * 1e-6, 3))
    sfr_col = Column(name="SFR",  data=np.round(10**(tbl["SFRa"]+tbl["KLum"]), 3))

    tbl.add_columns([dist_col, sfr_col])

    if write_to_disk:
        tbl.write(pth.join(PKG_DIR, "../", "combined_karachentsev.dat"),
                  format="ascii.fixed_width", overwrite=True)

    return tbl
