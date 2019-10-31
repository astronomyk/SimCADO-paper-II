from os import path as pth
import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column, vstack

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

SUP_DATA_DIR = pth.abspath(pth.join(pth.dirname(__file__), "../"))


def combine_karachantsev_tables(write_to_disk=False):

    dname = pth.join(SUP_DATA_DIR, "karachantsev_sfr_catalogue/")
    sfr = ascii.read(dname+"table3.dat", readme=dname+"ReadMe.txt")
    dname = pth.join(SUP_DATA_DIR, "karachantsev_list_of_galaxies/")
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
        tbl.write(pth.join(SUP_DATA_DIR, "combined_karachentsev.dat"),
                  format="ascii.fixed_width", overwrite=True)

    return tbl


def read_combined_table():
    # "SFRu" is in log10(SFR/L_K), "Dist" in Mpc
    tbl = ascii.read(pth.join(SUP_DATA_DIR, "combined_karachentsev.dat"),
                     format="fixed_width")
    return tbl


def put_hera_table_in_karachantsev_format(write=False, mpc_or_kpc="kpc"):
    scale_factor = 1e-3 if "k" in mpc_or_kpc else 1e-6

    hera = ascii.read(pth.join(SUP_DATA_DIR,
                               "hearsarc_milky_way_cluster",
                               "hearsarc_milky_way_clusters_list.txt"),
                      format="fixed_width")
    mwoc = hera[hera["log_age"] < 8]
    mwoc_dict = {"Name": mwoc["name"].data,
                 "RAh": [int(val[0:2]) for val in mwoc["ra"]],
                 "RAm": [int(val[3:5]) for val in mwoc["ra"]],
                 "RAs": [float(val[6:8]) for val in mwoc["ra"]],
                 "DE-": [val[:1] for val in mwoc["dec"]],
                 "DEd": [int(val[1:3]) for val in mwoc["dec"]],
                 "DEm": [int(val[4:6]) for val in mwoc["dec"]],
                 "DEs": [int(val[7:8]) for val in mwoc["dec"]],
                 "TT": [99] * len(mwoc),
                 "SFRa": [-99.] * len(mwoc),
                 "KLum": [-99.] * len(mwoc),
                 "DM": 5 * np.log10(mwoc["distance"]) - 5,
                 "Dist": mwoc["distance"] * scale_factor,
                 "n_clust": [1] * len(mwoc),
                 "reference": ["heasarc_mwsc"] * len(mwoc),
                 }

    mwoc_tbl = Table(data=[mwoc_dict[key] for key in mwoc_dict],
                     names=mwoc_dict.keys())

    if write:
        mwoc_tbl.write(pth.join(SUP_DATA_DIR, "mwocs_in_galocs_format.dat"),
                       format="ascii.fixed_width", overwrite=True)

    return mwoc_tbl


def combine_hera_with_karachantsev(write=False, mpc_or_kpc="kpc"):
    """
    - The hera data is only for Dec<+35. Distances are in pc.
    - The karachantsev data is for all galaxies WITHOUT the number of cluster
      and references. These were added later from the table
      ``literature_review_cluster_numbers.dat``
    - Karachantsev distances in Mpc

    """

    scale_factor = 1e3 if "k" in mpc_or_kpc else 1

    mwocs = put_hera_table_in_karachantsev_format(False, mpc_or_kpc)
    galocs = read_combined_table()
    galocs["Dist"] = galocs["Dist"] * scale_factor

    print("", mwocs.dtype, "\n", galocs.dtype)

    comb_tbl = vstack([mwocs, galocs])

    if write:
        comb_tbl.write(pth.join(SUP_DATA_DIR, "clusters_plus_galaxies.dat"),
                       format="ascii.fixed_width", overwrite=True)

    return comb_tbl

# combine_karachantsev_tables(write=False)
# combine_hera_with_karachantsev(write=False)




