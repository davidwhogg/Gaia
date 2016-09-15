"""
This file is part of the davidwhogg/Gaia project.
See the LICENSE!

A playground for the new TGAS data.
"""

import numpy as np
from astropy.io import fits
import pylab as plt

def get_data_from_disk():
    for filenum in range(1):
        fn = "../data/TgasSource_000-000-{:03d}.fits".format(filenum)
        hdulist = fits.open(fn)
        hdu = hdulist[1]
        if filenum == 0:
            header = hdu.header
            data = hdu.data
        else:
            data.append(hdu.data)
        hdulist.close()
    return header, data

if __name__ == "__main__":
    tgas_header, tgas = get_data_from_disk()
    I = tgas.field("parallax") > 20.
    plt.clf()
    plt.plot(tgas[I].field("parallax"), tgas[I].field("phot_g_mean_mag"), "k.", alpha=0.5)
    plt.savefig("foo.png")
