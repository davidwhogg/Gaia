"""
This file is part of the Gaia project.
Copyright 2016 David W. Hogg (NYU).

# to-do
- add my name!
"""

from matplotlib import rc
rc("font", family="serif", size=12)
rc("text", usetex=True)

import daft

def makepgm(modelcomplexity):
    pgm = daft.PGM([9., 7.], origin=[-0.75, -0.5], observed_style="inner")
    normal = {"alpha": 0.60}
    notyet = {"alpha": 0.15}
    observed = {"facecolor": "0.90", "alpha": 0.70}

    # missions
    instrumentx = 7.
    gaiay = 3.
    keplery = 0.
    apogeey = 1.
    photoy = 2.
    interfy = 4.
    pgm.add_node(daft.Node("apogee", r"\footnotesize spectrographs", instrumentx, apogeey, fixed=True, plot_params=normal))
    pgm.add_node(daft.Node("gaia", r"\footnotesize \textsl{Gaia}", instrumentx, gaiay, fixed=True, plot_params=normal))
    if modelcomplexity > 0:
        pgm.add_node(daft.Node("photo", r"\footnotesize imagers", instrumentx, photoy, fixed=True, plot_params=normal))
        pgm.add_node(daft.Node("kepler", r"\footnotesize \textsl{Kepler}", instrumentx, keplery, fixed=True, plot_params=normal))
        pgm.add_node(daft.Node("interf", r"\footnotesize interferometers", instrumentx, interfy, fixed=True, plot_params=normal))

    # physics
    physicsx = 0
    pgm.add_node(daft.Node("atmos", r"\footnotesize photospheres", physicsx, 2., fixed=True, plot_params=normal))
    if modelcomplexity > 0:
        pgm.add_node(daft.Node("steve", r"\footnotesize stellar evolution", physicsx, 3., fixed=True, plot_params=normal))
        pgm.add_node(daft.Node("struc", r"\footnotesize stellar structure", physicsx, 1., fixed=True, plot_params=normal))
        pgm.add_node(daft.Node("galaxy", r"Milky Way", 3.5, 6., fixed=True, plot_params=notyet, label_params=notyet))

    # latents
    parametery = 4
    pgm.add_node(daft.Node("logg", r"$\log g$", 1.5, parametery, plot_params=normal))
    pgm.add_node(daft.Node("teff", r"$T_{\mathrm{eff}}$", 2.5, parametery, plot_params=normal))
    pgm.add_node(daft.Node("feh", r"$[X/\mathrm{H}]$", 3., parametery + 1, aspect=1.3, plot_params=normal))
    pgm.add_node(daft.Node("distance", r"$D$", 5., parametery + 1., plot_params=normal))
    pgm.add_node(daft.Node("I", r"$\{I_{\lambda}\}$", 2.5, parametery - 2., plot_params=normal))
    if modelcomplexity > 0:
        pgm.add_node(daft.Node("minit", r"$M_{\mathrm{init}}$", 1., parametery + 2, plot_params=notyet, label_params=notyet))
        pgm.add_node(daft.Node("mass", r"$M$", 1., parametery + 1, plot_params=normal))
        pgm.add_node(daft.Node("age", r"age", 2., parametery + 1., plot_params=notyet, label_params=notyet))
        pgm.add_node(daft.Node("dust", r"dust", 4., parametery + 1., plot_params=normal))
        pgm.add_node(daft.Node("L", r"$L_{\mathrm{bol}}$", 3.5, parametery - 1., plot_params=normal))
        pgm.add_node(daft.Node("R", r"$R$", 1., parametery - 1., plot_params=normal))

    # observables
    pgm.add_node(daft.Node("fnorm", r"$\{f^{\mathrm{normed}}_{\lambda}\}$", 3., apogeey, observed=True, aspect=2.1, plot_params=observed))
    pgm.add_node(daft.Node("parallax", r"$\varpi$", 5., gaiay, observed=True, plot_params=observed))
    if modelcomplexity > 0:
        pgm.add_node(daft.Node("seismo", r"$\Delta\nu,\nu_{\mathrm{max}}$", 2., keplery, observed=True, aspect=2.1, plot_params=observed))
        pgm.add_node(daft.Node("m", r"$\{m_j\}$", 4., photoy, observed=True, aspect=2.1, plot_params=observed))
        pgm.add_node(daft.Node("fringes", r"fringes", 6., interfy, observed=True, aspect=2.1, plot_params=observed))

    # pgm
    pgm.add_edge("teff", "I", **normal)
    pgm.add_edge("logg", "I", **normal)
    pgm.add_edge("feh", "I", **normal)
    pgm.add_edge("I", "fnorm", **normal)
    pgm.add_edge("distance", "parallax", **normal)
    if modelcomplexity > 0:
        pgm.add_edge("teff", "L", **normal)
        pgm.add_edge("R", "L", **normal)
        pgm.add_edge("L", "m", **normal)
        pgm.add_edge("dust", "m", **normal)
        pgm.add_edge("distance", "m", **normal)
        pgm.add_edge("logg", "seismo", **normal)
        pgm.add_edge("mass", "seismo", **normal)
        pgm.add_edge("logg", "R", **normal)
        pgm.add_edge("mass", "R", **normal)
        pgm.add_edge("R", "fringes", **normal)
        pgm.add_edge("distance", "fringes", **normal)
        pgm.add_edge("feh", "logg", **normal)
        pgm.add_edge("feh", "teff", **normal)
        pgm.add_edge("minit", "mass", **notyet)
        pgm.add_edge("age", "mass", **notyet)
        for b in ["logg", "teff"]:
            pgm.add_edge("mass", b, **normal)
            pgm.add_edge("age", b, **notyet)

    # physics
    pgm.add_edge("atmos", "I", **normal)
    if modelcomplexity > 0:
        pgm.add_edge("steve", "mass", **notyet)
        pgm.add_edge("steve", "logg", **normal)
        pgm.add_edge("steve", "teff", **normal)
        pgm.add_edge("struc", "seismo", **normal)
        for b in ["minit", "age", "feh", "dust", "distance"]:
            pgm.add_edge("galaxy", b, **notyet)

    # noise models
    pgm.add_edge("apogee", "fnorm", **normal)
    pgm.add_edge("gaia", "parallax", **normal)
    if modelcomplexity > 0:
        pgm.add_edge("photo", "m", **normal)
        pgm.add_edge("kepler", "seismo", **normal)
        pgm.add_edge("interf", "fringes", **normal)

    pgm.render()
    prefix = "parameters{:02d}".format(modelcomplexity)
    pgm.figure.savefig(prefix + ".pdf")
    pgm.figure.savefig(prefix + ".png", dpi=150)

if __name__ == "__main__":
    makepgm(0)
    makepgm(2)
