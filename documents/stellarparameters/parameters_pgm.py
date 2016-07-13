"""
This file is part of the Gaia project.
Copyright 2016 David W. Hogg (NYU).
"""

from matplotlib import rc
rc("font", family="serif", size=12)
rc("text", usetex=True)

import daft

pgm = daft.PGM([7, 5], origin=[-1., -0.5], observed_style="inner")

# physics
pgm.add_node(daft.Node("galaxy", r"\footnotesize Galaxy", 3.5, 4, fixed=True))
pgm.add_node(daft.Node("imf", r"\footnotesize IMF", 0.5, 4, fixed=True))
pgm.add_node(daft.Node("steve", r"\footnotesize stellar evolution", 0, 2, fixed=True))
pgm.add_node(daft.Node("atmos", r"\footnotesize atmospheres", 0, 1, fixed=True))

# missions
pgm.add_node(daft.Node("gaia", r"\footnotesize \textsl{Gaia} noise model", 4.5, 2, fixed=True))
pgm.add_node(daft.Node("photo", r"\footnotesize photometry noise model", 4.5, 1, fixed=True))
pgm.add_node(daft.Node("apogee", r"\footnotesize \textsl{APOGEE} instrument model", 4.5, 0, fixed=True))

# latents
pgm.add_node(daft.Node("mass", r"$M$", 0.5, 3))
pgm.add_node(daft.Node("age", r"age", 1.5, 3))
pgm.add_node(daft.Node("feh", r"$[X/\mathrm{H}]$", 2.5, 3, aspect=1.3))
dusty = {"alpha": 0.25}
pgm.add_node(daft.Node("extinction", r"dust", 2.5, 2, plot_params=dusty, label_params=dusty))
pgm.add_node(daft.Node("distance", r"$D$", 3.5, 3))
pgm.add_node(daft.Node("mrl", r"$\begin{array}{c}M,R,L_{\mathrm{bol}} \\ {\scriptstyle (\log g,\,T_{\mathrm{eff}})}\end{array}$", 1.5, 2, scale=1.6, aspect=1.5))
pgm.add_node(daft.Node("Ltrue", r"$\{L_{\lambda}\}$", 1.5, 1))

# observables
pgm.add_node(daft.Node("fnorm", r"$\{f^{\mathrm{normed}}_{\lambda}\}$", 1.5, 0, observed=True, aspect=2.1))
pgm.add_node(daft.Node("m", r"$\{m_{\mathrm{band}}\}$", 2.5, 1, observed=True, aspect=1.8))
pgm.add_node(daft.Node("parallax", r"$\varpi$", 3.5, 2, observed=True))

# pgm
pgm.add_edge("mass", "mrl")
pgm.add_edge("age", "mrl")
pgm.add_edge("feh", "mrl")
pgm.add_edge("mrl", "Ltrue")
pgm.add_edge("Ltrue", "fnorm")
pgm.add_edge("Ltrue", "m")
pgm.add_edge("extinction", "m", **dusty)
pgm.add_edge("distance", "m")
pgm.add_edge("distance", "parallax")

# physics
pgm.add_edge("galaxy", "age")
pgm.add_edge("galaxy", "feh")
pgm.add_edge("galaxy", "extinction", **dusty)
pgm.add_edge("galaxy", "distance")
pgm.add_edge("imf", "mass")
pgm.add_edge("steve", "mrl")
pgm.add_edge("atmos", "Ltrue")

# noise models
pgm.add_edge("gaia", "parallax")
pgm.add_edge("photo", "m")
pgm.add_edge("apogee", "fnorm")

pgm.render()
pgm.figure.savefig("parameters.pdf")
pgm.figure.savefig("parameters.png", dpi=150)
