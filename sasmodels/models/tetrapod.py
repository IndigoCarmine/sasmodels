r"""
Definition
----------

Calculates the scattering from a tetrapod-shaped structure. A tetrapod consists
of four cylindrical arms radiating from a central point, oriented along the
(1,1,1), (-1,-1,1), (-1,1,-1), and (1,-1,-1) directions.

The scattering intensity is calculated as:

.. math::

    I(q) = (\Delta \rho)^2 \left[ 4 I_{\text{arm}}(q) + 6 I_{\text{corr}}(q) \right]

where $\Delta\rho$ is the scattering length density contrast between the
tetrapod and the solvent.

The arm scattering $I_{\text{arm}}(q)$ is obtained by integrating over all
orientations of each cylindrical arm:

.. math::

    I_{\text{arm}}(q) = \int_0^{\pi/4} F_{\text{cylinder}}(q, \theta, L, R)^2 \sin(\theta) d\theta

where $F_{\text{cylinder}}$ is the form factor of a cylinder.

The correlation term $I_{\text{corr}}(q)$ accounts for coherent scattering
between different arms.

Geometry
--------

The four arms are positioned along the following directions:
- Arm 1: (1, 1, 1)
- Arm 2: (-1, -1, 1)
- Arm 3: (-1, 1, -1)
- Arm 4: (1, -1, -1)

Each arm has length $L$ and radius $R$.

References
----------

#. Seoki Kyoo Seo *Korean J. Chem. Eng.* 34(2017) 1192-1198

Authorship and Verification
----------------------------

* **Author:** Yuhei Yamada
* **Last Modified by:**
"""

import numpy as np
from numpy import cos, inf, pi, sin

name = "tetrapod"
title = "Tetrapod with four cylindrical arms"
description = """
    Calculates the scattering from a tetrapod structure with four cylindrical
    arms radiating from a central point.
"""
category = "shape:cylinder"

#             [ "name", "units", default, [lower, upper], "type", "description"],
parameters = [
    ["length", "Ang", 400, [0, inf], "volume", "Cylindrical arm length"],
    ["radius", "Ang", 30, [0, inf], "volume", "Cylindrical arm radius"],
    ["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld", "Tetrapod scattering length density"],
    ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld", "Solvent scattering length density"],
]

source = ["lib/polevl.c", "lib/sas_J1.c", "lib/gauss76.c", "tetrapod.c"]
have_Fq = False
