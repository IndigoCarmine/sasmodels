r"""
Definition
----------

This model computes the scattering from a torus with an elliptical tube
cross-section and a concentric shell.

The geometric parameters are:

* $R$: major (ring) radius of the torus centerline
* $a$: core minor radius in the radial direction
* $\nu$: aspect ratio of the elliptical cross-section
* $t$: shell thickness

.. figure:: img/torus_elliptical_shell_geometry.png

    Schematic geometry of the torus with elliptical tube cross-section.

The same aspect ratio $\nu$ is used for both the core and outer cross-sections,
so the outer semi-axes are $a+t$ and $\nu(a+t)$.

For a given orientation angle $\theta$ between the torus symmetry axis and
$\vec q$, the kernel evaluates the amplitude using a numerical integral over
the tube section:

.. math::

    F(q,\theta; x, \Delta\rho)
    = 4\pi\,\Delta\rho\int_{R-x}^{R+x}
      r\,J_0\!\left(qr\sin\theta\right)
      \frac{\sin\!\left(q\,\gamma\cos\theta\right)}{q\cos\theta}\,dr

with

.. math::

    \gamma = \nu\sqrt{x^2-(r-R)^2}

The core-shell amplitude is formed from outer and inner contributions:

.. math::

    F_{cs}(q,\theta) =
    F\!\left(q,\theta;a + t,\rho_{shell}-\rho_{solvent}\right)
    -F\!\left(q,\theta;a,\rho_{shell}-\rho_{core}\right)

and the orientationally averaged intensity is

.. math::

    I(q) = S(q)\int_0^{\pi/2}\!\left|F_{cs}(q,\theta)\right|^2\sin\theta\,d\theta

where

.. math::

    S(q)=1+\frac{\kappa}{1+(q\zeta)^2}

Here $S(q)$ represents the structure factor accounting for
interparticle correlations.

References
----------
#.  T. Kawaguchi, *J. Appl. Crystallogr*, 34(2001) 580-584
#.  S. Förster, *J. Phys. Chem.*, 103(1999) 6657-6668
#.  M. J. Hollamby, *Angew. Chem. Int. Ed.* 55(2016) 9890

Authorship and Verification
---------------------------

* **Author:** Itsuki Tajima and Yuhei Yamada (Github user name: Indigo Carmine, https://orcid.org/0009-0003-9780-4135)
* **Last Modified by:**
* **Last Reviewed by:**
"""

from numpy import inf

# ---- SasView model metadata -------------------------------------------------
name = "torus_elliptical_shell"
title = "core-shell torus with elliptical cross-section"
description = "Core-shell torus with elliptical tube cross-section"
category = "shape:cylinder"
parameters = [
    # name          units        default  [min,   max]  type     description
    ["radius", "Ang", 100.0, [0, inf], "volume", "Torus major radius R"],
    [
        "core_radius",
        "Ang",
        5.0,
        [0, inf],
        "volume",
        "Elliptical core minor radius a",
    ],
    ["thickness", "Ang", 2.0, [0, inf], "volume", "Shell thickness"],
    ["nu", "", 1.0, [0.1, 10.0], "volume", "Aspect ratio b/a"],
    ["sld_core", "1e-6/Ang^2", 0.0, [-inf, inf], "sld", "Core SLD"],
    ["sld_shell", "1e-6/Ang^2", 1.0, [-inf, inf], "sld", "Shell SLD"],
    ["sld_solvent", "1e-6/Ang^2", 0.0, [-inf, inf], "sld", "Solvent SLD"],
]

valid = "radius >= core_radius + thickness"

# -- tell sasmodels that a C kernel is provided -------------------------------
source = [
    "lib/polevl.c",
    "lib/sas_j0.c",
    "lib/gauss76.c",
    "torus_elliptical_shell.c",
]  # compiled together with default libs
