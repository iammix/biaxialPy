from math import pi
import numpy as np
from scipy.spatial import ConvexHull
import geometry
import section_calc as sc


def find_na(x, y, xr, yr, dia, P, Mx, My, fyd, Ec=30 * 10 ** 6, Es=200 * 10 ** 6):
    itr_yn = 0
    max_itr = 1
    tol = 0.1
    # Stiffness ration

    n = Es / Ec
    yn_error = 2 * tol
    yn = (max(y) - min(y)) / 2
    yn_previous = yn

    # Iterate intersection between neutral axis and y-axis until P is satisfactory
    while (abs(yn_error) > tol and itr_yn < max_itr):
        itr_yn += 1
        EtAt = sc.transformed_axial_stiffness(x, y, xr, yr, dia, P, Ec=Ec, Es=Es)

        # Compute distances from concrete vertices and rebars to neutral axis
        dv, dr = sc.compute_dist_to_na(x, y, xr, yr, 0, yn)

        # Get coordinates of most compressed concrete fiber
        # If section is only in tension this will be the most tensioned fibre

        xmax = x[dv.index(min(dv))]
        ymax = y[dv.index(min(dv))]

        # Compute elastic centroid
        xel, yel = sc.elastic_centroid(x, y, xr, yr, dia, Ec=Ec, Es=Es)

        # Compute transformed moment of inertia
        Itx = sc.InertiaX(yel, x, y, xr, yr, dia, Ec=Ec, Es=Es)
        Ity = sc.InertiaY(xel, x, y, xr, yr, dia, Ec=Ec, Es=Es)

        # Compute value for max compression strain within the section strain field
        # In section is only in tension this will be the min tension strain


        eps_max = sc.strain_field_eval(xmax, ymax, P, Mx, My, Ec, EtAt, Itx, Ity)

        # NOTE: test
        eps_max = 0.0035
        # NOTE: Test

        # Compute geometry of the concrete stress block

        # Get list of indices for stress block and extract corresponding x and y coordinates