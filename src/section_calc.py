from math import pi, cos, sin, tan, atan, atan2, sqrt, ceil, floor
import numpy as np
import matplotlib.path as pltPath
import geometry

# CONSTANTS
EC = 30 * 10 ** 6
ES = 200 * 10 ** 6


def rebars_in_stress_block(x_sb, y_sb, xr, yr) -> list:
    """
    :return: Return a list with entry 'True' for rebars located inside the stress block, 'False', otherwise
    """
    if xr and yr:
        rebar_coords = [[xr[i], yr[i]] for i in range(len(xr))]
    else:
        raise ValueError('No rebars in section')

    # Compute are of stress block
    Asb = geometry.polygon_area(x_sb, y_sb)

    if Asb != 0:
        # Arrange stress block coordinates
        sb_poly = [[x_sb[i], y_sb[i]] for i in range(len(x_sb))]

        # Check if rebars are inside the strees block
        path = pltPath.Path(sb_poly)
        # Return 'True' if rebar is inside stress block
        rebars_inside = path.contains_points(rebar_coords)
    else:
        # All rebars are in tension (all entries are 'False')
        rebars_inside = [False] * len(xr)

    return rebars_inside


def compression_tension_rebars(x, y, xr, yr, dia):
    """
    :return: Return lists of rebar coordinates and diameters for rebars in compression and tension, respectively.
    """
    # Evaluate if rebars are inside or outside stress block (return list with 'True' or 'False')
    rebar_eval = rebars_in_stress_block(x, y, xr, yr)

    # Number of rebars
    nb = len(dia)

    # Extract rebars in compression
    dia_comp = [dia[i] for i in range(nb) if rebar_eval[i]]
    xr_comp = [xr[i] for i in range(nb) if rebar_eval[i]]
    yr_comp = [yr[i] for i in range(nb) if rebar_eval[i]]

    # Extract rebars in tension
    dia_tens = [dia[i] for i in range(nb) if not rebar_eval[i]]
    xr_tens = [xr[i] for i in range(nb) if not rebar_eval[i]]
    yr_tens = [yr[i] for i in range(nb) if not rebar_eval[i]]

    return xr_comp, yr_comp, dia_comp, xr_tens, yr_tens, dia_tens


def InertiaX(yc, x, y, xr, yr, d, Ec=EC, Es=ES) -> float:
    """
    :param yc:
    :param x:
    :param y:
    :param xr:
    :param yr:
    :param d:
    :param Ec:
    :param Es:
    :return: Moment of inertia about the x-axis
    """
    # Get number of vertices
    n = Es / Ec  # Stiffness ratio
    nv = len(x)

    # Create closed polygons by adding the first point to the end of the coordinate lists
    x = x + [x[0]]
    y = y + [y[0]]

    # Convert y-coordinates to specified axis of rotation
    y = [yc - i for i in y]
    yr = [yc - i for i in yr]

    # Compute list of terms for summation

    Icx_list = [1 / 12 * (y[i] ** 2 + y[i] * y[i + 1] + y[i + 1] ** 2) * (x[i] * y[i + 1] - x[i + 1] * y[i]) for i in
                range(nv)]

    # Sum list elements and use absolute value so order can be clockwise or counter-clockwise
    Icx = abs(sum(Icx_list))

    # Create seperate lists for rebas in compression (c) and tension(t)

    _, yr_c, d_c, _, yr_t, d_t = compression_tension_rebars(x, y, xr, yr, d)

    Isx_c = [pi / 64 * d_t[i] ** 4 + n * pi * d_t[i] ** 2 / 4 * yr_t[i] ** 2 for i in range(len(d_c))]

    Isx_t = [pi / 64 * d_t[i] ** 4 + n * pi * d_t[i] ** 2 / 4 * yr[i] ** 2 for i in range(len(d_t))]

    return Icx + sum(Isx_c) + sum(Isx_t)

def InertiaY(xc, x, y, yr, d, Ec=EC, Es=ES):
    """
    :param xc:
    :param x:
    :param y:
    :param yr:
    :param d:
    :param Ec:
    :param Es:
    :return: Return moment of inertia about the y-axis
    """
