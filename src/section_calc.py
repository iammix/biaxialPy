# Built-in libraries
from math import pi, cos, sin, tan, atan, atan2, sqrt, ceil, floor
import logging
import numpy as np
import matplotlib.path as mpltPath
import geometry



# CONSTANTS
EC = 30 * 10 ** 6
ES = 200 * 10 ** 6


def Itx(yc, x, y, xr, yr, d, Ec=EC, Es=ES):
    """
    Return moment of inertia for bending about the x-axis of a transformed reinforced concrete
    cross section. The section can have the shape of any non self-intersecting polygon.
    """

    n = Es / Ec  # Stiffness ratio

    # CONCRETE SECTION
    # Number of vertices
    nv = len(x)

    # Create a closed polygon by adding the first point to the end of the coordinate lists
    x = x + [x[0]]
    y = y + [y[0]]

    # Convert y-coordinates to specified axis of rotation
    y = [yc - i for i in y]
    yr = [yc - i for i in yr]

    # Compute list of terms for summation
    Icx_list = [1 / 12 * (y[i] ** 2 + y[i] * y[i + 1] + y[i + 1] ** 2) *
                (x[i] * y[i + 1] - x[i + 1] * y[i]) for i in range(nv)]

    # Sum list elements and use absolute value so order can be clockwise or counter-clockwise
    Icx = abs(sum(Icx_list))

    # REBARS
    # Create separate lists for rebars in compression (c) and tension (t)
    _, yr_c, d_c, _, yr_t, d_t = compression_tension_rebars(x, y, xr, yr, d)

    # Rebars in compression (areas to be multiplied by 'n-1')
    Isx_c = [pi / 64 * d_c[i] ** 4 + (n - 1) * pi * d_c[i] ** 2 / 4 * yr_c[i] ** 2 for i in range(len(d_c))]

    # Rebars in tension (areas to be multiplied by 'n)
    Isx_t = [pi / 64 * d_t[i] ** 4 + n * pi * d_t[i] ** 2 / 4 * yr_t[i] ** 2 for i in range(len(d_t))]

    return Icx + sum(Isx_c) + sum(Isx_t)


def Ity(xc, x, y, xr, yr, d, Ec=EC, Es=ES):
    """
    Return moment of inertia about the y-axis of a reinforced concrete cross
    section. The section can have the shape of any non self intersecting polygon.
    """

    n = Es / Ec  # Stiffness ratio

    # CONCRETE SECTION
    # Number of vertices
    nv = len(y)

    # Create a closed polygon by adding the first point to the end of the coordinate lists
    x = x + [x[0]]
    y = y + [y[0]]

    # Convert y-coordinates to specified axis of rotation
    x = [xc - i for i in x]
    xr = [xc - i for i in xr]

    # Compute list of terms for summation
    Icy_list = [1 / 12 * (x[i] ** 2 + x[i] * x[i + 1] + x[i + 1] ** 2) *
                (x[i] * y[i + 1] - x[i + 1] * y[i]) for i in range(nv)]

    # Sum list elements and use absolute value so order can be clockwise or counter-clockwise
    Icy = abs(sum(Icy_list))

    # REBARS
    # Create separate lists for rebars in compression (c) and tension (t)
    xr_c, _, d_c, xr_t, _, d_t = compression_tension_rebars(x, y, xr, yr, d)

    # Rebars in compression (areas to be multiplied by 'n-1')
    Isy_c = [pi / 64 * d_c[i] ** 4 + (n - 1) * pi * d_c[i] ** 2 / 4 * xr_c[i] ** 2 for i in range(len(d_c))]

    # Rebars in tension (areas to be multiplied by 'n)
    Isy_t = [pi / 64 * d_t[i] ** 4 + n * pi * d_t[i] ** 2 / 4 * xr_t[i] ** 2 for i in range(len(d_t))]

    return Icy + sum(Isy_c) + sum(Isy_t)


def elastic_centroid(x, y, xr, yr, dia, Ec=EC, Es=ES):
    """
    Return elastic centroid of a transformed reinforced concrete sections.
    Rebars located outside of the concrete defined by x and y is assumed to be
    surrounded by ineffective/crakced concrete.
    """

    # Stiffness ratio
    n = Es / Ec

    # Number of rebars
    nb = len(dia)

    # Evaluate if rebars are inside or outside stress block (returns list with 'True' or 'False')
    rebar_eval = rebars_in_stress_block(x, y, xr, yr)

    # Extract rebars in compression
    dia_comp = [dia[i] for i in range(nb) if rebar_eval[i]]
    xr_comp = [xr[i] for i in range(nb) if rebar_eval[i]]
    yr_comp = [yr[i] for i in range(nb) if rebar_eval[i]]

    # Extract rebars in tension
    dia_tens = [dia[i] for i in range(nb) if not rebar_eval[i]]
    xr_tens = [xr[i] for i in range(nb) if not rebar_eval[i]]
    yr_tens = [yr[i] for i in range(nb) if not rebar_eval[i]]

    # Compute centroid and area of concrete polygon
    xc, yc, Ac = geometry.polygon_centroid(x, y, return_area=True)

    # Compute total transformced area of section
    A_comp = sum([n * pi * d ** 2 / 4 for d in dia_tens])
    A_tens = sum([(n - 1) * pi * d ** 2 / 4 for d in dia_comp])
    A = Ac + A_comp + A_tens

    # Compute total 'moment area', i.e. area times moment arm
    Acx = Ac * xc
    Asx_comp = sum([(n - 1) * pi * dia_comp[i] ** 2 / 4 * xr_comp[i] for i in range(len(dia_comp))])
    Asx_tens = sum([n * pi * dia_tens[i] ** 2 / 4 * xr_tens[i] for i in range(len(dia_tens))])

    Acy = Ac * yc
    Asy_comp = sum([(n - 1) * pi * dia_comp[i] ** 2 / 4 * yr_comp[i] for i in range(len(dia_comp))])
    Asy_tens = sum([n * pi * dia_tens[i] ** 2 / 4 * yr_tens[i] for i in range(len(dia_tens))])

    # Compute x- and y-coordinate of elastic centroid for transformed section
    xel = (Acx + Asx_comp + Asx_tens) / A
    yel = (Acy + Asy_comp + Asy_tens) / A

    return xel, yel


def transformed_axial_stiffness(x, y, xr, yr, dia, P, Ec=EC, Es=ES):
    """ Return axial stiffness EA of transformed concrete section. """

    # Stiffness ratio
    n = Es / Ec

    # Area of rebars
    As = sum([pi * d ** 2 / 4 for d in dia])

    if P <= 0:
        # Axial force is compressive

        # Compute area of section
        A = geometry.polygon_area(x, y)

        # Area of concrete
        Ac = A - As

        # Transformed stiffness weighed by actual area (not transformed area in this case)
        Et = (Ec * Ac + (n - 1) * Es * As) / A

        # Transformed arae
        At = Ac + (n - 1) * As

        return Et * At

    else:
        # Axial force is tensile

        # Stiffness and area contribution comes from rebars only
        E = Es

        return E * As


def strain_field_eval(x, y, P, Mx, My, E, EA, Itx, Ity):
    """
    Return the evaluation of the strain field equation given for external loads
    P, Mx and My in point (x, y).
    """

    # Axial strain
    eps_P = P / EA

    # Curvature from bending about x- and y-axis
    kappa_x = Mx / (E * Itx)
    kappa_y = My / (E * Ity)

    return eps_P + y * kappa_x + x * kappa_y


def compute_plastic_centroid(x, y, xr, yr, As, fck, fyk):
    """ Return plastic centroid of a reinforced concrete section. """
    Ac = geometry.polygon_area(x, y)
    eta = 0.85
    F = sum([As[i] * fyk for i in As]) + eta * (Ac - sum(As)) * fck

    # TODO: [biaxialPy] Find correct and general arm for concrete force
    # assignees: iammix
    # labels: todo

    F_dx = sum([As[i] * fyk * xr[i] for i in range(len(xr))]) + eta * (Ac - sum(As)) * fck * 500 / 2
    F_dy = sum([As[i] * fyk * yr[i] for i in range(len(yr))]) + eta * (Ac - sum(As)) * fck * 375 / 2

    xpl = F_dx / F
    ypl = F_dy / F

    return xpl, ypl


def compression_tension_rebars(x, y, xr, yr, dia):
    """
    Return lists of rebar coordinates and diameters for rebars in compression
    and tension, respectively.
    """
    # Evaluate if rebars are inside or outside stress block (returns list with 'True' or 'False')
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


def compute_dist_to_na(x, y, xr, yr, alpha_deg, na_y):
    """ Return distances from neutral axis to all concrete section vertices
        and rebars. """
    # Convert input angle from [deg] to [rad]
    alpha = alpha_deg * pi / 180

    # Define two known points on line representing neutral axis
    na_x0 = 0
    na_y0 = tan(alpha) * na_x0 + na_y
    na_x1 = 1
    na_y1 = tan(alpha) * na_x1 + na_y

    # Compute signed distances from neutral axis to each vertex (neg. value => vertex in compr. / pos. value => vertex in tension)
    dv = [geometry.point_to_line_dist(x[i], y[i], na_x0, na_y0, na_x1, na_y1) for i in range(len(x))]

    # Compute signed distances from neutral axis to each rebar
    dr = [geometry.point_to_line_dist(xr[i], yr[i], na_x0, na_y0, na_x1, na_y1) for i in range(len(xr))]

    # Reverse sign of the signed distances if slope of neutral axis becomes negative
    if alpha_deg > 90 and alpha_deg <= 270:
        dv = list(np.negative(dv))
        dr = list(np.negative(dr))

    # Change potential distances of '-0.0' to '0.0' to avoid getting the wrong cross section state later
    dv = [0.0 if x == -0.0 else x for x in dv]

    return dv, dr


def stress_block_geometry(x, y, dv, dr, alpha_deg, na_y, lambda_=0.8):
    """
    Returns stress block geometry.
    """


    # PURE TENSION CASE
    # NOTE Test if this is true! Does not account for gap btw. sb and tension zone
    if all(d >= 0 for d in dv):
        cross_section_state = 'PURE TENSION'

        # Distance from neutral axis to extreme tension bar (all distances will be positve)
        c = max([d for d in dr if d > 0])

        # Set vertices of stress block
        x_sb = None
        y_sb = None

        # Set stress block area
        Asb = 0

        sb_cog = None

    # PURE COMPRESSION CASE
    elif all(d <= 0 for d in dv):  # NOTE Test if this is true!
        cross_section_state = 'PURE COMPRESSION'

        # Distance from neutral axis to extreme compression fiber (all distances will be negative)
        c = min(dv)

        # Set vertices of stress block (entire section)
        x_sb = x
        y_sb = y

        Asb = geometry.polygon_area(x, y)
        sb_cog = geometry.polygon_centroid(x, y)

    else:
        cross_section_state = 'MIXED TENSION/COMPRESSION'

        # Distance from neutral axis to extreme compression fiber (pos. in tension / negative in compression)
        c = min(dv)
        a = lambda_ * c

        # Signed perpendicular distance between neutral axis and stress block
        delta_p = c - a

        # Vert. dist. in y-coordinate from neutral axis to inner edge of stress block
        delta_v = delta_p / cos(alpha_deg * pi / 180)

        # Intersection between stress block inner edge and y-axis (parallel with neutral axis)

        sb_y_intersect = na_y - delta_v

        # Intersections between inner edge of stress block (parrallel with neutral axis) and section
        sb_xint, sb_yint = geometry.line_polygon_collisions(alpha_deg, sb_y_intersect, x, y)

        # Find concrete section vertices that are in compression
        x_compr_vertices, y_compr_vertices = geometry.get_section_compression_vertices(x, y, na_y, alpha_deg, delta_v)

        # Collect all stress block vertices
        x_sb = sb_xint + x_compr_vertices
        y_sb = sb_yint + y_compr_vertices

        # Order stress block vertices with respect to centroid for the entire section
        x_sb, y_sb = geometry.order_polygon_vertices(x_sb, y_sb, x, y, counterclockwise=True)

        # Compute area of the stress block by shoelace algorithm
        Asb = geometry.polygon_area(x_sb, y_sb)

        # Compute location of centroid for stress block polygon
        sb_cog = geometry.polygon_centroid(x_sb, y_sb)

    return x_sb, y_sb, Asb, sb_cog, c


def compute_rebar_strain(dist_to_na, c, eps_cu):
    """ Return strain in each rebar as a list """
    return [ri / abs(c) * eps_cu for ri in dist_to_na]


def compute_rebar_stress(eps_r, Es, fyd):
    """ Return stress in each rebar as a list """

    sigma_r = []
    for i in range(len(eps_r)):
        # Linear elastic stress in i'th bar
        si = eps_r[i] * Es

        # Check if rebar yields
        if abs(si) <= fyd:
            # Computed stress does not exceed yield stress
            sigma_r.append(si)
        else:
            # Computed stress exceeds yield, use yield stress instead
            sigma_r.append(np.sign(si) * fyd)

    return sigma_r


def rebars_in_stress_block(x_sb, y_sb, xr, yr):
    ''' Returns a list with entry 'True' for rebars located inside the stress block, 'False' otherwise '''

    if xr and yr:
        # Arrange rebar coordinates
        rebar_coords = [[xr[i], yr[i]] for i in range(len(xr))]
    else:
        raise ValueError('No rebars in section.')

    # Compute area of stress block
    Asb = geometry.polygon_area(x_sb, y_sb)

    if Asb != 0:
        # Arrange stress block coordinates
        sb_poly = [[x_sb[i], y_sb[i]] for i in range(len(x_sb))]

        # Check if rebars are inside the stress block
        path = mpltPath.Path(sb_poly)
        # Returns 'True' if rebar is inside stress block
        rebars_inside = path.contains_points(rebar_coords)
    else:
        # All rebars are in tension (all entries are 'False')
        rebars_inside = [False] * len(xr)

    return rebars_inside


def compute_rebar_forces(xr, yr, As, sigma_r, rebars_inside, fcd, lambda_=0.80):
    """ Return rebar forces as list."""
    Fr = []  # Forces in each rebar

    for i in range(len(xr)):
        if rebars_inside[i]:
            # Rebar is inside stress block, correct for disp. of concrete
            Fi = (sigma_r[i] + lambda_ * fcd) * As
        else:
            # Rebar is outside stress block
            Fi = sigma_r[i] * As
        Fr.append(Fi)

    return Fr


def compute_concrete_force(fcd, Asb, lambda_=0.80):
    """ Return compression force in the concrete. """
    Fc = -lambda_ * fcd * Asb
    return Fc


def compute_moment_vector_angle(Mx, My):
    """    Return the angle (in degrees) of the moment vector with respect to the x-axis    """
    if Mx == 0:
        if My == 0:
            phi = None
        else:
            phi = 90
    else:
        phi = atan(My / Mx) * 180 / pi

    return phi


def compute_C_T_forces(Fc, Fr):
    """    Return Compression (C) and Tension (T) forces of the section    """
    Fr_compr = [p for p in Fr if p <= 0]
    Fr_tension = [p for p in Fr if p > 0]
    C = sum(Fr_compr) + Fc
    T = sum(Fr_tension)

    return C, T


def compute_moment_contributions(xr, yr, Asb, sb_cog, Fc, Fr):
    """ Return the moment contributions from concrete and rebars in the cross section. """

    if Asb == 0:
        Mcx = 0
        Mcy = 0
    else:
        # Moment contribution from concrete about x-axis
        Mcx = -Fc * sb_cog[1]
        # Moment contribution from concrete about y-axis
        Mcy = -Fc * sb_cog[0]

    # Moment contribution from rebars about x- and y-axis (according to moment sign convention)
    Mrx = [-Fr[i] * yr[i] for i in range(len(yr))]
    Mry = [-Fr[i] * xr[i] for i in range(len(xr))]

    return Mcx, Mcy, Mrx, Mry


def compute_C_T_moments(C, T, Mcx, Mcy, Mry, Mrx, Fr, alpha_deg):
    """
    Return total moments generated in the section by Compression (C) and Tension (T) resisting forces.
    """
    My_compr = []
    Mx_compr = []
    My_tension = []
    Mx_tension = []
    for i in range(len(Fr)):
        if Fr[i] <= 0:
            My_compr.append(Mry[i])
            Mx_compr.append(Mrx[i])
        if Fr[i] > 0:
            My_tension.append(Mry[i])
            Mx_tension.append(Mrx[i])

    # Total moment for compression resisting forces (adapted for LH sign convention)
    if 90 <= alpha_deg <= 270:
        My_C = sum(My_compr) + Mcy
        Mx_C = sum(Mx_compr) + Mcx
    else:
        My_C = -(sum(My_compr) + Mcy)
        Mx_C = -(sum(Mx_compr) + Mcx)

    # Total moment for tension resisting forces (adapted for LH sign convention)
    if 90 <= alpha_deg <= 270:
        My_T = sum(My_tension)
        Mx_T = sum(Mx_tension)
    else:
        My_T = -sum(My_tension)
        Mx_T = -sum(Mx_tension)

    return Mx_C, My_C, Mx_T, My_T


def compute_C_T_forces_eccentricity(C, T, My_C, Mx_C, Mx_T, My_T):
    """    Return eccentricity of Compression (C) and Tension (T) forces.    """
    # Eccentricities of tension and compression forces
    if C == 0:
        ex_C = np.nan
        ey_C = np.nan
    else:
        ex_C = My_C / C
        ey_C = Mx_C / C

    if T == 0:
        ex_T = np.nan
        ey_T = np.nan
    else:
        ex_T = My_T / T
        ey_T = Mx_T / T

    return ex_C, ey_C, ex_T, ey_T


def perform_section_analysis(x, y, xr, yr, fcd, fyd, Es, eps_cu, As, alpha_deg, na_y, lambda_=0.80):
    """ Perform cross section analysis """


    dv, dr = compute_dist_to_na(x, y, xr, yr, alpha_deg, na_y)
    x_sb, y_sb, Asb, sb_cog, c = stress_block_geometry(x, y, dv, dr, alpha_deg, na_y, lambda_=lambda_)
    eps_r = compute_rebar_strain(dr, c, eps_cu)
    sigma_r = compute_rebar_stress(eps_r, Es, fyd)
    rebars_inside = rebars_in_stress_block(x_sb, y_sb, xr, yr)
    Fr = compute_rebar_forces(xr, yr, As, sigma_r, rebars_inside, fcd, lambda_=lambda_)
    Fc = compute_concrete_force(fcd, Asb)


    return Fc, Fr, Asb, sb_cog, x_sb, y_sb
