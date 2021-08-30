# Built-in packages
from math import pi
import logging

# Third party packages
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Project specific packages
import section_calc as sc
import section_plots
from geometry import point_to_point_dist_3d
from geometry import line_hull_intersection
from scipy.optimize import linprog


# TODO: [biaxialPy] Neutral axis rotation should about plastic centroid, see 'Structural Analysis of Cross Sections', p. 190
# assignees: iammix
# labels: todo



np.set_printoptions(precision=2)


def compute_capacities(Fc, Fr, Mcx, Mcy, Mrx, Mry):
    """    Returns capacities P, Mx and My    """
    # Total capacities
    P = sum(Fr) + Fc
    Mx = sum(Mrx) + Mcx
    My = sum(Mry) + Mcy

    return P, Mx, My


def compute_capacity_surface(x, y, xr, yr, fcd, fyd, Es, eps_cu, As, lambda_=0.80,  rotation_step=5, vertical_step=10):
    """ Returns coordinates for capacity surface of cross section"""
    # TODO: [biaxialPy] Find a good way to define steps and loop over entire function
    # assignees: iammix
    # labels: todo
    vs = vertical_step
    h = max(x) - min(x)
    na_y_list = list(np.linspace((min(x)-h/3), 0, vs)) + list(np.linspace(0, (max(x)+h/3), vs))
    alpha_list = [alpha for alpha in range(0, 360, rotation_step)]

    P_list = []
    Mx_list = []
    My_list = []
    na_y_computed = []
    alpha_computed = []
    for na_y in na_y_list:
        for alpha_deg in alpha_list:

            # Perform cross section ULS analysis
            Fc, Fr, Asb, sb_cog, _, _ = sc.perform_section_analysis(x, y, xr, yr, fcd, fyd, Es, eps_cu, As, alpha_deg, na_y, lambda_=lambda_)

            # Compute individual moments generated in the section
            Mcx, Mcy, Mrx, Mry = sc.compute_moment_contributions(xr, yr, Asb, sb_cog, Fc, Fr)

            # Compute capacities
            P, Mx, My = compute_capacities(Fc, Fr, Mcx, Mcy, Mrx, Mry)

            # Update lists of calculated pairs of vertical local and anlge for neutral axis
            na_y_computed.append(na_y)
            alpha_computed.append(alpha_deg)

            # Store iteration results
            P_list.append(P)
            Mx_list.append(Mx)
            My_list.append(My)

    return P_list, Mx_list, My_list, na_y_computed, alpha_computed





def in_point_cloud(points, x):
    n_points = len(points)
    n_dim = len(x)
    c = np.zeros(n_points)
    A = np.r_[points.T, np.ones((1, n_points))]
    b = np.r_[x, np.ones(1)]
    lp = linprog(c, A_eq=A, b_eq=b)
    return lp.success


def utilization_ratio(Ped, Mxed, Myed, P_capsurf, Mx_capsurf, My_capsurf):
    """
    Return the utilization ratio as the ratio between the distance from the load
    combination point to Origo and the distance from Origo to capacity surface in
    dirrection given by the point.
    """

    # Compute convex hull of the capacity surface point cloud
    cap_surf = np.transpose(np.array([P_capsurf, Mx_capsurf, My_capsurf]))
    convex_hull = ConvexHull(cap_surf)

    # Compute distance from Origo to all load combination points
    comb_array = np.transpose(np.array([Ped, Mxed, Myed]))
    lc_dist = [point_to_point_dist_3d([0, 0, 0], lc) for lc in comb_array]

    # Compute distance from Origo to all capacity surface points in line with load combinations
    cap_intersections = [line_hull_intersection(lc, convex_hull) for lc in comb_array]
    cap_dist = [point_to_point_dist_3d([0, 0, 0], cap) for cap in cap_intersections]

    # A for loop is made here in favor of a list comprehension in order to catch
    # the case where all loads in combination is 0. In that case, the Vector
    # has no direction.
    # Compute utilization ratios for all load combinations
    ur = []
    for i in range(len(comb_array)):
        # if not np.all(np.nonzero(comb_array[i])):
        if np.all(comb_array[i] == 0):
            # All external loads are 0, i.e. P=0, Mx=0 and My=0
            ur.append(0.00)
        else:
            # Some loads are nonzero
            ur.append(lc_dist[i] / cap_dist[i])

    return ur

