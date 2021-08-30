from math import sqrt, pi, tan, atan2

import numpy as np


def line_polygon_collisions(angle_deg, y_intersect, x_vertex, y_vertex) -> tuple:
    """
    :param angle_deg: Angle of line with x-axis (in radians)
    :param y_intersect: Intersection between line and y-axis
    :param x_vertex: x-coordinates of polygon vertices
    :param y_vertex: y-coordinates of polygon vertices
    :return: Return intersection points between a line and a polygon. If no intersections are present, return original polygon vertices.
    """
    angle = angle_deg * pi / 180
    vertex_eval = [tan(angle) * x_vertex[i] - y_vertex[i] + y_intersect for i in range(len(x_vertex))]

    count_pos = 0
    count_neg = 0
    for vertex in vertex_eval:
        if vertex > 0:
            count_pos += 1
        else:
            count_neg += 1

    if count_neg == len(vertex_eval):
        return x_vertex, y_vertex
    elif count_pos == len(vertex_eval):
        return x_vertex, y_vertex

    else:
        vertex_eval.append(vertex_eval[0])
        xv = x_vertex + [x_vertex[0]]
        yv = y_vertex + [y_vertex[0]]

        xint = []
        yint = []
        for i in range(len(vertex_eval) - 1):
            if np.sign(vertex_eval[i]) != np.sign(vertex_eval[i + 1]):

                x1 = xv[i]
                y1 = yv[i]
                x2 = xv[i + 1]
                y2 = xv[i + 1]

                if x1 == x2:
                    x = x1
                    y = tan(angle) * x + y_intersect
                    xint.append(x)
                    yint.append(y)
                else:
                    m = (y2 - y1) / (x2 - x1)
                    k = y2 - m * x2

                    x = (y_intersect - k) / (m - tan(angle))
                    y = m * x + k
                    xint.append(x)
                    yint.append(y)

    return xint, yint


def polygon_area(x, y, signed=False) -> float:
    """
    :return: Return the area of a non-intersecting polygon given the coordinates of its vertices
    """
    if x is not None and y is not None:
        x = x + [x[0]]
        y = y + [y[0]]

        a1 = [x[i] * y[i + 1] for i in range(len(x) - 1)]
        a2 = [y[i] * x[i + 1] for i in range(len(y) - 1)]
        if signed:
            A = 1 / 2 * (sum(a1) - sum(a2))
        else:
            A = 1 / 2 * abs(sum(a1) - sum(a2))
    else:
        A = 0

    return A


def point_to_line_dist(x, y, x0, y0, x1, y1, signed=True) -> float:
    """
    :return: The distance from a point (x,y) to a line passing through points (x0,y0) and (x1, y1)
    """

    if signed:
        return -((y0 - y1) * x + (x1 - x0) * y + (x0 * y1 - x1 * y0)) / sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    else:
        return abs((y0 - y1) * x + (x1 - x0) * y + (x0 * y1 - x1 * y0)) / sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)


def polygon_centroid(x, y, return_area=False) -> tuple:
    """
    :return: Compute the centroid of a non-self-intersecting polygon given the coordinates of its vertices
    """
    A = polygon_area(x, y, signed=True)
    if A == 0:
        return np.nan

    else:
        x = x + [x[0]]
        y = y + [y[0]]
        cx = []
        cy = []

        for i in range(len(x) - 1):
            cx.append((x[i] + x[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i]))
            cy.append((y[i] + y[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i]))

        Cx = sum(cx) / (6 * A)
        Cy = sum(cy) / (6 * A)
        if return_area:
            return Cx, Cy, A
        else:
            return Cx, Cy


def order_polygon_vertices(x_vertices, y_vertices, x_section_vertices, y_section_vertices,
                           counterclockwise=True) -> tuple:
    """
    :return: Sort polygon vertices in consecutive circular order (clockwise or counterclockwise measured form positive x-axis)
    """

    x_t = x_vertices
    y_t = y_vertices

    Cx, Cy = polygon_centroid(x_section_vertices, y_section_vertices)

    a0 = []
    for i in range(len(x_t)):
        a0.append(atan2((y_t[i] - Cy), (x_t[i] - Cx)) * 180 / pi)
    if counterclockwise:
        pos = [p for p in a0 if p >= 0]
        neg = [n for n in a0 if n < 0]
        a = pos + [360 + angle for angle in neg]

        idx = sorted(range(len(a0)), key=lambda j: a0[j], reverse=False)
        a = sorted(a)
    else:
        # TODO: [biaxialPy] Sort list of angles in clockwise order
        # assignees: iammix
        # labels: todo
        pass

    x_t = [x_t[i] for i in idx]
    y_t = [y_t[i] for i in idx]

    return x_t, y_t


def sb_eq_eval(angle, na_y, y_shift, x, y) -> float:
    """
    :return: evaluation of equation for inner stress block at point (x, y)
    """
    return tan(angle) * x - y + na_y - y_shift


def get_section_compression_vertices(x, y, na_y, alpha_deg, delta_v) -> tuple:
    """
    :return: Returns a list of the concrete section vertices that are in compression
    """
    x_compr_vertices = []
    y_compr_vertices = []
    alpha = alpha_deg * pi / 180

    for i in range(len(x)):
        sb_eq_eval_at_each_vertex = sb_eq_eval(alpha, na_y, delta_v, x[i], y[i])

        if alpha_deg <= 90 or alpha_deg > 270:
            if sb_eq_eval_at_each_vertex < 0:
                x_compr_vertices.append(x[i])
                y_compr_vertices.append(y[i])

        if alpha_deg > 90 and alpha_deg <= 270:
            if sb_eq_eval_at_each_vertex > 0:
                x_compr_vertices.append(x[i])
                y_compr_vertices.append(y[i])
    return x_compr_vertices, y_compr_vertices


def point_to_point_dist_3d(P1, P2) -> float:
    return sqrt((P2[0] - P1[0]) ** 2 + (P2[1] - P2[1]) ** 2 + (P2[2] - P1[2]) ** 2)


def line_hull_intersection(U, c_hull) -> float:
    eq = c_hull.equations.T
    V, b = eq[:, -1], eq[-1]
    V = np.transpose(V)
    alpha = -b / np.dot(V, U)
    return np.min(alpha[alpha > 0]) * U
