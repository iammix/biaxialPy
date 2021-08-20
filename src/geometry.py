from math import sqrt, pi, cos, sin, tan, atan, atan2
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
