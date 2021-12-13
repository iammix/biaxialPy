import pytest
import src.geometry as geometry


def test_line_polygon_intersections():
    angle = 30
    y_intersect = -2

    x_vertex = [-8, 8, 8, -8]
    y_vertex = [8, 8, -8, -8]

    x_int, y_int = geometry.line_polygon_collisions(angle, y_intersect, x_vertex, y_vertex)

    x_int = [round(i, 2) for i in x_int]
    y_int = [round(i, 2) for i in y_int]

    assert x_int == [8.0, -8.0] and y_int == [2.62, -6.62]
