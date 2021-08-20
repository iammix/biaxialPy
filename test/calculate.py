# Built-in packages
from math import pi, cos, sin, tan, atan, atan2, sqrt, ceil, floor
import logging

# Third party packages
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Wedge, Polygon
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Project specific packages
import section_calc as sc
import section_plot_uls
from geometry import point_to_point_dist_3d
from geometry import line_hull_intersection