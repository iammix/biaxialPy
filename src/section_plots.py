from math import pi, atan
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import geometry
import section_calc as sc


def plot_capcity_surface(X, Y, Z, plot_type="scatter", labels=['Mx', 'My', 'P']):
    """ Plots capacity surface as 3D graph"""

    fig_surface = plt.figure()
    ax = Axes3D(fig_surface)

    X = [i / 10 ** 6 for i in X]
    Y = [i / 10 ** 6 for i in Y]
    Z = [i / 10 ** 6 for i in Z]

    if plot_type == 'trisurf':
        surf = ax.plot_trisurf(X, Y, Z, linewidth=0.2, antialiased=True)

    if plot_type == 'wireframe':
        wire = ax.plot_wireframe(X, Y, Z, linewidth=0.2, antialiased=True)
    else:
        scat = ax.scatter(X, Y, Z, linewidth=0.2, antialiased=True)

    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])

    plt.title('Capacity surface')
    plt.show()


def plt_ULS_section(x, y, xr, yr, x_sb, y_sb, Asb, sb_cog, Fc, Fr, Mcx, Mcy, Mrx, Mry, Mx, My, alpha_deg, na_y):
    """    Returns a plot of ULS section state for given neutral axis location    """
    phi = sc.compute_moment_vector_angle(Mx, My)
    C, T = sc.compute_C_T_forces(Fc, Fr)
    Mx_C, My_C, Mx_T, My_T = sc.compute_C_T_moments(C, T, Mcx, Mry, Mrx, Fr, alpha_deg)
    ex_C, ey_C, ex_T, ey_T = sc.compute_C_T_forces_eccentricity(C, T, My_C, Mx_C, Mx_T, My_T)

    # Find collision points between neutral axis and concrete sections
    na_xint, na_yint = geometry.line_polygon_collisions(alpha_deg, na_y, x, y)

    fig, ax = plt.subplots()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.style.use('seaborn-white')

    # Full concrete section
    x_plot = x
    x_plot.append(x_plot[0])
    y_plot = y
    y_plot.append(y_plot[0])
    plt.plot(x_plot, y_plot, '-', color='k', linewidth=1)

    plt.plot(na_xint, na_yint, color='k', linewidth=1)
    plt.annotate('$C$', (ex_C, ey_C), color='b', horizontalalignment='left', verticalalignment='center', fontsize=12)
    plt.annotate('$T$', (ex_T, ey_T), color='b', horizontalalignment='left', verticalalignment='center', fontsize=12)

    # Plot stress blocks
    if Asb != 0:
        sb_coords = []
        for i in range(len(x_sb)):
            sb_coords.append([x_sb[i], y_sb[i]])
        ax.add_patch(patches.Polygon((sb_coords), facecolor='silver', edgecolor='k', linewidth=1))

    # Plot centroid of stress block
    if Asb != 0:
        plt.plot(sb_cog[0], sb_cog[1], 'x', color='grey', markersize='4')
        # TODO: [
