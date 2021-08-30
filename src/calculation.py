from math import pi
import numpy as np
from scipy.spatial import ConvexHull
import geometry
import section_calc as sc


def find_na(x, y, xr, yr, dia, P, Mx, My, fyd, Ec=30*10**6, Es=200*10**6):
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
        EtAt = sc.transformed_axial_stiffness(x, y, xr, yr, dia, P, Ec, Es)

        # Compute distances from concrete vertices and rebars to neutral axis
        dv, dr = sc.compute_dist_to_na(x, y, xr, yr, 0, yn)

        # Get coordinates of most compressed concrete fiber
        # If section is only in tension this will be the most tensioned fibre

        xmax = x[dv.index(min(dv))]
        ymax = y[dv.index(min(dv))]

        # Compute elastic centroid
        xel, yel = sc.elastic_centroid(x, y, xr, yr, dia, Ec=Ec, Es=Es)

        # Compute transformed moment of inertia
        Itx = sc.Itx(yel, x, y, xr, yr, dia, Ec=Ec, Es=Es)
        Ity = sc.Ity(xel, x, y, xr, yr, dia, Ec=Ec, Es=Es)

        # Compute value for max compression strain within the section strain field
        # In section is only in tension this will be the min tension strain

        eps_max = sc.strain_field_eval(xmax, ymax, P, Mx, My, Ec, EtAt, Itx, Ity)

        # NOTE: test
        eps_max = 0.0035
        # NOTE: Test

        # Compute geometry of the concrete stress block

        # Get list of indices for stress block and extract corresponding x and y coordinates
        idx_sb = np.where(np.array(dv) <= 0)[0].tolist()
        dv_sb = [dv[i] for i in idx_sb]
        x_sb = [x[i] for i in idx_sb]
        y_sb = [y[i] for i in idx_sb]

        dv_sb += [0, 0]
        intersections = geometry.line_polygon_collisions(0, yn, x, y)
        x_sb += [intersections[0][0], intersections[0][1]]
        y_sb += [intersections[1][0], intersections[1][1]]

        # Compute force resultant of stress block (volume of stress block)

        # Stresses and strains at vertices

        eps_sb = [di / min(dv_sb) * eps_max for di in dv_sb]
        sigma_sb = [eps_i * Ec for eps_i in eps_sb]

        # Construct points defining volume of stress block
        x_vol = x_sb * 2
        y_vol = y_sb * 2
        z_vol = [0] * len(x_sb) + sigma_sb

        # Compute force resultant (volume)
        points = np.transpose(np.array([x_vol, y_vol, z_vol]))
        Fc = ConvexHull(points).volume

        # Compute rebar strain
        eps_r = sc.compute_rebar_strain(dr, min(dv), eps_max)

        # Compute rebar stress
        sigma_r = sc.compute_rebar_stress(eps_r, Es, fyd)

        # Compute rebar forces
        Fr_each = [sigma_r[i] * pi * dia[i] ** 2 / 4 for i in range(len(sigma_r))]

        Fr = sum(Fr_each)
        print("eps_sb= ", eps_sb)
        print("Fc= ", Fc)
        print("eps_r= ", eps_r)
        print("sigma_r= ", sigma_r)
        print("Fr_each= ", Fr_each)
        print("Fr =", Fr)

        # Compute total force P

        P_guess = Fc + Fr

        # Error between computed axial load P and extrernally applied P
        yn_error = P - P_guess

        if yn_error >= 0:
            yn += tol
        else:
            yn = yn - 0.01

        yn_previous = yn

        print('')
        print("max strain = ", eps_max)
        print("Number of iterations: ", itr_yn)
        print("P = ", P)
        print("P_guess = ", P_guess)
        print("yn_error = ", yn_error)
        print("yn = ", yn)
        print("TRUTH CHECK: ", abs(yn_error) > tol and itr_yn < max_itr)

    return yn


if __name__ == '__main__':
    b = 0.250
    h = 0.500
    dia = [pi * 0.020 ** 2 / 4] * 3  # Rebar diameters [in] (No. 7 bars)
    c = 0.040
    Ec = 33
    Es = 200
    fyd = 500 / 1.15

    x = [0, b, b, 0]
    y = [0, 0, h, h]
    xr = [b / 3, b / 2, 2 / 3 * b]
    yr = [c, c, c]

    P = -80
    Mx = 91
    My = 0

    print('yn =', find_na(x, y, xr, yr, dia, P, Mx, My, fyd, Ec=Ec, Es=Es))
