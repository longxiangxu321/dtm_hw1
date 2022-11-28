# -- my_code_hw01.py
# -- hw01 GEO1015.2022
# -- [Longxiang Xu]
# -- [5722918]

import random
import math

import numpy as np


def nn_xy(dt, kd, all_z, x, y):
    """
    !!! TO BE COMPLETED !!!
     
    Function that interpolates with nearest neighbour method.
     
    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points 
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
    Output:
        z: the estimation of the height value, 
           (raise Exception if outside convex hull)
    """
    # -- kd-tree docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    # -- you are *not* allowed to use the function for the nn interpolation that I wrote for startinpy
    # -- you need to write your own code for this step
    z = random.uniform(0, 100)
    try:
        index_nn = dt.closest_point(x, y)
        z = dt.points[index_nn][2]
        return z
    except Exception:
        raise Exception("Outside convex hull")


def idw_xy(dt, kd, all_z, x, y, power, radius):
    """
    !!! TO BE COMPLETED !!!
     
    Function that interpolates with IDW
     
    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points 
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
        power:  power to use for IDW
        radius: search radius
¨    Output:
        z: the estimation of the height value, 
           (raise Exception if (1) outside convex hull or (2) no point in search radius
    """
    # -- kd-tree docs: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.html#scipy.spatial.KDTree
    z = random.uniform(0, 100)
    try:
        index_nn = dt.closest_point(x, y)  # determine whether input point is inside the convex hull
    except Exception:
        raise Exception("Outside convex hull")
    nn_s = kd.query_ball_point([x, y], radius)
    if len(nn_s) == 0:
        raise Exception("No point in search radius")
    weights = []
    nns_z = []
    for nn in nn_s:
        weight = pow(math.dist([x, y], kd.data[nn]), -power)
        weights.append(weight)
        nns_z.append(all_z[nn])

    nns_z = np.array(nns_z)
    weights = np.array(weights)
    z = np.average(nns_z, weights=weights)
    return z


def tin_xy(dt, kd, all_z, x, y):
    """
    !!! TO BE COMPLETED !!!
     
    Function that interpolates linearly in a TIN.
     
    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points 
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
    Output:
        z: the estimation of the height value, 
           (raise Exception if outside convex hull)
    """
    # -- startinpy docs: https://startinpy.rtfd.io/
    # -- you are *not* allowed to use the function for the TIN interpolation that I wrote for startinpy
    # -- you need to write your own code for this step
    z = random.uniform(0, 100)
    try:
        index_nn = dt.closest_point(x, y)  # determine whether input point is inside the convex hull
    except Exception:
        raise Exception("Outside convex hull")

    ind = dt.locate(x, y)  # indexes of vertices of the triangle containing input point
    x0, y0, z0 = dt.points[ind[0]][0], dt.points[ind[0]][1], dt.points[ind[0]][2]
    x1, y1, z1 = dt.points[ind[1]][0], dt.points[ind[1]][1], dt.points[ind[1]][2]
    x2, y2, z2 = dt.points[ind[2]][0], dt.points[ind[2]][1], dt.points[ind[2]][2]

    Ae = [[x0, x1, x2], [y0, y1, y2], [1, 1, 1]]
    be = [[x], [y], [1]]
    weights = scipy.linalg.solve(Ae, be)

    z = z0 * weights[0][0] + z1 * weights[1][0] + z2 * weights[2][0]
    return z


def get_area(*args):
    t_area = 0
    for i in range(len(args) - 1):
        t_area = t_area + (args[i][0] * args[i + 1][1]) - (args[i + 1][0] * args[i][1])

    t_area = t_area + args[-1][0] * args[0][1] - args[0][0] * args[-1][1]
    t_area = abs(t_area) / 2
    return t_area


def get_circumcenter(x, y, z):
    ax = x[0]
    ay = x[1]
    bx = y[0]
    by = y[1]
    cx = z[0]
    cy = z[1]
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d
    return [ux, uy]


def nni_xy(dt, kd, all_z, x, y):
    """
    !!! TO BE COMPLETED !!!
     
    Function that interpolates with natural neighbour interpolation method (nni).
     
    Input:
        dt:     the DT of the input points (a startinpy object)
        kd:     the kd-tree of the input points 
        all_z:  an array with all the z values, same order as kd.data
        x:      x-coordinate of the interpolation location
        y:      y-coordinate of the interpolation location
    Output:
        z: the estimation of the height value, 
           (raise Exception if outside convex hull)
    """
    # -- startinpy docs: https://startinpy.rtfd.io/
    # -- you are *not* allowed to use the function for the nni interpolation that I wrote for startinpy
    # -- you need to write your own code for this step
    z = random.uniform(0, 100)
    try:
        index_nn = dt.closest_point(x, y)  # determine whether input point is inside the convex hull
    except Exception:
        raise Exception("Outside convex hull")
    tri = dt.locate(x, y)  # return index of the three vertices (the three vertices forms the triangle contains the
    # input)

    p0 = dt.points[tri[0]]
    p1 = dt.points[tri[1]]
    p2 = dt.points[tri[2]]

    p0_n = dt.adjacent_vertices_to_vertex(tri[0])
    p1_n = dt.adjacent_vertices_to_vertex(tri[1])
    p2_n = dt.adjacent_vertices_to_vertex(tri[2])

    p01 = dt.points[np.setdiff1d(np.intersect1d(p0_n, p1_n), tri[2])[0]]
    p02 = dt.points[np.setdiff1d(np.intersect1d(p0_n, p2_n), tri[1])[0]]
    p12 = dt.points[np.setdiff1d(np.intersect1d(p1_n, p2_n), tri[0])[0]]

    v_01 = get_circumcenter(p0, p1, p01)
    v_02 = get_circumcenter(p0, p2, p02)
    v_12 = get_circumcenter(p1, p2, p12)
    v_012 = get_circumcenter(p0, p1, p2)

    new_pt = dt.insert_one_pt(x, y, z)

    v_0_01 = get_circumcenter(p0, p01, [x, y])
    v_1_01 = get_circumcenter(p1, p01, [x, y])

    v_0_02 = get_circumcenter(p0, p02, [x, y])
    v_2_02 = get_circumcenter(p2, p02, [x, y])

    v_1_12 = get_circumcenter(p1, p12, [x, y])
    v_2_12 = get_circumcenter(p2, p12, [x, y])

    a_01 = get_area(v_0_01, v_1_01, v_01)
    a_02 = get_area(v_0_02, v_2_02, v_02)
    a_12 = get_area(v_1_12, v_2_12, v_12)

    a0 = get_area(v_0_01, v_01, v_012, v_02, v_0_02)
    a1 = get_area(v_01, v_1_01, v_1_12, v_12, v_012)
    a2 = get_area(v_012, v_12, v_2_12, v_2_02, v_02)

    area = get_area(v_0_01, v_1_01, v_1_12, v_2_12, v_2_02, v_0_02)
    area_1 = a_01 + a_02 + a_12 + a0 + a1 + a2
    print(area - area_1)
    # breakpoint()

    weights = [a0 / area, a1 / area, a2 / area, a_01 / area, a_02 / area, a_12 / area]
    z_val = [p0[2], p1[2], p2[2], p01[2], p02[2], p12[2]]
    z = np.average(z_val, weights=weights)
    # 首先，找到输入点包含的三角形的三个顶点locate(), 找到与该三角形直接相接的三个三角形顶点adjacent_vertices_to_vertex()
    # 输出当中重复的点即为顶点三角形的每条边及其相接的三角形的顶点，求外接圆的圆心，即为voronoi vertex
    # 对新点，以及先前的三个点连线分别计算外接圆的圆心，即为新的voronoi vertex
    # 用方程求解新旧voronoi vertex连成的图形
    # 重复下一个


    # 先找到插入点的所有的相接点，再分别计算相接点的所有相接点的外接圆圆心，构成voronoi cell的到Ai
    # 插入点，找到所有相接点，计算所有外接圆圆心，计算面积A0
    # 找到插入点的所有相接点，再分别计算相接点的所有相接点的外接圆圆心，构成voronoi cell的到Bi
    # Ai - Bi / A = weight
    dt.remove(new_pt)
    # std = dt.interpolate_nni(x, y)
    # error = z - std
    # print(error)
    return z
