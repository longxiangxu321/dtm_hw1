# -- my_code_hw01.py
# -- hw01 GEO1015.2022
# -- [Longxiang Xu]
# -- [5722918]

import random
import math

import numpy as np
import scipy


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


def get_area(args):
    t_area = 0
    for i in range(len(args) - 1):
        t_area = t_area + (args[i][0] * args[i + 1][1]) - (args[i + 1][0] * args[i][1])

    t_area = t_area + args[-1][0] * args[0][1] - args[0][0] * args[-1][1]
    # t_area = abs(t_area) / 2
    t_area = t_area / 2
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


def get_vor_cell_size(dt, center):
    adj_triangles = dt.incident_triangles_to_vertex(center)
    v1_cell = []
    for triangle in adj_triangles:
        # breakpoint()
        x1 = dt.points[triangle][0][0:2]
        y1 = dt.points[triangle][1][0:2]
        z1 = dt.points[triangle][2][0:2]
        vd_vertex = get_circumcenter(x1, y1, z1)
        v1_cell.append(vd_vertex)
    vor_cell_size = get_area(v1_cell)
    return vor_cell_size


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

    # 先找到插入点的所有的相接点，再分别计算相接点和所有相接点的三角形的外接圆圆心，构成voronoi cell的到Ai
    # 插入点，找到所有相接点，计算所有外接圆圆心，计算面积A0
    # 找到插入点的所有相接点，再分别计算相接点的所有相接点的外接圆圆心，构成voronoi cell的到Bi
    # (Ai - Bi) / A = weight
    # 实际顺序要反过来一下,因为得先insert再remove

    # Insert point, find the adjacent vertices and calculate their voronoi cell size, Ai
    # Calculate the input point's voronoi cell size A
    # remove input point, then calculate the adjacent voronoi cell size Bi
    # weight = (Bi - Ai) / A
    in_pt = dt.insert_one_pt(x, y, z)
    adj_ver = dt.adjacent_vertices_to_vertex(in_pt)
    vo_cell_size = get_vor_cell_size(dt, in_pt)  # A

    # obtain the cell_size of the adjacent cells of the input point cell
    area_after_insertion = {}
    for pt in adj_ver:
        vor_cell_size = get_vor_cell_size(dt, pt)  # Bi
        area_after_insertion[pt] = vor_cell_size

    dt.remove(in_pt)

    # obtain the cell_size of the adjacent cells of the input point cell (before insertion)
    area_before_insertion = {}
    for pt in adj_ver:
        vor_size = get_vor_cell_size(dt, pt)  # Ai
        area_before_insertion[pt] = vor_size

    weights = {}
    for pt in adj_ver:
        weights[pt] = (area_before_insertion[pt] - area_after_insertion[pt]) / vo_cell_size

    wei = list(weights.values())
    z_val = []
    for pt in adj_ver:
        z_val.append(dt.points[pt][2])
    z = np.average(z_val, weights=wei)
    return z
