#-- geo1015_hw01.py
#-- hw01 GEO1015.2022
#-- Hugo Ledoux <h.ledoux@tudelft.nl>
#-- 2022-09-22

#------------------------------------------------------------------------------
# DO NOT MODIFY THIS FILE!!!
#------------------------------------------------------------------------------

import sys
import math
import csv
import json 
import time

import numpy as np
import scipy.spatial  #-- for kd-tree & nearest neighbour(s) queries
import startinpy      #-- for a Delaunay triangulation

import raster
import my_code_hw01

NO_DATA_VALUE = -99999

def main():
    #-- read the needed parameters from the file 'params.json' (must be in same folder)
    try:
        jparams = json.load(open('params.json'))
    except:
        print("ERROR: something is wrong with the params.json file.")
        sys.exit()
    #-- store the input 3D points in a DT
    dt = startinpy.DT()
    #-- cleaning of duplicates done in the process with tolerance of 1cm
    dt.snap_tolerance = 0.01
    with open(jparams['input-file']) as csvfile:
        r = csv.reader(csvfile, delimiter=' ')
        header = next(r)
        totall = 0
        for line in r:
            p = list(map(float, line)) #-- convert each str to a float
            assert(len(p) == 3)
            dt.insert_one_pt(p[0], p[1], p[2])
            totall += 1
        if totall > dt.number_of_vertices():
            print("INFO: {} duplicate points were removed".format(totall - dt.number_of_vertices()))
    #-- fetch all the (clean) points (see https://startinpy.readthedocs.io/en/latest/api.html#startinpy.DT.points)
    pts = dt.points[1:] 
    #-- construct a KD-tree also, for fast nearest neighbours queries
    kd = scipy.spatial.KDTree(pts[:,:2]) 
    all_z = pts[:,-1]
    #-- find bbox, we get bbox[minx,miny,maxx,maxy]   
    bbox = dt.get_bbox()
    # breakpoint()
    for key in jparams:
        if key != 'input-file':
            print("=== {} interpolation ===".format(key))
            start_time = time.time()
            myraster = raster.Raster(jparams[key]['cellsize'], bbox)
            i = 0
            for row in range((myraster.height - 1), -1, -1):
                j = 0
                y = myraster.bbox[1] + (row * myraster.cellsize) + (myraster.cellsize / 2)
                for col in range(myraster.width):
                    x = myraster.bbox[0] + (col * myraster.cellsize) + (myraster.cellsize / 2)
                    if   key == 'nn':
                        try:
                            myraster[i][j] = my_code_hw01.nn_xy(dt, kd, all_z, x, y)
                        except Exception as e:
                            myraster[i][j] = NO_DATA_VALUE
                    elif key == 'tin':
                        try:
                            myraster[i][j] = my_code_hw01.tin_xy(dt, kd, all_z, x, y)
                        except Exception as e:
                            myraster[i][j] = NO_DATA_VALUE
                    elif key == 'nni':
                        try:
                            myraster[i][j] = my_code_hw01.nni_xy(dt, kd, all_z, x, y)
                        except Exception as e:
                            myraster[i][j] = NO_DATA_VALUE
                    elif key == 'idw':
                        try:
                            myraster[i][j] = my_code_hw01.idw_xy(dt, kd, all_z, x, y, 
                                                            jparams[key]['power'], 
                                                            jparams[key]['radius'])
                        except Exception as e:
                            myraster[i][j] = NO_DATA_VALUE
                    j += 1
                i += 1 
            myraster.save(jparams[key]['output-file'])
            print("-->%ss" % round(time.time() - start_time, 2))        


if __name__ == '__main__':
    main()





