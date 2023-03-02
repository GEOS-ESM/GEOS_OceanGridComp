#!/usr/bin/env python3 

'''
To use this script load appropriate python and ncl modules, edit top of MAPL_Tripolar.cdl file,
then run it as:

./mkMAPLTripolar.py MAPL_Tripolar.cdl

'''

import os, sys
import numpy as np
import netCDF4 as nc

mapl_file='MAPL_Tripolar.nc'
cdl_file=sys.argv[1]
os.system(f'ncgen -o {mapl_file} {cdl_file}')

with nc.Dataset(mapl_file) as ff:
    grid_spec_file=ff['grid_spec_file'][:].data.tobytes().decode().replace('\x00','')
    basin_file=ff['basin_file'][:].data.tobytes().decode().replace('\x00','')

with nc.Dataset(grid_spec_file) as mom:
    lon_centers=mom['x_T'][:]
    lat_centers=mom['y_T'][:]
    lon_corners=mom['x_C'][:]
    lat_corners=mom['y_C'][:]
    lev=mom['zt'][:]
    mask=mom['wet'][:]
    htn=mom['ds_02_22_T'][:]*100
    hte=mom['ds_20_22_T'][:]*100
    hus=mom['ds_00_20_C'][:]*100
    huw=mom['ds_00_02_C'][:]*100
    areat=mom['area_T'][:]
    areau=mom['area_U'][:]

with nc.Dataset(basin_file) as ff:
    basin=ff['AMOC_MASK'][:]

with nc.Dataset(mapl_file, mode='r+') as mapl:
    mapl['lon_centers'][:] = lon_centers[:]
    mapl['lat_centers'][:] = lat_centers[:]
    mapl['lon_corners'][:] = lon_corners[:]
    mapl['lat_corners'][:] = lat_corners[:]
    mapl['lev'][:] = lev[:]
    mapl['mask'][:] = mask[:]
    mapl['htn'][:] = htn[:]
    mapl['hte'][:] = hte[:]
    mapl['hus'][:] = hus[:]
    mapl['huw'][:] = huw[:]
    mapl['areat'][:]=areat[:]
    mapl['areau'][:]=areau[:]
    mapl['atl_mask'][:]=atl_mask[:]

