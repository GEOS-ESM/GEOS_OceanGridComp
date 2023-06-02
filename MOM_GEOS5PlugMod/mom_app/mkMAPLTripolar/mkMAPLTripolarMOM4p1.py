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
    lon_centers=mom['geolon_t'][:]
    lat_centers=mom['geolat_t'][:]
    lon_corners=mom['geo_lonc'][:]
    lat_corners=mom['geo_latc'][:]
    lev=mom['zt'][:]
    mask=mom['wet'][:]
    htn=mom['dxtn'][:]*100
    hte=mom['dyte'][:]*100
    hus=mom['dxte'][:]*100
    huw=mom['dytn'][:]*100

with nc.Dataset(basin_file) as ff:
    basin=ff['tmask'][:]
    atl_mask=np.where(basin==2.0, 1.0, 0.0)

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
    mapl['basin'][:]=basin[:]
    mapl['atl_mask'][:]=atl_mask[:]

