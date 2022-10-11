#!/usr/bin/env python3 

'''
To use this script load appropriate python and ncl modules, edit top of MAPL_Tripolar.cdl file,
then run it as:

./mk_MAPL_Tripolar.py MAPL_Tripolar.cdl

'''

import os, sys
import numpy as np
import netCDF4 as nc

mapl_file='MAPL_Tripolar.nc'
cdl_file=sys.argv[1]
os.system(f'ncgen -o {mapl_file} {cdl_file}')

with nc.Dataset(mapl_file) as ff:
    mom_file=ff['mom_file'][:].data.tobytes().decode().replace('\x00','')
    vgrid_file=ff['vgrid_file'][:].data.tobytes().decode().replace('\x00','')
    topo_file=ff['topo_file'][:].data.tobytes().decode().replace('\x00','')
    basin_file=ff['basin_file'][:].data.tobytes().decode().replace('\x00','')

with nc.Dataset(mom_file) as ff:
    # Grid dimension: (jm,im) 
    xx=ff['x'][:]              # (2*jm+1, 2*im+1)
    yy=ff['y'][:]              # (2*jm+1, 2*im+1)
    dx=ff['dx'][:]             # (2*jm+1, 2*im)
    dy=ff['dy'][:]             # (2*jm, 2*im+1)
    angle_dx=ff['angle_dx'][:] # (2*jm+1, 2*im+1)
    area=ff['area'][:]         # (2*jm, 2*im)

with nc.Dataset(vgrid_file) as ff:
    dz=ff['dz'][:]

with nc.Dataset(topo_file) as ff:
    mask=ff['wet'][:]

with nc.Dataset(basin_file) as ff:
    basin=ff['basin'][:]
    atl_mask=np.where(basin==2.0, 1.0, 0.0)

with nc.Dataset(mapl_file, mode='r+') as ff:
    ff['lon_centers'][:]=xx[1::2,1::2]
    ff['lat_centers'][:]=yy[1::2,1::2]
    ff['lon_corners'][:]=xx[0::2,0::2]
    ff['lat_corners'][:]=yy[0::2,0::2]
    ff['lev'][:]=np.cumsum(dz)-0.5*dz
    ff['mask'][:]=mask[:]
    
    ff['htn'][:]=dx[2::2,0::2]+dx[2::2,1::2]
    ff['hte'][:]=dy[0::2,2::2]+dy[1::2,2::2]
    
    ff['huw'][0:-1,:]=dy[1:-1:2,1::2]+dy[2::2,1::2]
    ff['huw'][-1,:]=dy[-1,1::2]+dy[-1,-1:1:-2]
    
    ff['hus'][:,0:-1]=dx[1::2,1:-1:2]+dx[1::2,2::2]
    ff['hus'][:,-1]=dx[1::2,-1]+dx[1::2,0]
    
    ff['areat'][:]=area[0::2,0::2]+area[1::2,0::2]+area[1::2,1::2]+area[0::2,1::2]

    ll=ul=ur=lr=np.zeros(ff['areau'].shape)
    ll=area[1::2,1::2]; # lower left corner
    ul[0:-1,:]=area[2::2,1::2]; ul[-1,:]=area[-1,-1::-2] # upper left corner
    lr[:,0:-1]=area[1::2,2::2]; lr[:,-1]=area[1::2,0] # lower right corner
    ur[0:-1,0:-1]=area[2::2,2::2]
    ur[0:-1,-1]=area[2::2,0]
    ur[-1,0:-1]=area[-1,-2:1:-2]
    ur[-1,-1]=area[-1,0] # upper right corner
    ff['areau'][:]=ll+ul+ur+lr
    
    ff['angleu'][:]=angle_dx[2::2,2::2]
    ff['anglet'][:]=angle_dx[1::2,1::2]

    ff['basin'][:]=basin[:]
    ff['atl_mask'][:]=atl_mask[:]
