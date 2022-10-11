`mk_MAPL_Tripolar.py` utily creates `MAPL_Tripolar.nc` file from MOM6 mosaic grid. It takes as an argument `MAPL_Tripolar.cdl` file. `MAPL_Tripolar.cdl` file should edited before running `mk_MAPL_Triporal.py`. Metadata section in `MAPL_Tripolar.cdl` defines grid dimentions, variable names, precision and attributes. Example of metadata section:

```
dimensions:
	string = 255 ;
        n_center_x = 1440 ;
        n_center_y = 1080 ;
        n_corner_x = 1441 ;
        n_corner_y = 1081 ;
	n_levels = 75 ;

variables:
        float lon_centers(n_center_y,n_center_x) ;
                lon_centers:long_name = "longitude of cell center" ;
                lon_centers:units = "degrees_E" ;
        float lat_centers(n_center_y,n_center_x) ;
                lat_centers:long_name = "latitude cell center" ;
                lat_centers:units = "degrees_N" ;
.................................................................
	double areat(n_center_y,n_center_x) ;
	        areat:long_name = "area of T-Cell" ;
		areat:units = "m2" ;
	double areau(n_center_y,n_center_x) ;
	        areau:long_name = "area of U-Cell" ;
		areau:units = "m2" ;
```

Data section contains path to input files. Example of data section

```
data:
	mom_file = "/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_newtopo/INPUT/ocean_hgrid.nc" ;
	vgrid_file = "/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_newtopo/INPUT/vgrid.nc" ;
	topo_file = "/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_newtopo/INPUT/ocean_topog.nc" ;
	basin_file = "/discover/nobackup/projects/gmao/ssd/aogcm/ocean_bcs/MOM6/1440x1080_newtopo/INPUT/basin_codes.nc" ;
```

