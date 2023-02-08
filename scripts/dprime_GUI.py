# -*- coding: utf-8 -*-
"""
Description: This script enables a GUI interface for dprime_raster.py, which is
a script used for calculating audibility using 1/3 octave spectra

Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py, dprime_raster.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For dprime_GUI.py:
Copyright (C) 2/17/2017 - 4/12/2017 A.C. Keyel <skeyel@gmail.com>

Alexander "Sasha" Keyel
Postdoctoral Researcher
1474 Campus Delivery
Colorado State University
Fort Collins, CO 80523-1474

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License 2.0 as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""

import sys

# Determine whether input is coming from GUI or not. If name == main would probably work for this too!
try:
    sys.argv[1]
    debug = 0
except:
    debug = 1

# Read in values from GUI
if debug == 0:
    point_lst = sys.argv[1]
    point_id = sys.argv[2]
    band_input = sys.argv[3]

    # Convert band input to actual bands
    if band_input == "NMSIMGIS":
        bands1 = [50, 63, 80, 100, 125, 160, 200, 250]
        bands2 = [315, 400, 500, 630, 800, 1000, 1250, 1600]
        bands3 = [2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000]
        bands = bands1 + bands2 + bands3
    
    elif band_input == "SPREADGIS":
        bands = [125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000]
    elif band_input == "HARRISON":
        bands = [400, 500, 630, 800, 1000, 1250, 1600, 2000]
    else:
        bands = [int(band_input)]
    
    onethird_dir = sys.argv[4]    + "/" 
    band_fill = int(sys.argv[5])
    point_fill = int(sys.argv[6])
    ambient_dir = sys.argv[7]
    dprime_path = sys.argv[8]
    ambient_dir = ambient_dir + "/"
    dprime_path = dprime_path + "/"
    n_points = "all"


# Read in values indicated below        
if debug == 1:
    import archook
    archook.get_arcpy()    
        
    point_lst = "path/to/shapefile.shp"
    point_id = "FID"
    #bands = [125]
    bands1 = [50, 63, 80, 100, 125, 160, 200, 250]
    bands2 = [315, 400, 500, 630, 800, 1000, 1250, 1600]
    bands3 = [2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000]
    bands = bands1 + bands2 + bands3
    
    
    band_fill = 5
    point_fill = 3
    n_points = "all"
    onethird_dir = "C:/data/compressor_daily/nmsimgis/frequency_propagation/"
    ambient_dir = "C:/data/dprime_test/ambient/"
    dprime_path = "C:/data/dprime_test/"
    

# Run analysis
# Normally import statements are at the top, but if debug == 1, I need to import archook first.
# if debug == 0, archook is not necessary, and this will cause problems for people missing this module
import dprime_raster as dpr 
# from dprime_raster import *
dpr.calculate_dprime(point_lst, point_id, n_points, point_fill, bands, band_fill, onethird_dir, ambient_dir, dprime_path)
#point_source_file, source_id_field, n_points, point_fill, bands, band_fill, onethird_dir, ambient_dir, dprime_path


