# -*- coding: utf-8 -*-
'''
Description: Interface that can be used to run the NMSIMGIS sound propagation
             model
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py, nmsimhlpr.py
              nmsimgis_analysis.py, NumPy

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For D_nmsimgis.py:
Copyright (C) 2016, A.C. Keyel <skeyel@gmail.com>

Alexander "Sasha" Keyel
Postdoctoral Researcher
1474 Campus Delivery
Colorado State University
Fort Collins, CO 80523-1474
  

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License 2.0 as published by the Free
Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.

'''

#import multiprocessing
#http://stackoverflow.com/questions/18204782/runtimeerror-on-windows-trying-
#python-multiprocessing
# If re-enabling multiprocessing you need if __name__== '__main__':

# Import system module
import sys, ntpath        


# If sys.argv[1] is undefined, then assume running from Python window
try:
    sys.argv[1]
    debug = 0
except:
    debug = 1

if debug == 0:
    import soundprophlpr as sp
    import arcpy
 
   # Input parameters (optional parameters will come in as '#' if left blank)
    model = "nmsimgis"
    tbx_root = sys.argv[1]        
    #tbx_root = "C:/smt/"
    model_dir = sys.argv[2]
    results_dir = model_dir + "/nmsimgis/"
    point_source_file = sys.argv[3] 
    # Results file will have the same name as point_source_file
    results_label = ntpath.basename(point_source_file)[:-4] 
    source_id_field = sys.argv[4]
    Model_Extent = sys.argv[5]
    elevation = sys.argv[6]
    landcover = sys.argv[7]
    srcfile = sys.argv[8]  # Full path to the source file
    #head = float(sys.argv[9]) # Vehicle heading (deg. True North)
    # Vehicle roll angle (deg., right wing down is positive roll)
    #roll = float(sys.argv[10])
    #pitch = float(sys.argv[11])  # Vehicle pitch angle (deg., nose up is positive
    #vel = float(sys.argv[12])  # Velocity (m/s)
    #engpow = float(sys.argv[13])  # Engine power (%)
    temp = sys.argv[9]
    rh = sys.argv[10]
    #source_offset = float(sys.argv[16])
    receiver_offset = float(sys.argv[11])
    allow_defaults = sys.argv[12]    
    summarize_string = sys.argv[13]
    summarize = sp.unpack_summary(summarize_string)    
    ambient = sys.argv[14]
    keep_intermediates_string = sys.argv[15]
    keep_intermediates = sp.unpack_ki(keep_intermediates_string)
    drop_down_weight = sys.argv[16]      

    #Correct for ArcGIS's path import practices
    tbx_root = tbx_root + "/"        
    model_dir = model_dir + "/"
    if ambient == "#" or ambient == "":
        ambient = "NA"
    # No longer needed, as ambient is now always a file name (in fact this breaks the downstream code!)
    #if ambient != "NA":
    #    ambient = ambient + "/"
    
    # Add extra variable that does not apply
    timelog = "none"
    
    # Add variable that controls how GUI works
    is_GUI = "YES"
   
# Run from Python        
if debug == 1:    
    # This step imports a function to get my IDE (Spyder text editor)
    # to talk with ArcGIS's install of Python
    # For troubleshooting purposes, requires the archook function
    # https://github.com/JamesRamm/archook
    import archook
    archook.get_arcpy()    
    import soundprophlpr as sp
    import arcpy
    import nmsimhlpr
    reload(sp)
    reload(nmsimhlpr)

    import platform
    cpu = platform.processor()
    pathbit = "smt" # Path on my desktop, default path    
    if cpu != 'Intel64 Family 6 Model 44 Stepping 2, GenuineIntel':
        pathbit = "SPreAD-GIS" # Path on my laptop
    
    # For trouble-shooting/testing. Note: inputs must be in text format
    spath = "C:/%s/source_data/" % pathbit #Source path
    #spath = "C:/%s/source_data_old/" % pathbit #Source path
    
    # Specify a base directory
    model = "nmsimgis" #"all"        
    tbx_root = "C:/%s/toolbox/" % pathbit
    model_dir = "C:/data/"
    model_dir = "C:/data/"
    #model_dir = "C:/data/nmsimgis_sensitivity/"
    results_dir = model_dir + "nmsimgis/"
    results_label = "OHV"
    #results_label = "C182"
    #point_source_file = spath + "OHV_trail_point.shp"
    point_source_file = spath + "ATV2.shp"
    #point_source_file = spath + "C182.shp"
    source_id_field = "POINT"
    Model_Extent = spath + "model_extent.shp"
    #Model_Extent = spath + "model_extent_really_small.shp"
    Model_Extent = sp.get_extent(Model_Extent)
    elevation = spath + "dem"    
    #elevation = spath + "dem_m"
    landcover = spath + "nlcd"
    ambient = "NA" # spath + "ambient/"
    #ambient = spath + "ambient/ambient00400"
    # Full path to the source file
    srcfile = tbx_root + "Sources/HarrisonSources/ATV2.src"  
    #srcfile = tbx_root + "Sources/AirTourFixedWingSources/C182.src"
    temp = "28"    
    #temp = "15" 
    #temp = "10;20;5"
    rh = "26"    
    #rh = "70"
    #rh = "50;70;10"
    receiver_offset = 1.22
    allow_defaults = "YES" # Allow ArcGIS to fill in any missing fields in the shapefile
    
    #timelog = "none"
    timelog = model_dir + "nmsimgis_timing.csv"
    summarize = 1
    keep_intermediates = 1
    drop_down_weight = "A"
    is_GUI = "NO"

# Set a point to truncate the results
truncate = 0 # default truncate at 0 dB

# Convert weighting to appropriate code for SoundPropagation
weighting = sp.convert_to_weighting(drop_down_weight)         

# Add unused spreadgis parameters to make function-call happy
wind_dir = wind_sp = seas_cond = "NA"
 
# Now this is read in from the shapefile  
#source_info = [head, roll, pitch, vel, engpow, srcfile]
n_points = "all" # Option to run for fewer points. Only an option in the script

# Frequency values must be taken from .src file
#multiprocessing.freeze_support()

# Check whether a sensitivity analysis is desired
test = temp.split(";")

# If a single value is given for temperature, just run the main analysis
if len(test) == 1:
    arcpy.AddMessage("Running normal nmsimgis analysis")
    sp.SoundPropagation(model, point_source_file, source_id_field,
                        Model_Extent, ambient, elevation, landcover, temp, rh,
                        wind_dir, wind_sp, seas_cond, model_dir, results_dir,
                        results_label, srcfile, 
                        receiver_offset, n_points, tbx_root, timelog,
                        summarize = summarize,
                        keep_intermediates = keep_intermediates,
                        weighting = weighting, truncate = truncate,
                        landcover_type = "nlcd", is_GUI = is_GUI,
                        allow_defaults = allow_defaults)


# Otherwise, run the sensitivity analysis #**# STILL UNDER DEVELOPMENT                         
'''
if len(test) == 3:
    min_temp = float(test[0])
    max_temp = float(test[1])
    temp_increment = float(test[2])
    temp_extremes = [min_temp, max_temp]
    rh_values = rh.split(';')
    min_rh = float(rh_values[0])
    max_rh = float(rh_values[1])
    rh_increment = float(rh_values[2])
    rh_extremes = [min_rh, max_rh]        
    seas_cond_lst = [seas_cond]
    base_output_dir = model_dir
    timelog = "none"
    delete_existing_timelog = 1
    results_dir = model_dir # Results will go in the main model directory
    
    #**# FIX THIS
    d_max = 4000
    d_increment = 1000    
    
    arcpy.AddMessage("Running nmsimgis sensitivity analysis")
    sp.Sensitivity_Analysis(model, point_source_file, source_id_field,
                            Model_Extent, ambient, elevation, landcover,
                            d_max, d_increment,
                            temp_extremes, rh_extremes, temp_increment,
                            rh_increment, wind_sp, seas_cond_lst,
                            base_output_dir, results_dir, srcfile,
                            receiver_offset,
                            n_points, tbx_root, timelog,
                            delete_existing_timelog, keep_intermediates,
                            weighting, truncate = truncate)   
'''