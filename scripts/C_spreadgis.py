# -*- coding: utf-8 -*-
'''
Description: Interface for a function call that runs the SPreAD-GIS sound 
             propagation model
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py, spreadgishlpr.py
              NumPy

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For C_spreadgis.py:
Copyright (C) 2010 - 2016 Sarah E. Reed and A.C. Keyel <skeyel@gmail.com>

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

# Import system module
import sys, ntpath     

# If sys.argv[1] is undefined, then assume running from Python window
try:
    sys.argv[1]
    debug = 0
except:
    debug = 1
    
    try:
        import arcpy
    except:
        # This step imports a function to get my IDE (Spyder text editor) to
        # talk with ArcGIS's install of Python
        # archook is from https://github.com/JamesRamm/archook    
        # If arcpy is on the PYTHONPATH, do not need archook.
        # Otherwise, use archook to add arcpy to PYTHONPATH (temporary)
        import archook
        archook.get_arcpy()    

import soundprophlpr as sp

if debug == 0:
    # Input parameters (optional parameters will come in as '#' if left blank)
    model = "spreadgis"
    tbx_root = sys.argv[1]        
    model_dir = sys.argv[2]
    results_dir = model_dir + "/spreadgis/"
    point_source_file = sys.argv[3] 
    # ntpath.basename just gets the file name, [:-4] scrubs the .shp extension
    # Results will have the same name as input points file
    results_label = ntpath.basename(point_source_file)[:-4] 
    source_id_field = sys.argv[4]
    Model_Extent = sys.argv[5]
    elevation = sys.argv[6]
    landcover = sys.argv[7]
    in_freq = sys.argv[8]    
    #in_Sound_Level = sys.argv[16] # Often included in the frequency table
    #in_Measurement_Distance = sys.argv[17] #Often in the frequency table
    temp = sys.argv[9]
    rh = sys.argv[10]
    wind_dir = sys.argv[11]
    wind_sp = sys.argv[12]
    seas_cond = sys.argv[13]
    summarize_string = sys.argv[14]
    summarize = sp.unpack_summary(summarize_string)
    ambient = sys.argv[15]
    keep_intermediates_string = sys.argv[16]
    keep_intermediates = sp.unpack_ki(keep_intermediates_string)
    drop_down_weight = sys.argv[17]
    
    #Correct for ArcGIS's path import practices
    tbx_root = tbx_root + "/"        
    model_dir = model_dir + "/"

    # Correct for missing optional inputs
    if ambient == "#" or ambient == "NA":
        ambient = "NA"
    # else is no longer needed, as ambient is now ALWAYS a file name.
    #else:
    #    ambient = ambient + "/"
    

# Run in Python        
if debug == 1:    
    reload(sp)
    

    # For trouble-shooting/testing. Note: inputs must be in text format

    #**# To change the absolute paths between my laptop and desktop
    cpu = "desktop"
    #cpu = "laptop"
    pathbit = "smt"
    if cpu == "laptop":
        pathbit = "SPreAD-GIS"

    spath = "C:/%s/source_data/" % pathbit #Source path
    tbx_root = "C:/%s/toolbox/" % pathbit
    
    # Specify a base directory
    model_dir = "C:/data/"
    #model_dir = "C:/data/spreadgis_sensitivity/"
    results_dir = model_dir + "spreadgis/"
    results_label = "default"
    
    model = "spreadgis" #"all"        
    point_source_file = spath + "OHV_trail_point.shp"
    source_id_field = "POINT"
    Model_Extent = spath + "model_extent.shp"
    Model_Extent = sp.get_extent(Model_Extent)
    in_freq = spath + "Frequency_Table.csv"
    elevation = spath + "dem"
    landcover = spath + "nlcd"
    temp = "28" #.3333333333"
    rh = "30"
    #temp = "20;30;5"
    #rh = "20;40;10"
    wind_dir = "105"  
    wind_sp = "12" #.8748" 
    seas_cond = "clear, windy summer day"
    # A path-based approach to ambient conditions makes automating over
    # multiple frequencies easier.
    ambient = spath + "ambient/ambient00400" #"NA"
    summarize = 1
    keep_intermediates = 1
    drop_down_weight = "A"
        
#Inputs needed for nmsimgis but not SPreAD-GIS
receiver_offset = 1
   
# Add option to truncate at a specific dB level
truncate = 0 # Set to truncate final output at 0 dB   
   
# Convert weighting to the format required by sp.SoundPropagation            
weighting = sp.convert_to_weighting(drop_down_weight)         

# Option to run for fewer points. Only an option in the script            
n_points = "all" 

# Check whether a sensitivity analysis is desired
test = temp.split(";")

# If a single value is given for temperature, just run the main analysis
if len(test) == 1:
    arcpy.AddMessage("Running normal SPreAD-GIS analysis")
    sp.SoundPropagation(model, point_source_file, source_id_field,
                        Model_Extent, ambient, elevation, landcover, temp, rh,
                        wind_dir, wind_sp, seas_cond, model_dir, results_dir,
                        results_label, in_freq,
                        receiver_offset, n_points, tbx_root,
                        summarize = summarize,
                        keep_intermediates = keep_intermediates,
                        weighting = weighting, truncate = truncate,
                        landcover_type = "nlcd", use_old_barrier = 0)
         
# Otherwise, run a sensitivity analysis (#**# still under development)
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
    results_dir = model_dir # Keep the default simple

    #**# NEED TO THINK ABOUT HOW BEST TO DERIVE THESE - 1/2 extent for distance, and what about for increment???
    #**# Need to implement that derivation
    d_max = 4000
    d_increment = 1000

    arcpy.AddMessage("Running SPreAD-GIS sensitivity analysis")
    sp.Sensitivity_Analysis(model, point_source_file, source_id_field,
                            Model_Extent, ambient, elevation, landcover,
                            d_max, d_increment,
                            temp_extremes, rh_extremes, temp_increment,
                            rh_increment, wind_sp, seas_cond_lst, 
                            base_output_dir, results_dir, in_freq,
                            receiver_offset,
                            n_points, tbx_root, timelog,
                            delete_existing_timelog, keep_intermediates,
                            weighting, truncate = truncate)
'''