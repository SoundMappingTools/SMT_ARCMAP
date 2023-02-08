# -*- coding: utf-8 -*-
'''
Description: This script tests that the present version of ArcGIS is compatible
             with the core sound propagation models
             
Dependencies: ArcGIS >=10.3, Python 2.7, validationhlpr,py, soundprophlpr.py,
              nmsimhlpr.py, nmsimgis_analysis.py, spreadgishlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For B_Validation.py:
Alexander C. Keyel, 3/4/2016 - 4/12/2017
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

# IN ORDER TO RUN THIS SCRIPT FROM PYTHON:
# YOU NEED TO CHANGE tbx_root and testing_dir to match the paths on your computer

# Set toolbox directory
# These will get overwritten if running from the GUI or if the paths are not changed
tbx_root = "C:/Your/Path/Here/" # e.g. C:/smt/toolbox/
testing_dir = "C:/Your/Path/Here/" # e.g., C:/smt/testing/


# Import system module
import sys        

# If sys.argv[1] is undefined, then assume running from Python windowHacky way could be to wrap in try: except:, but that would miss legitimate errors here.
try:
    sys.argv[1]
    debug = 0
except:
    debug = 1

if debug == 0:
    testing_dir = sys.argv[1]    
    tbx_root = sys.argv[2]        
    
    #Correct for ArcGIS's path import practices
    testing_dir = testing_dir + "/" 
    tbx_root = tbx_root + "/"   
    skip_test_run = 0 # Code to skip generating the data. For debugging purposes     
    test_spread = 1 # Code to allow skipping spread-gis test
    test_nmsim = 1  #Code to allow skipping nmsimgis test
    test_iso = 1
    
    # Indictator to get NMSIM to perform properly
    is_GUI = "YES"

if debug == 1:    
    # This step imports a function to get my IDE (Spyder text editor)
    # to talk with ArcGIS's install of Python
    # For troubleshooting purposes, requires the archook function
    # https://github.com/JamesRamm/archook

    try:
        import arcpy
    except:
        import archook
        archook.get_arcpy()    

    # This can be ignored if you set toolbox_dir and testing_dir
    # This is used to automatically assign the paths on the computers used for development
    if tbx_root == "C:/Your/Path/Here/":
        # Assign the path correctly
        import platform
        cpu = platform.processor()
        pathbit = "smt" # Path on my desktop, default path    
        #if cpu != 'Intel64 Family 6 Model 44 Stepping 2, GenuineIntel':
        if cpu == "Intel64 Family 6 Model 61 Stepping 4, GenuineIntel":
            pathbit = "SPreAD-GIS" # Path on my laptop
    
        testing_dir = "C:/%s/testing/" % pathbit
        tbx_root = "C:/%s/toolbox/" % pathbit
        

    skip_test_run = 0
    test_spread = 1
    test_nmsim = 1 
    test_iso = 1
    
    # Indicator to get NMSIMGIS to work properly
    is_GUI = "NO"
                    
# Import & set global settings
import validationhlpr as vh
reload(vh)

# Validate that code runs correctly on this version of ArcGIS
vh.validate_arcgis_version(testing_dir, tbx_root, skip_test_run, test_spread, test_nmsim, test_iso, is_GUI)

    
