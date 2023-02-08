# -*- coding: utf-8 -*-
'''
Description: Interface script to convert from Sound Pressure Levels to an
             Leq measure
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For ConvertToLeq.py:
Copyright (C) 9/23/2016 A.C. Keyel <skeyel@gmail.com>

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

# Import standard library modules
import sys
# Import custom soundprophlpr module
import soundprophlpr as sp
import arcpy

arcpy.AddWarning("This tool is still being developed. Use at your own risk")

# Get parameters from ArcGIS toolbox
list_of_rasters = sys.argv[1]
sound_source_duration = sys.argv[2]        
Leq_duration = sys.argv[3]
results_path = sys.argv[4]

# Convert (weighted) SPL data to Leq
sp.convert_to_Leq_batch(list_of_rasters, sound_source_duration, Leq_duration,
                  results_path) 