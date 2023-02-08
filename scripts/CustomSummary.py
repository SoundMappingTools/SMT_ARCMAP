# -*- coding: utf-8 -*-
'''
Description: Interface for a tool that takes the output from NMSIMGIS or
SPreAD-GIS for a list of model run results and creates a weighted overall
summary based on the traffic volume for each model result. For example, if the 
noise impacts of a single car traveling down a road were assessed, this tool
could be used to scale up the result to 500 cars. As a second example, if the
model was run for a single car and a single heavy truck, this tool could assess
the impact (SPL) of 200 cars and 10 heavy trucks.
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For CustomSummary.py:
Copyright (C) 9/22/2016 - 4/12/2017 A.C. Keyel <skeyel@gmail.com>

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

# Get parameters from ArcGIS toolbox
list_of_model_results = sys.argv[1]        
list_of_weights = sys.argv[2]
results_path = sys.argv[3] + "/" # Convert to path
results_name = sys.argv[4] 

list_of_model_results = list_of_model_results.split(";") # Split into raster paths
list_of_weights = map(int, list_of_weights.split(";")) # Split & convert to float

#import arcpy
#arcpy.AddMessage(list_of_model_results)
#arcpy.AddMessage(list_of_weights)

# Run custom summary function
sp.custom_summary(list_of_model_results, list_of_weights, results_path,
                  results_name)