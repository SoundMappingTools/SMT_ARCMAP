# -*- coding: utf-8 -*-
'''
Description: Interface for a tool  that assesses focal impacts on a selected
             area or line
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For Focal_Impacts.py:
Copyright (C) 1/12/2017 - 4/11/2017, A.C. Keyel <skeyel@gmail.com>

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

import sys
import soundprophlpr as sp
import arcpy

# Read in script arguments
analysistype = sys.argv[1]
output_file_path = sys.argv[2] + "/"
output_file_name = sys.argv[3]
output_file = output_file_path + output_file_name

# Make sure it has a .shp ending
if output_file[-4:] != ".shp":
    output_file = output_file + ".shp"

# Set focal area/line and threshold
focal_area = sys.argv[4]
threshold = sys.argv[5]


# These are not used/required for the line analysis
point_source_file = sys.argv[6]
source_id_field = sys.argv[7]
data_dir = sys.argv[8] + "/"
data_label = sys.argv[9]

# These are used/required for the line analysis but not the others
dprime_raster = sys.argv[10]
temp_dir = sys.argv[11] + "/"

# Drop excess_adjustment from function call.
n_points = "all"

if analysistype == "sound pressure levels":

    # Run script
    #sp.AssessImpact(point_source_file, source_id_field, n_points, focal_area,
    #                analysis_dir, weighting_id, excess_adjustment = 0,
    #                threshold = 55)
    sp.AssessImpact(output_file, point_source_file, source_id_field, n_points,
                    focal_area, data_dir, data_label, threshold)

if analysistype == "audibility":
    # Check if focal area is a line or not
    info = arcpy.Describe(focal_area)

    if info.shapeType == "Polyline":
        # If line, do the line-based assessment #focal area in this case is the line
        sp.AssessLengthAudible(focal_area, output_file, dprime_raster, temp_dir, threshold = 7)    
    
    elif info.shapeType == "Polygon":
        # If area, do the area-based assessment
        sp.AssessAreaAudible(output_file, point_source_file, source_id_field, n_points, threshold, focal_area, data_dir, data_label)
    else:
        raise ValueError("The code currently only supports Polygons and Polylines. You input %s" % info.shapeType)