# -*- coding: utf-8 -*-
'''
Description: Interface for the tool designed to take source information and
convert it to a NMSim .src file. Source information can be in the Frequency
Table format used by SPreAD-GIS, or can be raw field data from a Larson Davis
recorder (Larson Davis input is still in alpha testing)
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For ConvertToSRC.py:
Copyright (C) 9/25/2016 A.C. Keyel <skeyel@gmail.com>

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
csv_lst = sys.argv[1].split(';')
file_type = sys.argv[2]
analysis_type = sys.argv[3]
outpath = sys.argv[4]
source_name = sys.argv[5]
source_description = sys.argv[6]
avg_name = sys.argv[7]
source_type = sys.argv[8]
AGL_m = float(sys.argv[9])
ref_dist_m = float(sys.argv[10])

outpath = outpath + "/"

#**# SPEED VEC & THINK ABOUT OTHER ASPECTS TO DIRECTIONALITY

# Call soundprophlpr function to actually make the .src file
sp.make_srcfile(csv_lst, file_type, analysis_type, outpath, source_name,
                source_description, avg_name, source_type, AGL_m, ref_dist_m)
