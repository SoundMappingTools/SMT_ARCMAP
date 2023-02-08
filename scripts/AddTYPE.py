'''
Description: This tool adds the SPREADTYPE and/or NMSIMTYPE field for an input
             landcover
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For AddTYPE.py:
Alexander C. Keyel, 10/4/2016 - 4/12/2017
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
landcover = sys.argv[1]        
landcover_type = sys.argv[2]
model = sys.argv[3]

if str.lower(model) == "spread-gis" or str.lower(model) == "both":
    # Run the Add_SPREADTYPE tool
    sp.Add_SPREADTYPE(landcover, landcover_type)
    
if str.lower(model) == "nmsimgis" or str.lower(model) == "both":
    sp.Add_NMSIMTYPE(landcover, landcover_type)