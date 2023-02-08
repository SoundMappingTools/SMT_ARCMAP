'''
Description: Interface for a function call that creates a raster dataset of
              ambient sound conditions by land cover type
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For A_AmbientSoundConditions.py:
Alexander C. Keyel, 2016
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

if __name__ == "__main__":
    import sys

    # Import the file containing the Create_Ambient function
    import soundprophlpr as sp

    try:
        sys.argv[1]
        debug = 0
    except:
        debug = 1
    
    if debug == 0:
        # Input parameters
        single_freq = sys.argv[1]
        freq_table = sys.argv[2]
        veg = sys.argv[3]
        in_table = sys.argv[4]
        wind_speed = sys.argv[5]
        snow = sys.argv[6]
        ambient_path = sys.argv[7]
        
        ambient_path = ambient_path + "/" #Add forward slash    

    if debug == 1:
        cpu = "desktop"
        cpu = "laptop"
        pathbit = "smt"
        if cpu == "laptop":
            pathbit = "SPreAD-GIS"
        
        spath = "C:/%s/source_data/" % pathbit
        single_freq = "400"
        freq_table = "#"
        veg = spath + "nlcd"
        in_table = spath + "Test_Background.csv"
        wind_speed = "#"
        snow = "#"
        ambient_path = spath

    
    if single_freq != "#" and freq_table != "#":
        raise ValueError("Either single frequency or frequency table should be"
                         "specified, not both. You entered %s" % single_freq +
                         "for single frequency and %s for fruency table" % freq_table)
    
    if single_freq == "#" and freq_table == "#":
        raise ValueError("Either single frequency or frequency table MUST be specified")
    
    if single_freq != "#":
        in_freq = single_freq
    if freq_table != "#":
        in_freq = freq_table
    
    sp.CreateAmbient(in_freq, veg, in_table, wind_speed, snow, ambient_path)

