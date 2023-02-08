# -*- coding: utf-8 -*-
'''
Description: This script calls a function from nmsimhlpr.py. The code was
placed in a separate file so that it could be called using subprocess.call.
This was necessary to bypass an error where the program ran out of memory as it
progressed through the loop, due to inadequate garbage collection in Python.
The functions it uses are all in nmsimhlpr, this function just calls that file.
             
Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py, nmsimhlpr.py

@author Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For nmsimgis_analysis.py:
Copyright (C) 1/22/2016 - 4/12/2017 A.C. Keyel <skeyel@gmail.com>

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

# Give code access to nmsimgis_analysis function from nmsimhlpr
import sys
import cPickle as pickle
from nmsimhlpr import nmsimgis_analysis as nmsim

pickle_in = sys.argv[1]
args = pickle.load(open(pickle_in, 'rb'))

p = args[0]
n_chunks = args[1]
chunk_paths = args[2]
#print chunk_paths[0]
col_remainder = args[3]
row_remainder = args[4]
col_subset_start = args[5]
row_subset_start = args[6]
freq = args[7]
numpy_dem = args[8]
numpy_impedance = args[9]
source_info = args[10]
temp = args[11]
rh = args[12]
receiver_offset = args[13]
model_dir = args[14]
point_dir = args[15]
dem_name = args[16]
xyzsrc = args[17]
zsrc = args[18]
dem_clip = args[19]
point_id = args[20]
point_fill = args[21]
freq_fill = args[22]
keep_intermediates = args[23]
n_processes = args[24]
nmsim_dll = args[25]
cell_x = args[26]
cell_y = args[27]
xmin = args[28]
xmax = args[29]
ymin = args[30]
ymax = args[31]
no_extent_settings = args[32]
do_timing = args[33]
my_times = args[34]
my_time_labels = args[35]
model = args[36]
out_pickle = args[37]

chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, my_times, my_time_labels = nmsim(p, n_chunks, chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, freq, numpy_dem, numpy_impedance,
                                         source_info, temp, rh, receiver_offset, 
                                         model_dir, point_dir, dem_name,
                                         xyzsrc, zsrc, dem_clip, point_id,
                                         point_fill, freq_fill,
                                         keep_intermediates, n_processes,
                                         nmsim_dll, cell_x, cell_y, xmin,
                                         xmax, ymin, ymax, no_extent_settings, do_timing,
                                         my_times, my_time_labels, model) 


# Pickle results to return them to the main code
# Write results to text file
#print chunk_paths[0]
return_args = [chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, my_times, my_time_labels]
pickle.dump(return_args, open(out_pickle, 'wb'))
#with open(temp_text, 'w') as tt:
#    print >> tt, "%s;%s;%s;%s;%s;%s;%s" % (chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, my_times, my_time_labels)
