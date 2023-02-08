# -*- coding: utf-8 -*-
"""
Description: Provides NMSim support to soundprophlpr.py

Dependencies: ArcGIS >=10.3, Python 2.7, soundprophlpr.py, nmsimgis_analysis.py
              NumPy

@author Alexander "Sasha" Keyel <skeyel@gmail.com> and Bruce Ikelheimer

AWeight, CWeight, FWeight and some content of other functions by B. Ikelheimer.
Based on code provided by Bruce Ikelheimer and NMSIM Standard Operating
Procedure by Kathryn Nuessly and R post-processing scripts by Damon Joyce

Copyright Information For nmsimhlpr.py:
Copyright (C) 1/22/2016 - 4/12/2017 A.C. Keyel <skeyel@gmail.com>

Alexander "Sasha" Keyel
Postdoctoral Researcher
1474 Campus Delivery
Colorado State University
Fort Collins, CO 80523-1474

   
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License 2.0 as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""

# Import modules & ArcGIS settings
#import multiprocessing # Import multiprocessing package for running nmsimgis
import sys, os, numpy, copy, ctypes, math, subprocess
import cPickle as pickle
from ctypes import cdll
import arcpy, arcpy.sa
from soundprophlpr import seq
from soundprophlpr import calculate_speed_of_sound
arcpy.env.overwriteOutput = True # Make this a setting for all processes
arcpy.CheckOutExtension("Spatial") # Check out the spatial analyst extension
arcpy.env.addOutputsToMap = False  # Prevent intermediates from being added to the map mid-processing #**# Doesn't work


# Core function for running NMSim process
def nmsimgis_setup(model_dir, point_dir, freq_dir, Sound_Source, dem_clip, land_clip, Model_Extent, source_offset,
             receiver_offset, source_info, temp, rh, freq, point_id, point_fill,
             freq_fill, tbx_root, my_times, my_time_labels, keep_intermediates,
             n_processes, is_GUI, model):

    # For troubleshooting
    '''
    elevation = dem_clip
    landcover = land_clip
    head, roll, pitch, vel, engpow, srcfile = source_info
    '''
   
    do_timing = 1
    if my_time_labels == "NA":
        do_timing = 0

    # Record time of nmsimgis start
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Initialization nmsimgis")

    # Set up inputs
    #from nmsimhlpr import * # For testing purposes    
    #freq = freq_lst # For testing purposes
    nmsim_dll = tbx_root + "dlls/NMSim_Libraries.dll"
        
    # Source comes in as a source point from a shapefile
    #Sound_Source = "C:/SPreAD-GIS/csu/data/update_ms/spreadgis/outputs\source_point.shp"

    # Extract elevation from dem.
    xyz_table = point_dir + "xyztable.dbf"
    arcpy.sa.Sample(dem_clip, Sound_Source, xyz_table, "BILINEAR") #Formerly NEAREST

    # Get field name for extracting z value
    dem_name = dem_clip.split("/")[-1]
    dem_name = dem_name[0:10]
    # Check if dem_name has a file extension. If so, drop the file extension from the file name (hope they did not name it with extraneous periods!)
    dem_name = dem_name.split(".")[0] # Split on period & discard everything after the first period.

    # Use a search cursor to extract x,y & z #**# should only be one point - but... you could save time and do this all at once for all points
    with arcpy.da.SearchCursor(xyz_table, ["X", "Y", dem_name]) as srows:
        for row in srows:
            x_source = row[0]
            y_source = row[1]
            z_ground = row[2] #**# This doesn't use bilinear interpolation - this could be a problem.
    
    xyzsrc = [x_source, y_source, z_ground] # Ground positions.
    
    # Extract observer offset from source point, if none, apply default of 1
    #**# Think about using the format used by SPOT & OFFSET from ArcGIS observer points tool
    zsrc = source_offset

    # Get lower left corner from dem for back-conversion (x,y)
    raster_info = arcpy.Describe(dem_clip)
    raster_extent = raster_info.Extent
    #coord_system = raster_info.spatialReference
    cell_x = raster_info.meanCellWidth
    cell_y = raster_info.meanCellHeight
    xmin = raster_extent.XMin
    xmax = raster_extent.XMax
    ymin = raster_extent.YMin
    ymax = raster_extent.YMax

    # Convert landcover to impedance
    impedance = landcover_to_impedance(land_clip, model_dir, Model_Extent)

    dem_correction, impedance_correction = check_rasters(dem_clip, impedance)

    # Bring clipped dem in as a numpy array for further processing
    #**# Consider improving memory management here
    numpy_dem = arcpy.RasterToNumPyArray(dem_clip)
    numpy_impedance = arcpy.RasterToNumPyArray(impedance) # convert impedence to a numpy array as well

    # Correct numpy arrays to force landscape extents to match
    numpy_dem = correct_array(numpy_dem, dem_correction)
    numpy_impedance = correct_array(numpy_impedance, impedance_correction)

    # Check that the correction was applied correctly
    if numpy_dem.shape != numpy_impedance.shape:
        m1 = "Raster arrays must have the same number of rows and columns"
        m2 = " Elevation has %s rows and %s columns, while impedance has %s rows and %s columns"  % (numpy_dem.shape[0], numpy_dem.shape[1], numpy_impedance.shape[0], numpy_impedance.shape[1])
        m3 = " something went wrong with clipping and conversion of %s and %s to numpy arrays."  % (dem_clip, impedance)
        m4 = "Impedance is based on landcover"
        
        raise ValueError("%s%s%s%s" % (m1,m2,m3,m4) )

    # Check that there is enough memory to run the analysis
    n_chunks = check_memory_needs(xmin, xmax, ymin, ymax, cell_x, cell_y, is_GUI)
    arcpy.AddMessage("processing output array in %s pieces" % n_chunks)
    
    # chunk paths will be two nested lists, one for data type, and one for frequencies
    nfreq = len(freq)
    freq_holder = ["NA"] * nfreq

    chunk_paths = [freq_holder] #create a nested list, to support the structure when intermediates are present
    if keep_intermediates == 1:
        chunk_paths = [copy.deepcopy(freq_holder), copy.deepcopy(freq_holder), copy.deepcopy(freq_holder), copy.deepcopy(freq_holder), copy.deepcopy(freq_holder)] # Create a nested list for each intermediate type #freq_holder * 5 did not give an appropriate result

    col_remainder = 0
    row_remainder = 0
    col_subset_start = 0
    row_subset_start = 0

    # Create two sets of environmental settings for processing
    # Save main environmental settings
    main_settings = model_dir + "main_settings.xml"    
    # Delete any pre-existing settings to prevent extent errors
    if os.path.exists(main_settings):
        os.remove(main_settings)
    arcpy.SaveSettings(main_settings)
    # Clear extent setting, to allow for appropriate extent rasters
    arcpy.ClearEnvironment("extent")
    arcpy.env.addOutputsToMap = False  # Prevent intermediates from being added to the map mid-processing
    no_extent_settings = model_dir + "no_extent.xml"
    arcpy.SaveSettings(no_extent_settings)

    for p in xrange(n_chunks):

        # Set up run inputs
        #temp_text = point_dir + "temp_text.txt"
        #args_text = point_dir + "args_text.txt" # To contain arguments with spaces in the name that interfere with the function call
        in_pickle = point_dir + "in_pickle_%s.p" % p
        
        # Set up subprocess.call
        out_pickle = point_dir + "out_pickle_%s.p" % p
        nmsimgis_analysis_py = "%s/scripts/nmsimgis_analysis.py" % tbx_root
        args1 = [p, n_chunks, chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, freq, numpy_dem, numpy_impedance]
        args2 = [source_info, temp, rh, receiver_offset, model_dir, point_dir, dem_name, xyzsrc, zsrc] #, lowerleft] #lowerleft is unpicklable - need to re-create it!
        args3 = [dem_clip, point_id, point_fill, freq_fill, keep_intermediates, n_processes, nmsim_dll, cell_x, cell_y, xmin]
        args4 = [xmax, ymin, ymax, no_extent_settings, do_timing, my_times, my_time_labels, model, out_pickle]
        args = args1 + args2 + args3 + args4
        pickle.dump(args, open(in_pickle, 'wb'))
        
        commands = "python %s %s" % (nmsimgis_analysis_py, in_pickle)
        test = subprocess.Popen(commands, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        output, error_output = test.communicate()
        if error_output != "":
            raise ValueError(error_output)
            
        arcpy.AddMessage("p = %s" % p)

        # Provide a more intuitive error message if this step fails #**# Ideally, you could pipe the error message from subprocess call, or at least write an error log.
            # read in values that were pickled
        pickle_out = pickle.load( open(out_pickle, 'rb') ) # outpickle is defined later & will be defined by this point in the loop.
            #except:
                #raise ValueError("Something went wrong during the analysis (line 183 in nmsimhlpr.py; subprocess.call function)")
            
        chunk_paths = pickle_out[0]
        col_remainder = pickle_out[1]
        row_remainder = pickle_out[2]
        col_subset_start = pickle_out[3]
        row_subset_start = pickle_out[4]
        my_times = pickle_out[5]
        my_time_labels = pickle_out[6]
    

    # if only one run, then no need to remerge
    if p != 0:
        arcpy.LoadSettings(no_extent_settings) # Prevent mosaic raster from applying the environmental settings creating NoData in the output        
        point_lbl = point_id.zfill(point_fill)
        merge_rasters(chunk_paths, keep_intermediates, freq_dir, point_lbl, freq, freq_fill)
        arcpy.LoadSettings(main_settings) # Restore original settings before return
        
    return(my_times, my_time_labels)

# Main nmsimgis analysis, can be split into pieces & looped by memory handler function above.
def nmsimgis_analysis(p, n_chunks, chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, freq, numpy_dem, numpy_impedance, source_info, temp, rh,
                      receiver_offset, model_dir, point_dir, 
                      dem_name, xyzsrc, zsrc, elevation,
                      point_id, point_fill, freq_fill,
                      keep_intermediates, n_processes, nmsim_dll,
                      cell_x, cell_y, xmin, xmax, ymin, ymax, no_extent_settings,
                      do_timing, my_times, my_time_labels, model, ambient_raster = "none"):

    # Create the lower left corner for the main raster
    # note: this variable is not pickleable.
    lowerleft = arcpy.Point(xmin, ymin)

    # get number of frequencies
    nfreq = len(freq)        
    head, roll, pitch, vel, engpow, srcfile = source_info
    
    # Create a numpy array of the same size as numpy_dem
    no_data_value = -9999  # shrunk by two digits to try to keep bit-size down  
    dem_dimensions = numpy_dem.shape
    db_array = numpy.zeros(dem_dimensions) # starting with 0's

    if p == 0:
        arcpy.AddMessage("Input array has dimensions: %s rows and %s columns" % (db_array.shape[0], db_array.shape[1]))
        #sys.stdout.flush()

    # Below, we create to sets of cell-based coordinates. xyzsrc_cell and the original rasters have the original extent
    # db_array may be subset into an array based coordinate system.
    # Both are necessary, as the rasters must retain their original extent, as the source may not fall within one of the 
    # raster subsets. But in order to use the subsets (and have them be smaller), they need their own coordinate system.
    # Apologies for any confusion from the paired coordinate_systems

    # Convert raster source location to a cell based coordinate system (with extent of the original rasters)
    xyzsrc_cell = convert_to_cell_coords(xyzsrc, cell_x, cell_y, xmin, xmax, ymin, ymax, dem_dimensions)    

    # If 0, just run as normal
    if n_chunks == 0:
        lowerleft_adj = lowerleft
        x_offset = 0
        y_offset = 0
    if n_chunks > 0:
        db_array, xy_extent, x_offset, y_offset, lowerleft_adj, col_remainder, row_remainder, col_subset_start, row_subset_start =  db_array_subset(db_array, n_chunks, p, col_remainder, row_remainder, col_subset_start, row_subset_start, dem_dimensions, xmin, xmax, ymin, ymax, cell_x, cell_y, lowerleft)

    print "sub-array dimensions are %s rows and %s columns" % (db_array.shape[0], db_array.shape[1]) #**# Move to top? Or standardize based on an attempt to optimize?

    db_array = db_array + no_data_value #**# Consider whether an alternate no data type would take less space?

    array_lst = []
    # Make a db_array for each frequency
    for i in range(nfreq):
        array_lst.append(copy.deepcopy(db_array)) # deepcopy is needed to avoid having nfreq references to the same list!

    # if creating intermediate files, copy the array list
    #http://stackoverflow.com/questions/2612802/how-to-clone-or-copy-a-list-in-python
    if keep_intermediates == 1:
        source_level_array = copy.deepcopy(array_lst)
        attenuation_array = copy.deepcopy(array_lst)
        absorption_array = copy.deepcopy(array_lst)
        ssl_array = copy.deepcopy(array_lst)

    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Setup for NMSim complete")

    array_dimensions = db_array.shape
    m_max = array_dimensions[1] #3 
    n_max = array_dimensions[0] #3

    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Prior to entering loop over raster")

    # Try once more to get ArcGIS to honor this envionmental setting
    arcpy.env.addOutputsToMap = False  # Prevent intermediates from being added to the map mid-processing

    count = 0
    for m in xrange(0, m_max):
        for n in xrange(0, n_max):
            if count == 0:
                my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Start of cell 1")

            #m = n = 0 # For testing purposes
            # x & y are in the raster coordinates, not the array coordinates
            x = int(m + x_offset) #xmin + m * cell_x + 0.5 * cell_x # Add a half-cell offset to place the point at the midpoint of the cell
            y = int(n + y_offset) #ymax - n * cell_y - 0.5 * cell_y # ditto
            #y = ymin + n * cell_y + 0.5 * cell_y # ditto

            # Get receiver height via bilinear interpolation
            # the dem is oriented differently than is intuitive for me. n is row, m is column
            z = do_bilinear_interpolation(numpy_dem, int(m + x_offset), int(n + y_offset))
        
            # Receiver locations should correspond to the mid-point of every cell with the appropriate elevation
            xyzrec_cell = [x,y,z] # Needs to be in cell units of the rasters. Formerly meters
            zrec = receiver_offset # observer elevation (m AGL)

            if count == 0:
                my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Post-transect; prior to adding task to queue")

            # Compute transect 
            npts, max_dist, dist_vec, hgt_vec, imp_vec, my_times, my_time_labels = get_transect_v2(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, numpy_dem, numpy_impedance, dem_name, point_dir, do_timing, my_times, my_time_labels, count)

            if model == "nmsimgis":
            
                core_out = nmsimgis_core(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, zsrc, zrec, head, roll, pitch,
                                         vel, engpow, srcfile, temp, rh, npts, max_dist,
                                         dist_vec, hgt_vec, imp_vec, freq, nfreq, nmsim_dll, point_dir)

            if model == "iso9613":
                if n == 0 and m == 0:
                    arcpy.AddWarning("Assuming down-wind propagation")    
                    arcpy.AddWarning("Assuming single refraction") 
                    arcpy.AddWarning("Section 7.5 Reflections was omitted from the calculations")
                    arcpy.AddMessage("Annex A (Foliage, Industrial Sites, and Housing) was omitted from the calculations")

                core_out = iso_core(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, zsrc, zrec, head, roll, pitch,
                                         vel, engpow, srcfile, temp, rh, npts, max_dist,
                                         dist_vec, hgt_vec, imp_vec, freq, nfreq, nmsim_dll, point_dir)
                

            db_lst = core_out[0]
            dbband_out = core_out[1]
            attenuation = core_out[2]
            absorption = core_out[3]
            spherical_spreading_loss = core_out[4]
    
            #out_vals = []        
            # extract values out of db_lst
            #for db in db_lst:
            #    out_vals.append(db)
            for i in xrange(len(db_lst)):                
                value = db_lst[i]
                array_lst[i][n, m] = value

                # Adding rounding to prevent a bug where the Sample tool can't extract from 64-bit rasters (try to keep bit size down)
                if keep_intermediates == 1:
                    source_level_array[i][n, m] = round(dbband_out[i],1)
                    attenuation_array[i][n, m] = round(attenuation[i], 1)
                    absorption_array[i][n, m] = round(absorption[i], 1)
                    ssl_array[i][n, m] = round(spherical_spreading_loss, 1)
            
            count += 1
            if count == 1 or count == 10 or count == 100 or count == 1000 or count == 10000 or count == 1000000:
                my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Looped through %s cells" % count)

    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Ended loop over %s cells" % count)
  
    # change environmental settings to avoid extra NoData rows & columns based on original extent size
    arcpy.LoadSettings(no_extent_settings)

    # Get projection from the input elevation file (only input because geoprocessing spatial reference objects can't be passed through the multiprocessing module! i.e. unpickleable)
    proj = arcpy.Describe(elevation).SpatialReference    
    
    # Loop through array list & write to file
    for i in xrange(nfreq):

        # When breaking rasters into pieces, make the first one have the base raster name (it will have the other rasters merged to it in a later step)
        if p == 0:
            p = "" 
        
        this_freq = freq[i]
        point_lbl = str(point_id).zfill(point_fill)
        freq_lbl = str(int(this_freq)).zfill(freq_fill) # int prevents a decimal point, which is then treated as an unknown file extension. So right now, it assumes you won't want to try two frequencies that differ by less than 1
        #big_raster_file = point_dir + "pt%s_pr%s%s_64bit.tif" % (point_lbl, freq_lbl, p) # Why it's saving these in 64 bit format is beyond me!
        raster_file = point_dir + "pt%s_pr%s%s.tif" % (point_lbl, freq_lbl, p) 
        chunk_paths= chunk_paths_update(chunk_paths, 0, i, raster_file) #Add to list for mosaicking
        this_array = array_lst[i]
        this_raster = arcpy.NumPyArrayToRaster(this_array, lowerleft_adj, cell_x, cell_y, no_data_value) # formerly lowerleft
        this_raster.save(raster_file)
        arcpy.DefineProjection_management(raster_file, proj)
        
        # Convert from 64-bit raster to smaller pixel depth to avoid a bug with the sample tool (which is fixed in 10.5)
        #arcpy.CopyRaster_management(big_raster_file, raster_file, "","", no_data_value, "","", "32_BIT_FLOAT")        
        #arcpy.Delete_management (big_raster_file)
        
             
        if keep_intermediates == 1:
            source_level_file = point_dir + "pt%s_sr%s%s.tif" % (point_lbl, freq_lbl, p)
            attenuation_file = point_dir + "pt%s_at%s%s.tif" % (point_lbl, freq_lbl, p)
            absorption_file = point_dir +  "pt%s_ab%s%s.tif" % (point_lbl, freq_lbl, p)
            ssl_file = point_dir + "pt%s_ss%s%s.tif" % (point_lbl, freq_lbl, p)
            chunk_paths = chunk_paths_update(chunk_paths, 1, i, source_level_file)
            chunk_paths = chunk_paths_update(chunk_paths, 2, i, attenuation_file)
            chunk_paths = chunk_paths_update(chunk_paths, 3, i, absorption_file)
            chunk_paths = chunk_paths_update(chunk_paths, 4, i, ssl_file)
            
            source_level_raster = arcpy.NumPyArrayToRaster(source_level_array[i], lowerleft_adj, cell_x, cell_y, no_data_value)
            source_level_raster.save(source_level_file)
            arcpy.DefineProjection_management(source_level_file, proj)
            attenuation_raster =  arcpy.NumPyArrayToRaster(attenuation_array[i], lowerleft_adj, cell_x, cell_y, no_data_value)
            attenuation_raster.save(attenuation_file)
            arcpy.DefineProjection_management(attenuation_file, proj)
            absorption_raster =   arcpy.NumPyArrayToRaster(absorption_array[i], lowerleft_adj, cell_x, cell_y, no_data_value)
            absorption_raster.save(absorption_file)
            arcpy.DefineProjection_management(absorption_file, proj)
            ssl_raster = arcpy.NumPyArrayToRaster(ssl_array[i], lowerleft_adj, cell_x, cell_y, no_data_value)
            ssl_raster.save(ssl_file)
            arcpy.DefineProjection_management(ssl_file, proj)
        
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Data written to file")
    
    return(chunk_paths, col_remainder, row_remainder, col_subset_start, row_subset_start, my_times, my_time_labels)


# simple function to handle updating of the chunk_paths variable
def chunk_paths_update(chunk_paths, index, i, new_value):
    old_list = chunk_paths[index][i]
    if old_list == "NA":
        chunk_paths[index][i] = [new_value]
    else:
        old_list.append(new_value)
        chunk_paths[index][i] == old_list
    
    return chunk_paths

# Function to determine a good number to use for running the landscape
# This function does not assess the computer's memory available for the task
# That will be a future code upgrade.
def check_memory_needs(xmin, xmax, ymin, ymax, cell_x, cell_y, is_GUI):

    #**# Limit to 2 pieces, as GUI throws an error in the mosaic raster step if
    # there are more than 2 pieces (can bypass this using copy raster, but
    # that will make the code messy and more computationally intensive for
    # python users)
    if is_GUI == "YES":
        n_chunks = 2
    else:
        rows = (ymax - ymin) / float(cell_y)
        cols = (xmax - xmin) / float(cell_x)
        
        ncells = rows * cols
        processing_window = 10000 #**# Change to an input? Then the user can 
        # adjust based on their computer's memory.
        
        n_chunks = int(math.ceil(ncells / processing_window))
    
        # Set minimum number of chunks due to an unfixed bug when only one piece is specified    
        if n_chunks < 2:
            n_chunks = 2
    
    return n_chunks

# Subset the area to be run
# Changing xmin & ymin need to be compensated for in the main raster through x & y offsets
def db_array_subset(db_array, n_chunks, p, col_remainder, row_remainder, col_subset_start, row_subset_start, dem_dimensions, xmin, xmax, ymin, ymax, cell_x, cell_y, lowerleft):

    # Get dimensions
    ncol = dem_dimensions[1]
    nrow = dem_dimensions[0]

    # Convet x & y to cell coords
    ymax_cell = 0
    ymin_cell = nrow - 1 # ymax starts at 0, not 1, so subtract 1 to compensate

    # subset first by columns
    col_chunks = n_chunks
    row_chunks = 0        
    if n_chunks > ncol:
        arcpy.Message("Your landscape is being split by columns and rows to make it more manageable for your computer's memory. Unfortunately, the code that does splitting by rows has not been tested, and therefore likely contains a bug. Please report any odd or irregular results or error messages")
        col_chunks = ncol
        row_chunks = n_chunks - ncol

    db_array_subset, new_xmin, new_xmax, x_offset, col_remainder, col_subset_start = subset_array(db_array, p, col_chunks, ncol, col_remainder, col_subset_start, 1) # 1 is for columns
    # Watch for if this misses any columns!    
    #print "col_subset_start: %s; xmin: %s; x_offset: %s; cell_x: %s" % (col_subset_start, xmin, x_offset, cell_x)
    new_x_UTM = xmin + x_offset * cell_x
    #print "new_x_UTM: %s" % new_x_UTM

    # then by rows if necessary  
    new_ymin = ymin_cell # default of no change (i.e. only columns split into pieces)
    new_ymax = ymax_cell 
    y_offset = 0 
    new_y_UTM = ymin
    if row_chunks > 0:
        db_array_subset, new_ymax, new_ymin, y_offset, row_remainder, row_subset_start = subset_array(db_array_subset, p, row_chunks, nrow, row_remainder, row_subset_start, 0) # axis = 0 is for rows
        # Note reversal in ymin & ymax. The start is the top, and the finish is the bottom, because numpy arrays start with the top, not the bottom.
        new_y_UTM = ymin + y_offset * cell_y

    #Finally, update extents (xy_extent now contains the new_min & new_max. Offset allows for correct selection from the original rasters)
    xy_extent = [new_xmin, new_xmax, new_ymin, new_ymax]    
    lowerleft_adj = arcpy.Point(new_x_UTM, new_y_UTM)
    
    return(db_array_subset, xy_extent, x_offset, y_offset, lowerleft_adj, col_remainder, row_remainder, col_subset_start, row_subset_start)

# Code to do subsetting
def subset_array(array, p, chunks, dim, remainder, subset_start, axis):
    import math
    
    # Calculate the size of the original raster that needs to be subset
    subset_size_initial = float(dim) / chunks # float to ensure that fractions are preserved
    subset_size = int(math.ceil(subset_size_initial))

    # calculate remainder
    remainder = (subset_size_initial - subset_size) + remainder    
    #remainder = -(dim % chunks) / float(chunks) + remainder #Get the difference between ceil and what was actually run. % gives you the remainder, chunks gives you the denominator

    # Apply remainder any time it gets to be greater than 1
    if remainder <= -1: # remember, remainder is negative
        subset_size = int(subset_size + math.ceil(remainder)) # ceil here - want -1.4 to round to -1 # int is to avoid a deprecation warning where the integer value in float form will not work in future versions
        remainder = remainder - math.ceil(remainder) # - here because subtracting a negative is addition e.g., -1.25 - -1 = -0.25

    # Calculate start and end points for subsetting (and for raster dimensions)
    offset = subset_start #p * subset_size # The old assumed equal sized subsets!
    new_min = 0 + offset
    new_max = new_min + subset_size
    subset_start = new_max

    #print "offset is %s" % offset

    if axis == 1:
        array_subset = array[0:, new_min:new_max] #0: keeps all rows
    
    if axis == 0:
        array_subset = array[new_min:new_max, 0:] # Here 0: keeps all columns

    return(array_subset, new_min, new_max, offset, remainder, subset_start)

# Code to remerge the different pieces of the output raster
#**# When running nmsimgis directly from the ArcGIS toolbox, more than 2 pieces
# causes the Mosaic Raster step to fail (unable to delete raster!). This can
# be worked around by creating a copy and using the copy in subsequent steps
# but that could add computational overhead. Right now, the GUI version is limited
# to 2 pieces, but this will limit the applications of the GUI interface.
def merge_rasters(chunk_paths, keep_intermediates, freq_dir, point_lbl, freq, freq_fill):

    # This loops through all intermediates if present
    for output_type in chunk_paths:
        # Loop through each frequency
        for frequency in output_type:
            
            first_raster = frequency[0]
            #count = 0
            for path in frequency[1:]:
                #count += 1
                #arcpy.AddMessage(first_raster)
                #arcpy.AddMessage(path)
                # Merge rasters together
                # This will make the first raster contain all the data
                
                # Patch to try to fix GUI problem where more than 3 pieces can't be combined
                #new_first_raster = first_raster[:-4] + "p" + str(count) + ".tif"
                #arcpy.CopyRaster_management(first_raster, new_first_raster)
                
                #rcpy.Mosaic_management([path], new_first_raster, "LAST") # here, raster values should be equivalent, using LAST to try to get around nodata problem. Later, some version of addition could be used, except the default methods don't allow dB addition.
                arcpy.Mosaic_management([path], first_raster, "LAST") # here, raster values should be equivalent, using LAST to try to get around nodata problem. Later, some version of addition could be used, except the default methods don't allow dB addition.

    #**# This section needs patching if going the new_first_raster route
    for i in xrange(len(chunk_paths[0])):
        frequency = chunk_paths[0][i]
        freq_lbl = str(int(freq[i])).zfill(freq_fill)
        Raster1a = frequency[0]
        Raster1b = freq_dir + "pt" + point_lbl + "_pr" + freq_lbl + ".tif"
        #print freq_fill        
        #print Raster1a
        #print Raster1b
        #sys.stdout.flush()
        # Change to 32 bit float shrinks the raster & patches around a bug in the sample tool.
        arcpy.CopyRaster_management(Raster1a, Raster1b, "","", -9999, "","", "32_BIT_FLOAT")

    arcpy.AddMessage("Finished Mosaic")

# Define function to run all NMSim functions
def nmsimgis_core(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, zsrc, zrec, head, roll, pitch, vel, engpow,
                  srcfile, temp, rh, npts, max_dist, dist_vec, hgt_vec,
                  imp_vec, freq, nfreq, nmsim_dll, point_dir, ambient_vals = ["none"]):

    # for testing purposes
    #testing = 1
    #temp = temp_s
    #rh = hum_s
    #xyz_rec = xyzrec_cell
    #srcfile

    import math
    temp = float(temp) # needs to be in float format
    rh = float(rh) # Needs to be in float format
    freq = map(float, freq) # Convert list to float format

    #Reconvert xyzsrc and xyzrec from cell units to m
    xsrc_m = xyzsrc_cell[0] * cell_x
    ysrc_m = xyzsrc_cell[1] * cell_y
    xyzsrc = [xsrc_m, ysrc_m, xyzsrc_cell[2]] # z coordinate was ok
    
    xrec_m = xyzrec_cell[0] * cell_x
    yrec_m = xyzrec_cell[1] * cell_y
    xyzrec = [xrec_m, yrec_m, xyzrec_cell[2]] # z coordinate was ok    

    # Calculate 3D distance between source and receiver
    dist_3d = math.sqrt((xsrc_m - xrec_m)**2 + (ysrc_m - yrec_m)**2 + (xyzsrc[2] + zsrc - xyzrec[2] - zrec)**2)
       
    # Set up source information for later propagation analysis
    refdist_out, distance_out, dbband_out = NMSim_Source(nmsim_dll, xyzsrc, xyzrec, zsrc,
                                           zrec, head, roll, vel, pitch, engpow,
                                           srcfile, nfreq)
   
    # Evaluate atmospheric absorption for a single receiver point
    absorption = NMSim_atmospheric_absorption(nmsim_dll, temp, rh, freq,
                                                            nfreq, dist_3d)
    
    # Evaluate sound propagation based on terrain, impedance, source height,
    # and receiver height
    attenuation = NMSim_propagation(nmsim_dll, zsrc, zrec, npts,
                                              dist_vec, hgt_vec, imp_vec, freq,
                                              nfreq)
    
    # Add spherical spreading loss (this is even across frequencies)    
    if distance_out == 0:
        distance_out = 1 # Switch to 1 m distance if the source point is exactly on the midpoint of the cell of origin
    spherical_spreading_loss = 20 * math.log10(distance_out * refdist_out**-1) #do = source reference distance (e.g., spl measured at 15m), d = distance between source & receiver
    
    
    # Combine frequency loss to get final decibel values (convert to function)
    # I'm sure there are better ways to code this
    db_lst = [0] * nfreq
    for i in range(nfreq):
        this_source_level = dbband_out[i]
        this_atmoabs = absorption[i]
        this_attenuation = attenuation[i]
    
        db_out = this_source_level - this_atmoabs - this_attenuation - spherical_spreading_loss

        # Cap out extreme negative values to avoid code crashes
        if db_out < -300:
            db_out = -300
        
        db_lst[i] = db_out

    # convert negative values to 0's
    #db_lst = remove_negatives(db_lst) # Not really appropriate to do this, especially before summing up for overall results.

    # Output results of sound propagation for use by end users
    return(db_lst, dbband_out, attenuation, absorption, spherical_spreading_loss)

# Convert from UTM to cell coordinates
def convert_to_cell_coords(xyzsrc, cell_x, cell_y, xmin, xmax, ymin, ymax, dimensions):
    
    '''
    test code:
    xmin = 1000
    xmax = 2000
    ymin = 4000
    ymax = 8000
    cell_x = 100
    cell_y = 400
    xyzsrc = [1550, 4900, 20] 
    '''
    
    # Start with xyzsrc in UTM coordinates
    x = xyzsrc[0]    
    y = xyzsrc[1]
    z = xyzsrc[2] # z will remain unchanged

    # Get number of cells and check against original raster to confirm no problems
    n_x_cells = (xmax - xmin) * cell_x**-1
    n_x_cells = int( round(n_x_cells, 0))
    n_y_cells = (ymax - ymin) * cell_y**-1 
    n_y_cells = int( round(n_y_cells, 0))
    
    
    # x corresponds to columns (m), y corresponds to rows (n) (except in reverse order, as y decreases as one goes down, and rows increase)    
    if n_x_cells != dimensions[1]:
        raise ValueError("Calculated number of cells in x direction (columns) (%s) does not match actual number of cells(%s)" % (n_x_cells, dimensions[1]))
    if n_y_cells != dimensions[0]:
        raise ValueError("Calculated number of cells in y direction (rows) (%s) does not match actual number of cells(%s)" % (n_y_cells, dimensions[0]))
              
    # Rescale to cell units
    new_x = cell_rescale_x(x, xmin, cell_x)
    new_y = cell_rescale_y(y, ymin, cell_y, n_y_cells)
    
    # Check that new_x, new_y are on the grid
    check_x_y(new_x,new_y, n_x_cells, n_y_cells)

    # Recompile cell
    xyzsrc_cell = [new_x, new_y, z]
    
    return(xyzsrc_cell)

# Actually do rescaling
def cell_rescale_x(x, xmin, cell_x):
    half_x = 0.5 * cell_x
    # Subtract minimum (which is adjusted to give cell midpoints)
    new_x = x - (xmin + half_x)
    # Divide by cell size
    new_x = new_x * cell_x**-1
    # Index begins with 0
    return(new_x)

# Actually do rescaling. y cell coordinates decrease from maximum
def cell_rescale_y(y, ymin, cell_y, n_y_cells):
    half_y = 0.5 * cell_y
    # Subtract minimum (which is adjusted to give cell midpoints)
    new_y = y - (ymin + half_y)
    # Divide by cell size
    new_y = new_y * cell_y**-1
    # Index begins with 0
    
    # subtract new_y from maximum number of cells to account for the fact that
    # the row numbering starts at top & goes down, while y traditionally starts at bottom and increases going up!
    new_y = n_y_cells - new_y - 1 # -1 adjusts for fact that index starts at 0, and not 1
    
    return(new_y)

# Check that source x & y are within the landscape boundaries
def check_x_y(new_x,new_y, n_x_cells, n_y_cells):
    is_error = 0    
    xproblem = yproblem = "ok"
    if new_x < 0:
        is_error = 1
        xproblem = "too far to the left (west)"
    if new_y < 0:
        is_error = 1
        yproblem = "too far to the top (north)"
    if new_x > n_x_cells:
        is_error = 1
        xproblem = "too far to the right (east)"
    if new_y > n_y_cells:
        is_error = 1
        yproblem = "too far to the bottom (south)"
    if is_error == 1:
        raise ValueError("Source points must be more than 1/2 cell length/width from the edge of the landscape. Source x-coordinate is %s and source y-coordinate is %s" % (xproblem, yproblem))

def check_rasters(elevation, impedance):

    elev = arcpy.Describe(elevation)
    imp = arcpy.Describe(impedance)
    
    # Check that cell sizes are the same
    if round(elev.meanCellWidth,10) != round(imp.meanCellWidth, 10): # round was necessary because ArcGIS wouldn't give exact results even when snap raster and cell size environment settings were set.
        raise ValueError("elevation and impedance (based on landcover) cell widths do not match (%s, %s" % (elev.meanCellWidth, imp.meanCellWidth))
    if round(elev.meanCellHeight, 10) != round(imp.meanCellHeight, 10):
        raise ValueError("elevation and impedance (based on landcover) cell heights do not match (%s, %s" % (elev.meanCellHeight, imp.meanCellHeight))

    # Identify whether to remove the top or bottom row from the chosen array
    # Put in context of cell size
    XMin_diff = (elev.Extent.XMin - imp.Extent.XMin) / float(elev.meanCellWidth)
    XMax_diff = (imp.Extent.XMax - elev.Extent.XMax) / float(elev.meanCellWidth)
    YMin_diff = (elev.Extent.YMin - imp.Extent.YMin) / float(elev.meanCellHeight)
    YMax_diff = (imp.Extent.YMax - elev.Extent.YMax) / float(elev.meanCellHeight)

    # Calculate offsets    
    dem_drop_left, imp_drop_left = comparison(XMin_diff)  
    dem_drop_right, imp_drop_right = comparison(XMax_diff)  
    dem_drop_bottom, imp_drop_bottom = comparison(YMin_diff)  
    dem_drop_top, imp_drop_top = comparison(YMax_diff)  

    '''
    # This is flawed & irrelevant
    # Check for rounding problems
    dem_drop_left, dem_drop_right = round_check(dem_drop_left, dem_drop_right)
    dem_drop_bottom, dem_drop_top = round_check(dem_drop_bottom, dem_drop_top)    
    imp_drop_left, imp_drop_right = round_check(imp_drop_left, imp_drop_right)
    imp_drop_bottom, imp_drop_top = round_check(imp_drop_bottom, imp_drop_top)    
    '''            
            
    dem_correction = [dem_drop_left, dem_drop_right, dem_drop_bottom, dem_drop_top]
    impedance_correction = [imp_drop_left, imp_drop_right, imp_drop_bottom, imp_drop_top]
    return(dem_correction, impedance_correction)

# Simplify the comparsion for getting drops out
def comparison(diff):    
    if diff < 0:
        dem_drop = round(abs(diff), 0)
        imp_drop = 0
    if diff > 0:
        dem_drop = 0
        imp_drop = round(diff, 0)
    if diff == 0:
        dem_drop = 0
        imp_drop = 0

    return(dem_drop, imp_drop)


# Correct arrays for differences in extent size (because the ArcGIS tools that should do this don't completely work!)
def correct_array(array, correction_factor):
    left_correction = int(correction_factor[0])
    right_correction = int(correction_factor[1])
    bottom_correction = int(correction_factor[2])
    top_correction = int(correction_factor[3])
    
    # Give syntax for rows & columns of interest out appropriate syntax to perform drops
    # array[:,0] # first column
    # array[:,-1] # Last column
    # array[0,:] # First row
    # array[-1,:] # Last row
    
    # Apply correction
    array = array[:,left_correction:] # if left correction is 0, no change is made. If it is 1, the first column is dropped, etc. Trailing colon says to use the entire array
    # :0 will remove all rows/columns! Since no change is desired, skip the processing step if 0 where this syntax is used!
    if right_correction != 0:
        array = array[:, :(-right_correction)] # if right correction is 0, no change is made. If it is 1, the value is -2, and the last column is dropped
    if bottom_correction != 0:
        array = array[:(-bottom_correction), :] # if bottom correction is 0, no change. If it is 1, the last row is dropped, etc.
    array = array[top_correction:,:] # if top correction is 0, no change, if 1, the first row is dropped, etc.

    return(array)

# Function to get the raster name used in field names by arcgis
def get_raster_name(raster): 
    raster_name = raster.split("/")[-1] #**# May cause an error if someone uses "\\" in their paths
    raster_name = raster_name[0:10]
    # Strip any file endings from the name
    raster_name = raster_name.split(".")[0]
    return(raster_name)


# Sample a raster & extract the results #This code is similar to that used in validation script
def sample_extract(raster, points, point_dir, my_times, my_time_labels, count):
    import arcpy, arcpy.sa
    arcpy.CheckOutExtension("Spatial")
    
    raster_name = get_raster_name(raster)
    
    dbf = point_dir + "temp%s.dbf" % count #.dbf
    if count == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "Pre-sampling %s" % raster_name)

    arcpy.sa.Sample(raster, points, dbf, "BILINEAR", "FID","") # Changed to bilinear, because Bruce said that that is what is used by NMSIM

    if count == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "Post-sampling %s" % raster_name)

    # use a cursor to extract values from dbf for hgt_vec
    #**# OUTDATED SEARCH CURSOR
    out_vec = []
    trows = arcpy.SearchCursor(dbf)
    for trow in trows:
        value = trow.getValue(raster_name)
        out_vec.append(value)
    del trows, trow

    if count == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "Post-search cursor %s" % raster_name)
    
    return(out_vec, my_times, my_time_labels)

# Values based on GroundImpedance.csv from Bruce Ikelheimer, which is used in NMSim
# Now no longer following the NLCD system but a more generalized system
# Default system
    #1 = WATER = water
    #2 = SNOW = Perennial Ice/Snow
    #3 = URBAN1 = developed, open space
    #4 = URBAN2 = developed, low intensity
    #5 = URBAN3 = developed, medium intensity
    #6 = URBAN4 = developed, high intensity
    #7 = BARREN = Barren
    #8 = SHORE = Unconsolidated Shore
    #9 = FOREST = Forest (all types)
    #10 = SHRUB or GRASS = grassland/shrubland/pastureland/cropland (all types)
    #11 = WETLAND = wetlands (all types)
    #12 = LICHEN = lichens (alaska only)
    #13 = MOSS = Moss (alaska only)
    #14 = AK = No Data (Alaska only)
def landcover_to_impedance(landcover, point_dir, Model_Extent):
    impedance_pre = point_dir + "imp_pre" #Clip should not be necessary, but ArcGIS ADDS A NON-EXISTENT ROW!!!
    impedance = point_dir + "impedance"

    # Delete any pre-existing layers (overwrite should take care of this, but the GUI is ignoring that setting)
    arcpy.Delete_management(impedance_pre)
    arcpy.Delete_management(impedance)

    water = "WATER 100000"    
    other = "AK 200; SNOW 50; SHORE 30; LICHEN 10000; MOSS 70"
    wetland = "WETLAND 100000"
    forest = "FOREST 70"
    shrub = "SHRUB 200; GRASS_HERB 200"
    barren = "BARREN 100000"
    urban = "URBAN1 200;URBAN2 300;URBAN3 5000;URBAN4 10000"
    custom1 = "C2 2; C4 4; C8 8; C16 16; C32 32; C64 64; C128 128; C256 256; C512 512; C1024 1024"
    custom2 = "C2048 2048; C4096 4096; C8192 8192; C16384 16384; C32768 32768; C65536 65536; C131072 131072"
    reclass_string = "%s;%s;%s;%s;%s;%s;%s;%s;%s" % (other, forest, shrub, barren, water, urban, wetland, custom1, custom2)
        
    reclass = arcpy.sa.Reclassify(landcover, "NMSIMTYPE", reclass_string, "DATA")

    reclass.save(impedance_pre)
    
    #**# Patched more cleverly above. Hopefully that works.
    # This try/except is here because the GUI version was not honoring the overwrite environmental setting.
    #try:        
    #    reclass.save(impedance_pre)
    #except:
    #    impedance_pre = point_dir + "imp_pre2"
    #    reclass.save(impedance_pre)
    #    # Delete the impedance file, as this will crash the GUI also
    #    arcpy.Delete_management(impedance)
    arcpy.Clip_management(impedance_pre, Model_Extent, impedance, "", "", "", "NO_MAINTAIN_EXTENT")
    
    return impedance



# Values based on GroundImpedance.csv from Bruce Ikelheimer, which is used in NMSim
#FORMERLY: Values based on ImpAttrib.csv from Bruce Ikelheimer, which is [NOT] used in NMSim 
def landcover_to_impedance_old(landcover, point_dir, Model_Extent):
    import arcpy   
    impedance_pre = point_dir + "imp_pre" #Clip should not be necessary, but ArcGIS ADDS A NON-EXISTENT ROW!!!
    impedance = point_dir + "impedance"
    other = "1 200; 12 50; 32 30; 73 10000; 74 70"
    wetland1 = "90 100000; 91 100000; 92 100000; 93 100000; 94 100000"
    #wetland2 = "95 100000; 96 100000; 97 100000; 98 100000; 99 100000" 
    forest = "42 70;43 70;41 70"
    shrub = "51 200; 52 200"
    grass = "71 200; 72 200; 81 200; 82 200; 95 40"
    barren = "31 100000"
    water = "11 100000"
    urban = "21 200;22 300;23 5000;24 10000"    
    reclass_string = "%s;%s;%s;%s;%s;%s;%s;%s" % (other, forest, shrub, grass, barren, water, urban, wetland1) #, wetland2;%s
    
    #**# For some reason I'm getting an error with the wetland2 values
    # arcgisscripting.ExecuteError: ERROR 000622: Failed to execute (Reclassify). Parameters are not valid.
    # ERROR 000628: Cannot set input into parameter remap.
    
    reclass = arcpy.sa.Reclassify(landcover, "VALUE", reclass_string, "DATA")
    reclass.save(impedance_pre)
    arcpy.Clip_management(impedance_pre, Model_Extent, impedance, "", "", "", "NO_MAINTAIN_EXTENT")
    
    return impedance

# Convert landcover values to impedances
def landcover_to_impedance_oldest(lc_vec):
       

    imp_vec = []
    
    for i in range(len(lc_vec)):
        landcover = str(int(lc_vec[i])) #Why I wrote these as strings, I do not know.
        
        if landcover == "42" or landcover == "43" or landcover == "41":        
            imp_value = 60 #CON or HWD
        elif landcover == "52":
            imp_value = 50 #SHB, Taken from orchard/hedgerows
        elif landcover == "71" or landcover == "95":
            imp_value = 40 #HEB, from Vegetation (General)
        elif landcover == "31":
            imp_value = 650 #BAR, for Soil & Ground surface
        elif landcover == "11":
            imp_value = 1000000 #WAT, from Fresh Water (General)
        elif landcover == "21" or landcover == "22" or landcover == "23" or landcover == "24":
            imp_value = 30000 #URB, value is for Aircraft Parking Area, Open Storage, Artificial Mountain, were closest thing to Urban I could find! 30000 is also true for Rocky, Rough Surface
        else:
            raise ValueError("A subset of NLCD 2011 codes are used by nmsimgis\n",
                             "Landcover type %s is not recognized/supported by the code." % landcover)
           
        '''
        Old version based on SPREADTYPE. BUT SAMPLE TAKES FROM THE VALUE FIELD
        so this version is based on the NLCD fields. We'll need to figure out
        a better way to get this working.
        if landcover == "CON":
            imp_value = 60
        if landcover == "HWD":
            imp_value = 60
        if landcover == "SHB": #Taken from orchard/hedgerows
            imp_value = 50
        if landcover == "HEB":
            imp_value = 40 # From Vegetation (General)
        if landcover == "BAR":
            imp_value = 650 # Soil & Ground surface
        if landcover == "WAT":
            imp_value = 1000000 # Fresh Water (General)
        if landcover == "URB":
            imp_value = 30000 # Aircraft Parking Area, Open Storage, Artificial Mountain, were closest thing to Urban I could find! 30000 is also true for Rocky, Rough Surface
        '''
        imp_vec.append(imp_value)
        
    return imp_vec

# Get a transect for calculating NMSim results. Previous code works, but is slow (0.7 sec). The Aim is to increase performance to tolerable levels! Ideal is to execute in <0.1 s
# cut , cell_x, cell_y, coord_system
def get_transect_v2(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, numpy_dem, numpy_impedance, dem_name, point_dir, do_timing, my_times, my_time_labels, count):
    #import time
    #this_time = time.time()
    import math

    if count == 0:
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "update this timing code or disable")

    # calculate max_dist as distance between xyzsrc & xyzrec (without respect to elevation)
    xsource = xyzsrc_cell[0]
    ysource = xyzsrc_cell[1]
    xreceiver = xyzrec_cell[0]
    yreceiver = xyzrec_cell[1]
    xdiff = (xreceiver - xsource) * cell_x # Convert back to m
    ydiff = (yreceiver - ysource) * cell_y # convert back to m
    max_dist = math.sqrt( (xdiff)**2 + (ydiff)**2 ) #**# can add z if needed

    # Identify transect points
    npts = 100  # Bruce Ikelheimer said that ONECUT algorithm is ALWAYS given 100 points.
    dist_vec = find_distances(max_dist, npts)

    # Get coordiantes for dist_vec locations
    x_coords = find_coords(xsource, xreceiver, npts)
    y_coords = find_coords(ysource, yreceiver, npts)
    
    
    # use coordinates to extract elevation
    hgt_vec = extract_values(numpy_dem, x_coords, y_coords)

    # use coordiantes to extract impedance
    imp_vec = extract_values(numpy_impedance, x_coords, y_coords)

    # Check that everything is the correct length
    if len(dist_vec) != npts or len(hgt_vec) != npts or len(imp_vec) != npts:
        raise ValueError("npts (%s), dist_vec (%s), hgt_vec (%s), and imp_vec (%s) all must be the same length." % (npts, dist_vec, hgt_vec, imp_vec))
    
    '''
    #**# Temporary code to troubleshoot problem with nmsimgis results
    diagnostic_file = "C:/smt/csu/supporting_docs/nmsim_vs_nmsimgis/diagnostics_sh.csv"
    write_header = 0    
    import os
    if not os.path.exists(diagnostic_file):
        write_header = 1

    with open(diagnostic_file, 'a') as f:
        # write header
        if write_header == 1:
            f.write("SourceX,SourceY,SourceZ,ReceiverX,ReceiverY,ReceiverZ,Distance,Elevation,Impedance\n")
        # write values
        for i in xrange(len(dist_vec)):
            f.write("%s,%s,%s,%s,%s\n" % (str(xyzsrc_cell).strip("[]"), str(xyzrec_cell).strip("[]"), dist_vec[i],hgt_vec[i],imp_vec[i]))
    '''
    #arcpy.AddMessage("Elapsed time: %s" % (time.time() - this_time))
    

    return(npts, max_dist, dist_vec, hgt_vec, imp_vec, my_times, my_time_labels)

def find_distances(max_dist, npts, n_digits = 8): #offset = 0, 
    increment = max_dist * (npts - 1)**-1 # - 1 accounts for the inclusion of the end points, **-1 does division
    increment = round(increment, n_digits) # Add rounding to ensure a uniform increment

    # If the increment is 0, then just create a vector of 0's
    #if round(increment, n_digits) == 0:
    if increment == 0:
        dist_vec = [0] * npts

    else:
        # create dist_vec
        dist_vec = seq(0,max_dist, increment)
        # No longer necessary with shift in rounding step
        #for i in range(len(dist_vec)):
        #    dist_vec[i] = round(dist_vec[i], n_digits)
    
        # check that the final point was added to dist_vec (if dist_vec is one short, add the final point)
        # Sometimes there is a mismatch due to rounding, and the last point is not added.
        if len(dist_vec) == (npts - 1):
            dist_vec.append(max_dist)

            
    return(dist_vec)

# Extract coordinates for a single direction
def find_coords(source, receiver, npts):
    direction = 1
    if source > receiver:
        direction = -1
    
    abs_diff = abs(receiver - source)
    # Get distances between source & receiver in this dimension    
    coords_1d = find_distances(abs_diff, npts)  # 0 corresponds to the offset input.

    # Convert distances to coordinates        
    for i in range(len(coords_1d)):
        this_dist = coords_1d[i]
        # Start at source location, then add distance in the appropriate direction
        new_dist = source + this_dist * direction
        coords_1d[i] = new_dist

    return(coords_1d)

# Extract values from a raster
def extract_values(in_array, x_coords, y_coords, extract_type = "BILINEAR"): #BILINEAR
    out_vec = []
    # x & y coords have the same length
    for i in range(len(x_coords)):
               
        if extract_type == "NEAREST":
            raise ValueError("There is an undiagnosed problem where the nearest neighbor calculation fails (e.g., from Source 1,8 to Receiver 0,8 in the designed landscape)")
            # Get value out of raster
            x = int(x_coords[i])
            y = int(y_coords[i])
            this_value = in_array[y, x] # array is rows, columns. rows are y (but inverted), columns are x

        if extract_type == "BILINEAR":
            xcoord = round(x_coords[i], 9) # Rounding is to try to correct a problem with floating points not being exact
            ycoord = round(y_coords[i], 9)
            this_value = do_bilinear_interpolation(in_array, xcoord, ycoord)

        # Append the value to the list
        out_vec.append(this_value)
    
    return(out_vec)

# Bilinear interopolation as described by Bruce Ikelheimer
# As for the bilinear interpolation, the formula assumes rectangular control points. X is the position along the x axis within a cell, and dx is the size of each cell. Therefore  x=0 is on the left hand side of the cell, and x = dx is on the right hand side. The same goes for Y. The cell is defined by an index, with i being the left side, i+1 being the right side, j being the bottom, and j+1 being the top. For this, the elevation data is assumed to be stored in a 2D array called zgrid.
#     fracx = x/dx
#     fracy = y/dy
#
#     z11=zgrid(i,j)
#     z21=zgrid(i+1,j)
#     z12=zgrid(i,j+1)
#     z22=zgrid(i+1,j+1)

# The following is a two-way linear interpolation formula:
#     zz = z11 + fracx*(z21-z11) + fracy*(z12-z11)
#                 + fracx*fracy*(z11+z22-z21-z12) 
# NOTE: This uses an inverted coordinate system for y (relative to the above), because of the row/column orientation
def do_bilinear_interpolation(in_array, xcoord, ycoord):
    import math
    # Get coordinates of upper left corner
    j = int(math.floor(xcoord))
    i = int(math.floor(ycoord))

    # dx & dy are set to 1, coordinate system is in cell units.
    fracx = round((xcoord - j), 8) # get distance past cell index. Rounding is to prevent minor floating point irregularities from crashing the programs
    fracy = round((ycoord - i), 8) # get distance past the cell index
    
    # assign corner values, check that i + 1 & j + 1 are in the array, if not, handle appropriately
    z11, z21, z12, z22 = check_i_j(in_array, i, j, fracx, fracy, xcoord, ycoord)
    
    zz = z11 + fracy*(z21-z11) + fracx*(z12-z11) + fracx*fracy*(z11+z22-z21-z12) # fracx & fracy switched from formula above because x = j and y = i, due to the row, column orientation of the matrix
    
    return(zz)


# Check that if a source or receiver falls on the endpoint, the correct action is taken
def check_i_j(in_array, i, j, fracx, fracy, xcoord, ycoord):
    # Added rounding to 4 digits, because I'm getting an error about overflow in short scalars! Which shouldn't be happening, but it is.

    is_error = 0

    dimensions = in_array.shape
    x_dim = dimensions[1]
    y_dim = dimensions[0]

    z11 = in_array[i,j]
    z11 = round(z11, 8)

    # My array is indexed by rows, columns. Rows increase as one goes down. But that inversion really shouldn't change things much.
    # - 1 is because python indices start with 0, so if the number of cells is reached, the index is already out of range.
    if j < (x_dim - 1):
        z12 = in_array[i, j + 1]
        z12 = round(z12, 8)
    if i < (y_dim - 1):
        z21 = in_array[i + 1, j]
        z21 = round(z21, 8)
    if j < (x_dim - 1) and i < (y_dim - 1):
        z22 = in_array[i + 1, j + 1]
        z22 = round(z22, 8)
    
    # Check if array indices are at edge of range
    if (x_dim - 1) == j:
        # if yes, check if fractions is 0. this is ok, but requires different calculations
        if fracx == 0:
            z12 = 0 # won't matter, because fracx will nullify it.
            z22 = 0 # won't matter, because fracx will make it 0
        else:
            is_error = 1
    if (y_dim - 1) == i:
        if fracy == 0:
            z21 = 0
            z22 = 0
        else:
            is_error = 1
                      
    if j > (x_dim - 1) or i > (y_dim - 1):
        is_error = 1

    # Test code
    if "z12" not in locals():
        arcpy.AddMessage("i:%s, j:%s" % (i,j))
        arcpy.AddMessage("fracx = %s" % fracx)
        #sys.stdout.flush()

    if is_error == 1:
        # This error should not ever be reached - there is an earlier check that should prevent these locations from being entered, with a more detailed error messsage
        raise ValueError("A source or receiver location is outside the supported landscape (must be >= 1/2 cell width/length from the edge of the landscape)")

    return(z11,z21,z12,z22)

# Get a transect for calculating NMSim results
def get_transect(xyzsrc, xyzrec, elevation, impedance, dem_name, point_dir, cell_x, cell_y, coord_system, do_timing, my_times, my_time_labels, count):
    import math

    # calculate max_dist as distance between xyzsrc & xyzrec (without respect to elevation)
    xsource = xyzsrc[0]
    ysource = xyzsrc[1]
    xreceiver = xyzrec[0]
    yreceiver = xyzrec[1]
    max_dist = math.sqrt( (xreceiver - xsource)**2 + (yreceiver - ysource)**2 ) #**# can add z if needed

    #Calculate minimum cell width
    #cell_width = min([cell_x, cell_y])

    # Set up number of points along transect # Changed on 2016-04-21 based on input from Bruce Ikelheimer about the ONECUT algorithm
    #points_per_cell_width = 1 # Provide a variable that is easy to change that controls the resolution of points along the transect    
    npts = 100  # Bruce Ikelheimer said that ONECUT algorithm is ALWAYS given 100 points.
    #npts = max_dist * cell_width**-1 * points_per_cell_width + 1 # Plus 1 adds a point at the origin
    #npts = int(math.ceil(round(npts,3))) # There will be variable spacing between the last point and the receiver. But we need to include the origin & receiver. But round is there to ensure that slight rounding differences does not cause the code to add extra points
    increment = max_dist * (npts - 1)**-1 # - 1 accounts for the inclusion of the end points, **-1 does division
    #print "max_dist = %s, increment = %s" % (max_dist, increment)

    if count == 0:
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "get_transect initial calculations completed")
     
    # create dist_vec #**# Watch for problems with rounding to 4 decimal places. This was to prevent really ugly repititions in numbers
    dist_vec = seq(0,max_dist, increment) #cell_width * points_per_cell_width**-1)
    for i in range(len(dist_vec)):
        dist_vec[i] = round(dist_vec[i], 8)
    
    # check that the final point was added to dist_vec (if dist_vec is one short, add the final point)
    # Sometimes there is a mismatch due to rounding, and the last point is not added.
    if len(dist_vec) == (npts - 1):
       dist_vec.append(max_dist)

    if count == 0:
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "get_transect dist_vec created")
    
    # Create points that correspond to dist_vec for purpose of extracting heights from the dem
    dist_points = make_transect_points(xyzsrc, xyzrec, npts, point_dir, coord_system, my_times, my_time_labels, count)    

    if count == 0:
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "get_transect dist_points created")

    # Sample dem at points to get hgt_vec #**# try to replace the arcgis steps with numpy code???
    hgt_vec, my_times, my_time_labels = sample_extract(elevation, dist_points, point_dir, my_times, my_time_labels, count)

    if count == 0:
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "get_transect hgt_vec created")

    # CHANGED THIS - have an earlier step convert landcover to impedance, then just sample the impedance.
    # Sample landcover at points to get landcover    
    imp_vec, my_times, my_time_labels = sample_extract(impedance, dist_points, point_dir, my_times, my_time_labels, count)    
    #lc_vec = sample_extract(landcover, dist_points, point_dir)

    if count == 0:
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "get_transect imp_vec created")
    
    # Get imp_vec from landcover at sample points
    #imp_vec = landcover_to_impedance(lc_vec)

    # Check that everything is the correct length
    if len(dist_vec) != npts or len(hgt_vec) != npts or len(imp_vec) != npts:
        raise ValueError("npts (%s), dist_vec (%s), hgt_vec (%s), and imp_vec (%s) all must be the same length." % (npts, dist_vec, hgt_vec, imp_vec))

    # Sampling appears to occur from top to bottom, and then in order of point creation
    # Thus, when top-to-bottom order reverses the vectors, we need to flip them back
    if yreceiver > ysource:
        hgt_vec = hgt_vec[::-1] # you can also write hgt_vec = list(reversed(hgt_vec))
        imp_vec = imp_vec[::-1]

    return(npts, max_dist, dist_vec, hgt_vec, imp_vec, my_times, my_time_labels)

# Get points along a transect
def make_transect_points(xyzsrc, xyzrec, npts, point_dir, coord_system, my_times, my_time_labels, tcount):

    import arcpy

    if tcount == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "Start make transect points")


    # get origin x & y
    x_source = xyzsrc[0]
    y_source = xyzsrc[1]
    x_receiver = xyzrec[0]
    y_receiver = xyzrec[1]

    # get difference between source & receiver
    xdiff = x_receiver - x_source
    ydiff = y_receiver - y_source    

    # if statement prevents dividing by zero (the loop that uses x_inc will not be entered with only one point)
    if round(xdiff, 3) == 0 and round(ydiff, 3) == 0:
        x_inc = 0
        y_inc = 0
    else:
        # get x & y increments (should include sign!)
        x_inc = xdiff * (npts - 1)**-1
        y_inc = ydiff * (npts - 1)**-1

    ''' No longer needed?
    # Need to know if adding cell size or subtracting it from x & y
    if xdiff >= 0:
        xdir = 1
    else:
        xdir = -1
        
    if ydiff >= 0:
        ydir = 1
    else:
        ydir = -1
    '''

    x_value = x_source
    y_value = y_source

    if tcount == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "About to write to file")


    #**# think about converting from numpy array (or figuring out how to write a shapefile in python)\
    points_table = point_dir + "point_table%s.csv" % tcount
    with open(points_table, 'w') as pt:
        count = 0
        sep = ","
        pt.write("X%sY\n" % sep) # Write header

        # -1 is because npts includes the final point, which is added manually below
        for i in range((npts - 1)):
            pt.write("%s%s%s\n" % (x_value, sep, y_value)) 
            x_value = x_value + x_inc # * xdir
            y_value = y_value + y_inc # * ydir   
            count += 1
            
        # Add final point
        pt.write("%s%s%s\n" % (x_receiver, sep, y_receiver))


    if tcount == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "File written")

    # Convert table to shapefile
    transect_lyr = "transect_points%s" % tcount
    transect_points = point_dir + "transect_points%s.shp" % tcount

    #**# Check if this is rate-limiting
    #**# Try replacing with pyshp package's tools, I bet it's MUCH faster.
    #**# No spatial reference for this point layer
    arcpy.MakeXYEventLayer_management (points_table, "X", "Y", transect_lyr) 
    arcpy.CopyFeatures_management(transect_lyr, transect_points)
    arcpy.DefineProjection_management(transect_points, coord_system)
    #arcpy.AddField_management(transect_points, "ID", "LONG")
    #arcpy.CalculateField_management(transect_points, "ID","FID", "PYTHON_9.3")

    if tcount == 0:
        #**# Do timing is fixed to 1, you'll need to fix that to not run timing when you move out of testing
        my_times, my_time_labels = sg_get_time(1, my_times, my_time_labels, "File converted to .shp")

    return(transect_points)


# Call NMSim propagation algorithm
def NMSim_propagation(nmsim_dll, zsrc, zrec, npts, dist_vec, hgt_vec, imp_vec,
                      freq, nfreq):
    nmsim = cdll.LoadLibrary(nmsim_dll) 
    
    #Load the data into the ctype variables
    npts_c = ctypes.c_int8(npts)  #Number of distances
    dist_c = (ctypes.c_float * npts)(*dist_vec)
    imp_c = (ctypes.c_float * npts)(*imp_vec)
    hgt_c = (ctypes.c_float * npts)(*hgt_vec)
    freq_c = (ctypes.c_float * nfreq)(*freq)
    nfreq_c = ctypes.c_int8(nfreq)
    zsrc_c = ctypes.c_float(zsrc)
    zrec_c = ctypes.c_float(zrec)
    
    dummy = [0] * nfreq
    #Output
    attbnd = (ctypes.c_float * nfreq)(*dummy)
    
    #Call to NMSim routine ONECUT
    nmsim.ONECUT(attbnd, ctypes.byref(npts_c), dist_c, hgt_c, imp_c, ctypes.byref(zsrc_c), ctypes.byref(zrec_c), ctypes.byref(nfreq_c),
                 freq_c)
    
    #Convert attbnd to python list
    # Extract data from modified absatm object
    attbnd_list = dummy
    for i in range(0,nfreq):
        attbnd_list[i]= -attbnd[i] 
    #attbnd_list = remove_negatives(attbnd_list)
    
    return attbnd_list                 
   

# Call NMSim atmospheric absorption algorithm
def NMSim_atmospheric_absorption(nmsim_dll, temp, rh, freq, nfreq, max_dist):
    nmsim = cdll.LoadLibrary(nmsim_dll) 
    dummy = [0] * nfreq
    
    # Copy the data over into ctype variables
    temp_c = ctypes.c_float(temp)
    rh_c = ctypes.c_float(rh)
    nfreq_c = ctypes.c_int8(nfreq)
    freq_c = (ctypes.c_float* nfreq) (*freq)  #This is the format necessary to correctly copy the list into the ctype array
    absatm = (ctypes.c_float* nfreq) (*dummy)

    ctypes.c_short

    # When calling the dll, floats and integers must be passed by reference. Everything else is a simple pass.
    # this changes absatm by list reference
    nmsim.SETATM(absatm, freq_c, ctypes.byref(nfreq_c), ctypes.byref(temp_c), ctypes.byref(rh_c))
    
    # Extract data (dB/m) from modified absatm object
    absatm_list = [0] * nfreq
    for i in range(0,nfreq):
        absatm_list[i]=absatm[i]    
    
    # Multiply by distance to receiver (m) to get total absorption (dB)
    absatm_out = [x * max_dist for x in absatm_list]
    #absatm_out = remove_negatives(absatm_out)
    
    return(absatm_out)


# Call NMSim source frequency generation algorithm
def NMSim_Source(nmsim_dll, xyzsrc, xyzrec, zsrc, zrec, head, roll, vel, pitch,
                 engpow, srcfile, nfreq):
    
    nmsim = cdll.LoadLibrary(nmsim_dll)

    #arcpy.AddMessage("%s" % srcfile)
    #srcfile_name = os.path.basename(srcfile)

    # Copy the data over into ctype variables
    #input variables
    xyzsrc_c = (ctypes.c_float * 3)(*xyzsrc)
    xyzrec_c = (ctypes.c_float * 3)(*xyzrec)
    zsrc_c = ctypes.c_float(zsrc)
    zrec_c = ctypes.c_float(zrec)
    head_c = ctypes.c_float(head)
    roll_c = ctypes.c_float(roll)
    vel_c = ctypes.c_float(vel)
    pitch_c = ctypes.c_float(pitch)
    engpow_c = ctypes.c_float(engpow)
    srcfile_c = ctypes.c_char_p(srcfile)
    nfreq_c = ctypes.c_int8(nfreq)
    str_len_c = ctypes.c_int8(len(srcfile)) #Fortran wants to know how long the string is
    dummy = [0] * nfreq

    #output variables
    refdist = ctypes.c_float(0) #The reference distance of this source
    distance = ctypes.c_float(0) #The distance between the source and the receiver
    dbband = (ctypes.c_float * nfreq)(*dummy) #The source levels at the given frequencies
    
    # set working directory to help out NMSIMSOURCE    
    main_working_directory = os.getcwd()
    temp_working_directory = os.path.dirname(os.path.realpath(srcfile)) # get path of .src file (os.path.realpath gets the path including the file name, os.path.dirname gets just the directory without a trailing slash.)
    os.chdir(temp_working_directory)    

    #Call the source data
    nmsim.NMSIMSOURCE(dbband, ctypes.byref(refdist), ctypes.byref(distance), xyzsrc_c, xyzrec_c, ctypes.byref(zsrc_c),
                      ctypes.byref(zrec_c), ctypes.byref(head_c), ctypes.byref(roll_c), ctypes.byref(vel_c), ctypes.byref(pitch_c),
                      ctypes.byref(engpow_c), srcfile_c, ctypes.byref(str_len_c),  #You need to pass the string, AND its length
                      ctypes.byref(nfreq_c))
    # Reset working directory to what it had been.    
    os.chdir(main_working_directory)    
    
    # Extract refdist, distance & dbband
    refdist_out = refdist.value
    distance_out = distance.value
    dbband_out = [0] * nfreq
    for i in range(0,nfreq):
        dbband_out[i] = dbband[i]    
    #dbband_out = remove_negatives(dbband_out)
    
    return(refdist_out, distance_out, dbband_out)    


# Define a function to remove negative values from a list
def remove_negatives(in_list):
    new_list = list(in_list)
    for i in range(len(new_list)):
        if new_list[i] < 0:
            new_list[i] = 0
    
    return(new_list)
    

def sg_get_time(do_timing, my_times, my_time_labels, label):
    import time    
    if do_timing == 1:
        this_time = time.time()
        my_times.append(this_time)
        my_time_labels.append(label)

    return(my_times, my_time_labels)


# Added core function for running the ISO 9613-2. Modified from nmsimgis_core function
def iso_core(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, zsrc, zrec, head, roll, pitch, vel, engpow,
                  srcfile, temp, rh, npts, max_dist, dist_vec, hgt_vec,
                  imp_vec, freq, nfreq, nmsim_dll, point_dir, ambient_vals = ["none"]):


    temp = float(temp) # needs to be in float format
    rh = float(rh) # Needs to be in float format
    freq = map(float, freq) # Convert list to float format

    #Reconvert xyzsrc and xyzrec from cell units to m
    xsrc_m = xyzsrc_cell[0] * cell_x
    ysrc_m = xyzsrc_cell[1] * cell_y
    xyzsrc = [xsrc_m, ysrc_m, xyzsrc_cell[2]] # z coordinate was ok
    
    xrec_m = xyzrec_cell[0] * cell_x
    yrec_m = xyzrec_cell[1] * cell_y
    xyzrec = [xrec_m, yrec_m, xyzrec_cell[2]] # z coordinate was ok    

    # Calculate 3D distance between source and receiver
    dist_3d = math.sqrt((xsrc_m - xrec_m)**2 + (ysrc_m - yrec_m)**2 + (xyzsrc[2] + zsrc - xyzrec[2] - zrec)**2)
       
    # Set up source information for later propagation analysis
    refdist_out, distance_out, dbband_out = NMSim_Source(nmsim_dll, xyzsrc, xyzrec, zsrc,
                                           zrec, head, roll, vel, pitch, engpow,
                                           srcfile, nfreq)
   
    # Evaluate atmospheric absorption for a single receiver point
    absorption = NMSim_atmospheric_absorption(nmsim_dll, temp, rh, freq,
                                                            nfreq, dist_3d)
    
    # Evaluate sound propagation based on terrain, impedance, source height,
    # and receiver height
    # Includes ground effect and barrier effect together
    attenuation = ISO_propagation(xyzsrc, xyzrec, zsrc, zrec, dist_vec, hgt_vec, imp_vec, freq, temp)
    
    # Reflections
    # Skipped - this will be complicated to code and I'm not sure how widely applicable it will be. It would be easier to code for special cases
    # than generally for every point in the landscape

    # Meteorological correction needs to be applied later - requires an A-weighted value
    
    # Meteorological correction
    #**# This is for the A-weighted sound level, and not for individual frequency bands - needs to be moved!
    # Need to add two more arguments - whether or not to apply the correction and the C0 coefficient.
    # They could just be one argument - NA if not to apply it, and C0 if you are to apply it.    
    #Cmet = ISO_meteorological_correction(C0, xyzsrc, xyzrec, zsrc, zrec)
    
    
    # Add spherical spreading loss (this is even across frequencies)    
    if distance_out == 0:
        distance_out = 1 # Switch to 1 m distance if the source point is exactly on the midpoint of the cell of origin
    spherical_spreading_loss = 20 * math.log10(distance_out * refdist_out**-1) #do = source reference distance (e.g., spl measured at 15m), d = distance between source & receiver
    
    
    # Combine frequency loss to get final decibel values (convert to function)
    # I'm sure there are better ways to code this
    db_lst = [0] * nfreq
    for i in range(nfreq):
        this_source_level = dbband_out[i]
        this_atmoabs = absorption[i]
        this_attenuation = attenuation[i]
        #this_ground = Aground[i] 
        #this_bar = Abarrier[i] 
    
        #db_out = this_source_level - this_atmoabs - this_ground - this_bar - spherical_spreading_loss
        db_out = this_source_level - this_atmoabs - this_attenuation - spherical_spreading_loss

        # Cap out extreme negative values to avoid code crashes
        if db_out < -300:
            db_out = -300
        
        db_lst[i] = db_out

    # convert negative values to 0's
    #db_lst = remove_negatives(db_lst) # Not really appropriate to do this, especially before summing up for overall results.

    # Output results of sound propagation for use by end users
    return(db_lst, dbband_out, attenuation, absorption, spherical_spreading_loss)

# Call ISO 9613-2 terrain/ground propagation algorithm
def ISO_propagation(xyzsrc, xyzrec, zsrc, zrec, dist_vec, hgt_vec, imp_vec, freq, temp_C):

    # Calculate ground effects
    Aground = ISO_ground_effects(xyzsrc, xyzrec, zsrc, zrec, dist_vec, imp_vec, freq)
    
    # Calculate barrier effects
    Abarrier = ISO_barrier_effects(Aground, freq, hgt_vec, dist_vec, xyzsrc, xyzrec, zsrc, zrec, temp_C)
    
    # Calculate attenuation due to ground and barrier effects
    attenuation = [0] * len(freq)
    for f in xrange(len(freq)):
        attenuation[f] = Aground[f] + Abarrier[f]
    
    return(attenuation)
    

#**# MOVE TO SOUNDPROPHLPR - NEEDS TO BE APPLIED AFTER SUMMARY STEPS
# ISO meteorological correction
def ISO_meteorological_correction(C0, xyzsrc, xyzrec, zsrc, zrec):
    if C0 == "NA":
        Cmet = 0
    else:
        # Get distance ("the distance between the source and receiver projected to the horizontal ground plane, in metres")
        #**# I don't think this is the 2D distance. This is the ground distance - so I think you use the ground elevation for this, as the ISO assumes flat (but not necessarily horizontal) ground.        
        #dp = math.sqrt( (xyzsrc[0] - xyzrec[0])**2 + (xyzsrc[1] - xyzrec[1])**2 )        
        dp = math.sqrt( (xyzsrc[0] - xyzrec[0])**2 + (xyzsrc[1] - xyzrec[1])**2 + (xyzsrc[2] - xyzrec[2])**2 )
        hs = zsrc # Source height above ground level
        hr = zrec # Receiver heigh above ground level
        test = 10 * (hs + hr)
        if dp <= test:
            Cmet = 0
        
        if dp > test:
            Cmet = C0 * ( 1 - test / dp)
    
    return Cmet

# Calculate ground effects according to ISO 9613-2
def ISO_ground_effects(xyzsrc, xyzrec, zsrc, zrec, dist_vec, imp_vec, freq):

    # Calculate straight line distance between source and receiver, ignoring source and receiver heights
    # I'm assuming this is what "the distance between the source and receiver projected to the horizontal ground plane, in metres" means
    dp = math.sqrt( (xyzsrc[0] - xyzrec[0])**2 + (xyzsrc[1] - xyzrec[1])**2 + (xyzsrc[2] - xyzrec[2])**2 )
    
    # Calculate distances that are within the source and receiver regions
    source_region = 30 * zsrc
    receiver_region = 30 * zrec
    
    # Calculate which distances from dist_vec correspond to the source region
    # Go through distances, until you find one outside the source region
    for i in xrange(len(dist_vec)):
        this_dist = dist_vec[i]
        
        if this_dist > source_region:
            break
    
    # Now that the index above has been identified, use it to extract impedances from imp_vec
    source_imp_vec = imp_vec[0:i] # in Python, the last index is excluded. This works to our advantage, as the last i value corresponds to a value outside the source region
    Gsource = calculate_G(source_imp_vec) 
    Asource = calculate_ground_effect(freq, Gsource, zsrc, dp)
    
    # Calculate the distances from dist_vec that correspond to the receiver region
    last_index = len(dist_vec) - 1 # -1 corrects for the last element not being included
    last_dist = dist_vec[-1]
    for j in xrange(last_index, -1, -1): # Go from the last index to 0. -1 is not included, as it is the last number in the index. The second -1 indicates the increment to decrease by
        this_dist = dist_vec[j]
        # Get distance from receiver, not distance from source
        corrected_dist = last_dist - this_dist

        if corrected_dist > receiver_region:
            break
    start_element = j + 1 # Need to add one, because j is 1 lower than the start of the range
    receiver_imp_vec = imp_vec[start_element:] # in Python, the last index is excluded. This works to our advantage, as the last i value corresponds to a value outside the source region
    Greceiver = calculate_G(receiver_imp_vec)  
    Areceiver = calculate_ground_effect(freq, Greceiver, zrec, dp)

    # identify middle ground, if any
    if i >= (j + 1):
        Amiddle = [0] * len(freq)
    else:
       middle_imp_vec = imp_vec[i:j+1]
       Gmiddle = calculate_G(middle_imp_vec)
       Amiddle = calculate_middle_ground_effect(freq, Gmiddle, zsrc, zrec, dp)

    Aground = ["NA"] * len(freq)
    for i in xrange(len(Aground)):
        Aground[i] = Asource[i] + Amiddle[i] + Areceiver[i]
        
    return Aground

# Calculate G coefficient
def calculate_G(imp_vec):
    # Hard ground G = 0
    # Porous ground G = 1
    # Mixed ground G = fraction porous

    # Loop through impedance vector and identify ground as hard or soft
    imp_sum = 0
    for imp in imp_vec:
        # >= 10,000 hard ground as this corresponds to high intensity developed lands according to Bruce Ikelheimer
        # <= 200 soft ground - corresponds to natural vegetation layers
        # < 10,000 > 200 corresponds to mixed ground
        #    # 300 is ~70% soft 30% hard, # 5000 is ~35 soft 65 hard
        # But to keep things simple, and because I don't have the data to model
        # the non-linearity described above, I'm just going to make it all act as 50% hard/soft
    
        if imp >= 10000:
            this_G = 0
        elif imp <= 200:
            this_G = 1
        else:
            this_G = 0.5
        imp_sum += this_G
        
    G = imp_sum * len(imp_vec)**-1 # divide the sum by the number of elements to get the fraction corresponding to ground hardness 
    
    return G 

# Calculation for getting a ground effect from the ground effect coefficient
def calculate_ground_effect(freq, G, h, dp):
    '''
    freq: The frequencies for consideration. This will be input as 1/3 octave, but the ISO gives calcualtions for octave bands
    G: Ground effect coefficient, which ranges from 0 (hard) to 1 (soft)
    location: source, middle, receiver. Source & receiver are equivalent, and really only middle vs. not middle matters
    height: required for source and receiver, not used for middle calculation
    '''
    # Create list to hold outputs
    Acomponent = ["NA"] * len(freq)
    # Loop through frequencies
    for fr in xrange(len(freq)):
        f = freq[fr]
        
        if f < 50:
            arcpy.AddError("ISO 9613-2 does not support frequency bands below 50 Hz")
        
        # Map 1/3 octave bands to their relevant octave band
        if f == 50 or f == 63 or f == 80:
            A = -1.5
        
        if f == 100 or f == 125 or f == 160:
            A = -1.5 + G * aprime(h, dp)
        
        if f == 200 or f == 250 or f == 315:
            A = -1.5 + G * bprime(h, dp)
            
        if f == 400 or f == 500 or f == 630:
            A = -1.5 + G * cprime(h, dp)
        
        if f == 800 or f == 1000 or f == 1250:
            A = -1.5 + G * dprime(h, dp)
        
        if f >= 1600:
            A = 1.5 * (1 - G)

        # Assign this A value to the Acompenent slot for this frequency
        Acomponent[fr] = A

    return Acomponent

# Calcualte aprime function
def aprime(h, dp):
    term1 = math.exp(-0.12 * (h - 5)**2)
    term2 = 1 - math.exp(-dp / float(50)) #float is necessary to get normal division
    term3 = math.exp(-0.09 * h**2)
    term4 = 1 - math.exp(-2.8 * 10**-6 * dp**2)
    
    value = 1.5 + 3.0 * term1 * term2 + 5.7 * term3 * term4
    return value

# Calculate bprime function
def bprime(h, dp):

    term1 = math.exp(-0.09 * h**2)
    term2 = 1 - math.exp(-dp / float(50)) #float is necessary to get normal division
    
    value = 1.5 + 8.6 * term1 * term2
    return value

# Calculate cprime function
def cprime(h, dp):
    
    term1 = math.exp(-0.46 * h**2)
    term2 = 1 - math.exp(-dp / float(50)) #float is necessary to get normal division
    
    value = 1.5 + 14.0 * term1 * term2
    return value

# Calculate dprime function
def dprime(h, dp):
    term1 = math.exp(-0.9 * h**2)
    term2 = 1 - math.exp(-dp / float(50)) #float is necessary to get normal division
    
    value = 1.5 + 5.0 * term1 * term2
    return value

# Calculate ground effect for the middle
def calculate_middle_ground_effect(freq, Gmiddle, hsource, hreceiver, dp):
    Amiddle = [0] * len(freq) # create a list of values
    
    # define q
    # Ignore the criterion for q = 0 - if we're at this point, that has already been satisfied
    top = 30 * (hsource + hreceiver)
    bottom = dp
    q = 1 - top / bottom    
    
    # Loop through frequencies
    for fr in xrange(len(freq)):
        f = freq[fr]
        
        # Different behavior for the 63 Hz octave band than for the rest
        if f <= 80:
            # NOTE the 2) footnote in the ISO document looks a lot like raising
            # q to the 21 power. But this is a footnote, not part of the equation.
            A = -3 * q
        else:
            A = -3 * q * (1 - Gmiddle)
        
        Amiddle[fr] = A
    
    return Amiddle

# Calculate barrier effects according to ISO 9613-2    
def ISO_barrier_effects(Aground, freq, hgt_vec, dist_vec, xyzsrc, xyzrec, zsrc, zrec, temp_C):
    # Assumes:
    # 1) Surface density is at least 10 kg/m2
    # 2) Object has closed surface without cracks or gaps
    # 3) Horizontal dimension of object is larger than acoustic wavelength lambda

    # 1 & 2 are met by most terrain barriers (earth is heavy, and with exceptions of arches and holes, not cracked or gapped typically)
    # 3 would be computationally intensive to test, but could be untrue.
    # 63 Hz has a wavelength of 5.46 m. This is smaller than the typical cell
    # resolution, so the assumption that the barrier is greater than the wavelength
    # is reasonable for any barriers encountered given the resolution of the data.
    
    # We omitted diffaction around vertical edges - this would be computationally
    # difficult to calculate for every cell on the landscape
        
    # Create a list to contain barrier loss output
    Abarrier = [0] * len(freq)    

    # identify if there is a barrier, and if so, where it is
    barrier_height, barrier_distance, slope = identify_barrier(dist_vec, hgt_vec, xyzsrc, xyzrec, zsrc, zrec)
    
    # If barrier height is 0, leave Abarrier at 0. Otherwise, calculate barrier loss
    if barrier_height != 0:
    
        # Get distance to the barrier from source & receiver
        dss, dsr = get_barrier_distances(barrier_height, barrier_distance, slope, dist_vec[-1], xyzsrc, xyzrec, zsrc, zrec)    
    
        # Loop  through frequencies
        for fr in xrange(len(freq)):
            f = freq[fr]
            c = calculate_speed_of_sound(temp_C) # speed of sound (m/s), calculation imported from soundprophlpr
            wavelength = c / f # c is a float, so this should divide correctly
        
            # Calculate z and C3 constant (depend on single vs. double refraction)
            z, C3 = calculate_z_C3(dss, dsr, dist_vec, wavelength)
                
            Kmet = calculate_Kmet(dss, dsr, z, dist_vec[-1])
                        
            Dz = 10 * math.log(3 + (20/wavelength) * C3 * z * Kmet)
            
            # ISO 9613-2 caps single refraction at 20 dB (double-refraction barriers at 25 dB)
            if Dz > 20:
                Dz = 20
            
            Abar = Dz - Aground[fr]
                            
            Abarrier[fr] = Abar
    
    return Abarrier

# Calculate z used in calculation of barrier effect by ISO 9613-2
def calculate_z_C3(dss, dsr, dist_vec, wavelength):
    
    # For now, assume single diffraction - computationally easier!    
    # identify whether to use single or double diffraction
    #diffraction_type, dss, dsr, e, a = compute_diffraction(hgt_vec, dist_vec)
    diffraction_type = "single"
    d = dist_vec[-1] # Get distance between source & receiver

    # The ISO contains a variable a.             
    # "a is the component distance parallel to the barrier edge between source and receiver, in metres"
    # In this implementation, a is already incorporated in the calculation of dss and dsr, so here it is set to 0
    a = 0 
    
    if diffraction_type == "single":
        C3 = 1
        z = math.sqrt((dss + dsr)**2 + a**2) - d
    
    # Not used - assuming single diffraction    
    #if diffraction_type == "double":
    #     top = 1 + ((5 * wavelength) / e)**2
    #     bottom = 1*3**-1 + (5 * wavelength / e)**2
    #     
    #     C3 = top * bottom**-1
    #     
    #     z = math.sqrt((dss + dsr + e)**2 + a**2) - d
         
    return(z, C3)

# Calculate meteorological coefficient
def calculate_Kmet(dss, dsr, z, distance):
    
    if z <= 0:
        Kmet = 1
    else:
        term1 = -(1*2000**-1)        
        term2 = math.sqrt((dss * dsr * distance) * (2 * z)**-1)
        Kmet = math.exp(term1 * term2)
    
    return Kmet

# Compute whether single or double refraction is more appropriate
def compute_diffraction(hgt_vec, dist_vec):
    raise ValueError("We decided to assume single diffraction for simplicity")

# Calculate dss and dsr   
# Code copied from spreadgishlpr.py 
def identify_barrier(dist_vec, hgt_vec, xyzsrc, xyzrec, zsrc, zrec):

    # Get elevation of the source & receiver    
    source_elevation = xyzsrc[2] + zsrc
    receiver_elevation = xyzrec[2] + zrec

    # Calculate the slope between the source and that cell
    distance = dist_vec[-1]
    # If the source and target cell are the same, there is no barrier
    if distance <= 0:
        max_height = 0
        bar_dist = 0        
        slope = 0
    # Otherwise, calculate the slope, barrier height, and barrier distance
    else:
        slope = (receiver_elevation - source_elevation) * distance**-1
        
        max_height = 0 # Negative heights do not obstruct sound propagation
        bar_dist = distance
        # Start at 1, source point cannot be a barrier
        for i in xrange(1, len(dist_vec)):
            slope_elevation = slope * dist_vec[i] + source_elevation # Slope * distance gives the increase or decrease in height between source and receiver
            height_above_slope = hgt_vec[i] - slope_elevation
            
            # Calculate the maximum height above the slope along the terrain cut
            if height_above_slope > max_height:
                max_height = height_above_slope
                bar_dist = dist_vec[i]
    
    return(max_height, bar_dist, slope)


# Compute distance between source, receiver and barrier
def get_barrier_distances(barrier_height, barrier_distance, slope, distance, xyzsrc, xyzrec, zsrc, zrec):

    # Get proportion along the path between source and receiver corresponding to the barrier location
    bar_prop = barrier_distance * distance**-1

    # Calculate x, y, z distances to barrier
    dss_x = (xyzsrc[0] - xyzrec[0]) * bar_prop # Get x distance from source
    dsr_x = (xyzsrc[0] - xyzrec[0]) * (1 - bar_prop) # Get x distance from receiver

    dss_y = (xyzsrc[1] - xyzrec[1]) * bar_prop  # Get y distance from the source
    dsr_y = (xyzsrc[1] - xyzrec[1]) * (1 - bar_prop) # Get y distance from receiver

    source_elevation = xyzsrc[2] + zsrc
    receiver_elevation = xyzrec[2] + zrec
    barrier_elevation = slope * barrier_distance + source_elevation # Slope * distance gives the increase or decrease in height between source and receiver

    dss_z = source_elevation - barrier_elevation
    dsr_z = receiver_elevation - barrier_elevation

    # Get distance from source to top of the barrier
    dss = math.sqrt( dss_x**2 + dss_y**2 + dss_z**2 ) 
    
    # Get distance from receiver to the top of the barrier
    dsr = math.sqrt( dsr_x**2 + dsr_y**2 + dsr_z**2 ) 
    
    return(dss, dsr)
    
    
# End of file