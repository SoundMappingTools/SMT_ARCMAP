# -*- coding: utf-8 -*-
"""
Description: Functions to support running SPreAD-GIS, NMSIMGIS, ISO 9613-2 and
the tools included in the Extra Tools section of the toolbox.

Dependencies: ArcGIS >=10.3, Python 2.7
Dependencies for some applications: spreadgishlpr.py, nmsimhlpr.py,
                                    nmsimgis_analysis.py, NumPy

Development History:
Sarah Reed developed the SPreAD-GIS model. Some of that code is retained and
modified in functions here, while other pieces of that code have been moved
to spreadgishlpr.py. It has been expanded by A.C. Keyel to provide support
to SPreAD-GIS, NMSIMGIS, and ISO 9613-2 with a minimum of code repitition.
Tools used by the Extra Tools are also included here. These tools are designed
to be accessed directly through Python or indirectly through interface scripts
that support the ArcGIS toolbox GUI. Code for calculating weights was modified
from code written by Bruce Ikelheimer, provided under the GPL license.

@author Alexander "Sasha" Keyel <skeyel@gmail.com>, Sarah E. Reed,
        Bruce Ikelheimer
Copyright (C) 2010 - 2016 Sarah E. Reed, A.C. Keyel <skeyel@gmail.com>, and
                          Bruce Ikelheimer
                          
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

# Import required modules once at start of file
import sys, os, shutil, math, subprocess
import arcpy, arcpy.sa
# Check out the spatial analyst extension
arcpy.CheckOutExtension("Spatial")
# Delete intermediate datasets, if they exist
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False  # Prevent intermediates from being added to the map mid-processing

# Main Sound Mapping Tools call, this can call either SPreAD-GIS or nmsimgis
# NOTE: inputs are now in metric, but many of the actual calculations are still based on ft
def SoundPropagation(model, point_source_file, source_id_field, Model_Extent,
                     ambient_path, elevation, landcover, temp, rh, wind_dir,
                     wind_sp, seas_cond, base_dir, results_dir, results_label,
                     in_freq, receiver_offset = 1, n_points = "all",
                     tbx_root = "C:/smt/toolbox/", timelog = "none",
                     delete_existing_timelog = 1, summarize = 0,
                     keep_intermediates = 0, weighting = ["A"],
                     n_processes = 1, propagation_speed_output = 1,
                     truncate = 0, landcover_type = "nlcd", is_GUI = "NO",
                     use_old_barrier = 1, allow_defaults = "YES"):
    # ambient_path = "NA" bypasses comparisons to the ambient environment
    
    '''
    # for troubleshooting
    from soundprophlpr import *
    from nmsimhlpr import *    

    # For ???
    in_freq = srcfile    
    point_source_file = sources
    Model_Extent = lc_extent
    base_dir = opath
    summarize = 0
    keep_intermediates = 1

    # For sensitivity analysis
    temp = max_weather[1]
    rh = max_weather[2]
    
    # For scaling analysis
    point_source_file = source1 
    source_id_field = source1_id
    source_info= drill_source_info
    source_offset = drill_AGL_ft
    in_freq = drill_src
    summarize = 2
    '''

    # Check that the base_dir does not exceed the ESRI limit. If it does, throw an error and request the user to use a shorter file path name
    #**# In the long run, convert the intermediates to .tif to increase or eliminate this limitation!
    #**# May be overly restrictive - some intermediates have been converted to .tif
    check_path_length(base_dir)

    # Check if inputs are in a consistent projected coordinate system
    #**# This does not check ambient inputs
    check_projections([point_source_file, elevation, landcover])

    # Set general parameters for rasters
    dscRD = arcpy.Describe(elevation)
    #RasterExtent = dscRD.Extent
    CellSize = dscRD.MeanCellHeight

    # Set Environment settings # Now set here for all tools instead of in each tool.
    arcpy.env.extent = Model_Extent #RasterExtent
    arcpy.env.cellSize = CellSize
    arcpy.env.snapRaster = arcpy.Raster(elevation)

    # Set up timing for analysis. If timelog == "none" no timing will be performed
    my_times, my_time_labels = setup_timing(timelog)

    # Set up paths
    model_dir, idir, pdir, freq_dir = dir_setup(base_dir, results_dir, model, summarize)

    # Clip dem to desired extent "NO_MAINTAIN_EXTENT" prevents resampling the rasters & may result in a slightly larger final extent
    dem_clip = model_dir + "dem_clip"
    if summarize < 3:
        arcpy.Clip_management(elevation, Model_Extent, dem_clip, "", "", "", "NO_MAINTAIN_EXTENT")

        # Check if model is anticipated to take a "long" time and give a warning
        estimate_model_run_time(dem_clip, CellSize, model)

        # Create dem_ft needed by spreadgishlpr
        if model == "spreadgis":
            dem_ft = idir + "dem_ft"
            this_dem = arcpy.Raster(dem_clip) * 3.28084 # Conversion values taken from Google on 2016-03-16
            this_dem.save(dem_ft)   
            
    # Change extent to match the clipped extent - see if this helps?? (error in nmsimgis where extents don't match, even though there is NO reason for them not to!)
    Model_Extent = get_extent(dem_clip)
    arcpy.env.extent = Model_Extent #RasterExtent

    # Set up points information
    # set up inputs for multipoints
    point_lst = get_points(point_source_file, source_id_field, n_points)
    Sound_Source = idir + "source_point.shp"

    # extract list of frequencies # note, minimum receiver height is also extracted from the source file when reading from a .src file
    # Added summarize option, as otherwise freq_lst is undefined, and for summarization, we just need a list of the frequencies.
    if model == "spreadgis" and in_freq[-4:] != ".src" or summarize >= 3:    
        #**# Value weighting is deprecated, not recommended for use and may be removed in future versions.
        freq_lst, sound_level_lst, measurement_distance_lst, min_source_height, value_weighting = extract_frequencies(in_freq)

    # Add constants for adding leading 0's, but not more than 3 to avoid problems with file names being too long: 13 char max
    #**# This needs future attention, as model runs with > 4 digits could crash after a long computation time.
    point_fill = len(str(len(point_lst))) # Get number of digits. len of str is number of characters. len of pointlist gets number of elements. So it's the number of digits in len(point_lst)     
    point_fill = min(point_fill, 3)
    freq_fill = 5

    if summarize < 3:

        # Reclassify landcover for the user, if desired
        if landcover_type != "default":
            reclass_landcover(landcover, landcover_type, model)

        #if model == "nmsimgis":
            
        # Clip landcover to match elevation raster (#Snap raster should ensure everything lines up)
        land_clip_pre = idir + "land_clip_pre"
        land_clip = idir + "land_clip"
        arcpy.Clip_management(landcover, Model_Extent, land_clip_pre, "", "", "", "MAINTAIN_EXTENT") # Formerly land_clip_pre, no_maintain_extent
        #The second clip step attempts to remove an extra row of added cells that should not have been there in the first place
        arcpy.Clip_management(land_clip_pre, Model_Extent, land_clip, "", "", "", "MAINTAIN_EXTENT")
    
        my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Begin loop through points")
        
        # Loop through points, then loop through frequencies (as applicable)
        # Loop through points
        point_counter = 1
        for point_id in point_lst:
            #point_id = point_lst[0] #for testing purposes

            # Set conditions under which to display messages
            display_messages = setup_message_display(point_counter)
    
            point_dir = pdir + "pt%s/" % point_id
        
            # Patch for problem where sometimes an existing source file crashes the code despite overwrite == T being on. Trying to manually delete it to see if this fixes the problem
            if arcpy.Exists(Sound_Source):
                arcpy.Delete_management(Sound_Source)
        
            arcpy.Select_analysis(point_source_file, Sound_Source, source_id_field + " = " + point_id)
    
            # Check that one point and only one point was selected
            check_selection(Sound_Source)    

            if model == "nmsimgis" or in_freq[-4:] == ".src" or model == "iso9613":
                # Read in attributes of Sound_Source
                #**# Consider adding an option to allow source file to vary by point
                source_info, source_offset = read_SoundSource(in_freq, Sound_Source, allow_defaults)
                in_file = in_freq
    
                # If model is spreadgis, but in_freq is a .src file, then create a .csv input for SPreAD-GIS
                if model == "spreadgis":
                    allowed_frequencies_only = 1
                    # Re-define in_freq to be the new csv to be created
                    in_file = results_dir + "%s_pt%s.csv" % (os.path.basename(in_freq)[:-4], point_id) # put new source .csv in the results directory
                    use_nmsimgis_source(in_file, source_info, source_offset, tbx_root, allowed_frequencies_only)    

                # Extract required info for later model steps 
                freq_lst, sound_level_lst, measurement_distance_lst, min_source_height, value_weighting = extract_frequencies(in_file)

    
            if display_messages == 1:
                arcpy.AddMessage("Running model for point %s... model to run is: %s" % (point_id, model))
    
            if model == "spreadgis": # or model == "all":
                import spreadgishlpr #**# This should be moved to the top
    
                my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Begin SPreAD-GIS for point %s" % point_id)

                s_offset_exists = 'source_offset' in locals() or 'source_offset' in globals()
                if s_offset_exists == False:
                    source_offset = 0.34 # use a default of 1 if it does not exist
                    if use_old_barrier == 0:
                        arcpy.AddMessage("A source offset of %s was applied to the data to minimize compuational errors" % source_offset)
                        arcpy.AddMessage("Unfortunately, the SPreAD-GIS code does not support any source offsets specified in the points shapefile")
                        arcpy.AddMessage("Unless an NMSim .src file is used as the input")
    
                #Main function call
                # Changed to work with land_clip instead of landcover
                my_times, my_time_labels = spreadgishlpr.spreadgis(point_dir, freq_dir, Sound_Source, point_id,
                          n_points, Model_Extent, freq_lst, sound_level_lst, measurement_distance_lst,
                          dem_ft, land_clip, idir, CellSize, source_offset, receiver_offset, temp, rh, wind_dir, wind_sp, seas_cond, tbx_root,
                          point_fill, freq_fill, my_times, my_time_labels, point_counter, keep_intermediates, use_old_barrier)
    
                my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "SPreAD-GIS model completed for point %s" % point_id)
    
    
            if model == "nmsimgis" or model == "iso9613": # or model == "all":
    
                # Check that source height is valid
                if min_source_height > source_offset:
                    raise ValueError("nmsim .src files specify a minimum source height.\n" +
                    "Your selected source height (%s)" % source_offset +
                    "must be greater than the minimum source height (%s)." % min_source_height)

                make(point_dir) # create directory
    
                '''for troubleshooting
                freq = freq_lst
                reload(nmsimhlpr)
                import nmsimhlpr
                from nmsimhlpr import * 
                '''
        
                my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Begin nmsimgis for point %s" % point_id)
    
                #Main nmsimgis function call
                '''
                #**# Turned multiprocessing off - it was adding its own set of problems
                # Wrapped in multiprocessing to prevent a memory error. Now using subprocess.call to bypass the memory error
                multiprocessing.freeze_support() # Windows-specific code to help prevent crashes
                nmsim_tasks = multiprocessing.Queue()
                nmsim_results = multiprocessing.Queue()

                num_laborers = 1 #**# If you move this part of the code OUTSIDE the loop, you could increase CPU's and use multiprocessing here (and for SPreAD-GIS)
                laborers= [ Laborer(nmsim_tasks, nmsim_results)
                            for i in xrange(num_laborers) ]
                             
                for L in laborers:
                    L.start()
                nmsim_tasks.put(Task_nmsimgis(model_dir, point_dir, freq_dir, Sound_Source, dem_clip, land_clip, Model_Extent, source_offset, receiver_offset, source_info,
                                  temp, rh, freq_lst, point_id, point_fill, freq_fill, tbx_root, my_times, my_time_labels, keep_intermediates, n_processes))
                for k in xrange(num_laborers):
                    nmsim_tasks.put(None)
                                
                my_times, my_time_labels = nmsim_results.get() #**# Need to change if more than 1 laborer - then have an extraction loop, but also move OUTSIDE points loop
                '''
                import nmsimhlpr #**# This should be moved to the top
                my_times, my_time_labels = nmsimhlpr.nmsimgis_setup(model_dir, point_dir, freq_dir, Sound_Source, dem_clip, land_clip, Model_Extent, source_offset, receiver_offset, source_info,
                                  temp, rh, freq_lst, point_id, point_fill, freq_fill, tbx_root, my_times, my_time_labels, keep_intermediates, n_processes, is_GUI, model)
                #arcpy.LoadSettings(settings_file) # restore original settings
    
                my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Completed %s for point %s" % (model, point_id))
    
    
            if model == "hybrid" or model == "nmsim.exe":
                raise ValueError("%s model is not yet supported." % model)

            # Compare propagation dB values to ambient dB values for each frequency
            # .tif or .tiff (more generally .xxx or .xxxx will process only at the end)
            # If ambient path is NA, skip this step
            if ambient_path != "NA":
                
                # Only run if intermediate files are desired
                if keep_intermediates == 1:
                    excess_freq_dir = model_dir + "frequency_excess/"
                    make(excess_freq_dir)
                    compute_excess_by_freqs(point_dir, excess_freq_dir, freq_lst, freq_fill, ambient_path, model, point_id)

            # Add propagation speed as an output
            arcpy.AddMessage(propagation_speed_output)
            if propagation_speed_output == 1:
                # Resetting here, because nmsimgis does not appear to be honoring them
                #arcpy.AddMessage(Model_Extent)
                # I had this try/except disabled, as it appeared to be working, but now someone is getting an error here.
                # I have no idea why ArcGIS is unable to set the extent. Does the Model_Extent variable get redefined at some point?                
                try:
                    arcpy.env.extent = Model_Extent #RasterExtent
                except:
                    arcpy.AddMessage("Could not set extent %s for propagation speed raster for unknown reasons" % Model_Extent)
                try:
                    arcpy.env.cellSize = CellSize
                except:
                    arcpy.AddMessage("Could not set cellsize %s for unknown reasons" % CellSize)
                try:
                    arcpy.env.snapRaster = arcpy.Raster(elevation)
                except:
                    arcpy.AddMessage("Could not set snapRaster for unknown reasons")
                create_propagation_time_raster(Sound_Source, temp, point_dir, model_dir, point_id, point_fill)
                
            # Delete point folders to keep from generating 2,000,000 files!
            if keep_intermediates == 0:
                shutil.rmtree(point_dir, ignore_errors=True) # ignore_errors = True is a temporary patch because not everything can be deleted due to an ArcGIS schema lock
                
            point_counter += 1
        
    # Analysis outputs:
    my_times, my_time_labels = create_summary(point_lst, freq_dir, results_dir, results_label, freq_lst, summarize, weighting, point_fill, freq_fill, timelog, truncate, my_times, my_time_labels, outpath = model_dir, keep_intermediates = keep_intermediates)

    # Set up for ambient analysis
    # Check if ambient is a raster for a single frequency or for an overall raster
    ambient_path, ambient_raster = prep_ambient(ambient_path, freq_lst, freq_fill)

    # Compare summary to overall dB levels
    if ambient_path != "NA" or ambient_raster != "NA":
        compute_excess(freq_lst, point_lst, ambient_path, ambient_raster, weighting, results_dir, results_label, freq_fill, point_fill, freq_dir, summarize)

    # Delete temporary workspace
    if keep_intermediates == 0:
        shutil.rmtree(pdir, ignore_errors=True) # Delete entire points directory (should be empty - point directories are deleted as the loop progresses)
        shutil.rmtree(idir, ignore_errors=True)

    # Print timing of results to timing log file
    #**# Wrap code in a try/except to allow the extraction of timing in the event of an error?
    my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "%s run completed" % model)
    output_timing(timelog, my_times, my_time_labels, delete_existing_timelog)

    if weighting[0] == "W":
        arcpy.AddWarning("Species-specific results are still being developed and vetted. Please do not use these for actual scientific purposes at this point.")        

# Main function call to create ambient sound conditions. Moved from A_AmbientSoundConditions.py on 2016-04-11
def CreateAmbient(in_freq, veg, in_table, wind_speed, snow, ambient_path = "C:/smt/source_data/", freq_fill = 5): 
 
    # Harrison et al. 1980 provides a table with defaults, which I loosely interpreted to create a .csv in the toolbox/Tables folder
                     
    # Import system modules
    #import arcpy, arcpy.sa
    
    # Delete intermediate datasets, if they exist
    arcpy.env.overwriteOutput = True
    
    #freq_only = 1 tells the downstream code not to check for dB field or measurement_distance fields     
    freq_lst = extract_frequencies(in_freq, freq_only = 1)[0] # extract_frequencies returns 5 objects. We just want the first & fifth. 
        
    # Loop through & create ambient conditions for each frequency
    for freq_s in freq_lst:

        freq_lbl = str(int(freq_s)).zfill(freq_fill)

        #give_error = 0
        #error_lst = []
        #try:
        # Add weighting option & reclassification option?
        BAR_sl, CON_sl, HEB_sl, HWD_sl, SHB_sl, URB_sl, WAT_sl = extract_ambient(in_table, freq_s, wind_speed, snow)
    
        # Reclassify land cover dataset by ambient sound level
        ambient = ambient_path + "ambient%s" % freq_lbl
        expression = "BAR " + BAR_sl + ";CON " + CON_sl + ";HEB " + HEB_sl + ";HWD " + HWD_sl + ";SHB " + SHB_sl + ";URB " + URB_sl + ";WAT " + WAT_sl
        reclass = arcpy.sa.Reclassify(veg, "SPREADTYPE", expression, "DATA")
        reclass.save(ambient)
        #except:
        '''
        import sys            
        give_error = 1
        error_lst.append(sys.exc_info()[0:1])  
            
        if give_error == 1:
            arcpy.AddWarning("There was one or more errors with creating ambient background conditions. Ambient backgrounds were not created for all frequencies")
            for error in error_lst:
                arcpy.AddMessage("%s\n%s\n" % (error[0], error[1]))
        '''

'''
#**# Formerly used when there was multiprocessing in the code
# copied from Worker class of nmsimhlpr.
class Laborer(multiprocessing.Process):
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        
    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            # kill process when no more tasks (Is this really necessary? Shouldn't this just happen automatically?)
            if next_task is None:
                print "%s: Exiting" % proc_name
                break
            answer = next_task() 
            self.result_queue.put(answer) #Put answer into results queue
        return
    
class Task_nmsimgis(object):
    def __init__(self, model_dir, point_dir, freq_dir, Sound_Source, dem_clip, land_clip, Model_Extent, source_offset,
                 receiver_offset, source_info, temp, rh, freq_lst, point_id,
                 point_fill, freq_fill, tbx_root, my_times, my_time_labels,
                 keep_intermediates, n_processes):
        self.model_dir = model_dir
        self.point_dir = point_dir
        self.freq_dir = freq_dir
        self.Sound_Source = Sound_Source
        self.dem_clip = dem_clip
        self.land_clip = land_clip
        self.Model_Extent = Model_Extent
        self.source_offset = source_offset
        self.receiver_offset = receiver_offset
        self.source_info = source_info
        self.temp = temp
        self.rh = rh
        self.freq_lst = freq_lst
        self.point_id = point_id
        self.point_fill = point_fill
        self.freq_fill = freq_fill
        self.tbx_root = tbx_root
        self.my_times = my_times
        self.my_time_labels = my_time_labels
        self.keep_intermediates = keep_intermediates
        self.n_processes = n_processes
        
    def __call__(self):
        try:
            import nmsimhlpr
            my_times, my_time_labels = nmsimhlpr.nmsimgis_setup(self.model_dir, self.point_dir, self.freq_dir, self.Sound_Source, self.dem_clip, self.land_clip, self.Model_Extent, self.source_offset,
                                                          self.receiver_offset, self.source_info, self.temp, self.rh, self.freq_lst,
                                                          self.point_id, self.point_fill, self.freq_fill, self.tbx_root, self.my_times, self.my_time_labels,
                                                          self.keep_intermediates, self.n_processes)
        except:
            arcpy.AddMessage("An error occurred: %s\n%s" % (sys.exc_info()[0], sys.exc_info()[1]))
        
        # Return information on code timing (perhaps a terrible idea? This could get confusing in a multiprocessing context)
        return(my_times, my_time_labels)
'''     

# Summarize over model runs
def create_summary(point_lst, path_lst_input, results_dir, results_label, freq_lst, summarize, weighting, point_fill, freq_fill, timelog, truncate, my_times, my_time_labels, outpath = "none", keep_intermediates = 1):
   
   #Set stuff up, then call the existing code.
    # 0: do not do any summarizing
    #if summarize == 0:
        # Nothing happens.

    # 1: summarize over both points and frequencies
    if summarize == 1 or summarize == 3:
        
        path_lst = setup_path_lst(path_lst_input, point_lst)
        main_outpath = setup_outpath(path_lst_input, outpath, summarize)
        points_outpath = main_outpath + "point_summaries/"
        weights_outpath = main_outpath + "weightings/"
        points_temp = points_outpath + "temp/"
        weights_temp = weights_outpath + "temp/"
        sum_temp = main_outpath + "sum_temp/"
        make(points_temp)
        make(weights_temp)
        make(sum_temp)
        
        # Summarize over multiple points
        #e.g. for multiple simultaneous noise sources
        sum_sound_from_multiple_points(freq_lst, point_lst, path_lst, point_fill, freq_fill, points_outpath, points_temp)
        my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Sum Noise from multiple points completed")

        # Re-weight value for each frequency (e.g., A-weighting)
        #weighting = ["A"] # ["A"] for humans; ["Z"] or ["F"] for flat, ["C"] for C-weighting (loud sound sources). ["W",*args] is for wildlife, with *args specifying the species-specific reduction
        # Check that weighting value is allowed
        
        #**# Disable wildlife-specific weighting pending agreement on the correct way to resolve it
        #if weighting[0] != "A" and weighting[0] != "F" and weighting[0] != "Z" and weighting[0] != "C" and weighting[0] != "W":
        #    raise ValueError("Weighting approach %s is not supported by the code. A,C,Z, and W (for wildlife) are the only supported weighting types" % weighting[0] )
        if weighting[0] != "A" and weighting[0] != "F" and weighting[0] != "Z" and weighting[0] != "C":
            raise ValueError("Weighting approach %s is not supported by the code. A, C, and Z are the only supported weighting types" % weighting[0] )

        weight_sound_from_multiple_frequencies(freq_lst, weighting, points_outpath, weights_outpath, weights_temp, freq_fill)
        my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "%s-weighting completed for all frequencies" % weighting[0])

        #print weighting[0]
        #raise ValueError("A-weighting is behaving wonky")
    
        # Sum energy across all frequencies to get overall energy level (should this require all octave bands to be represented?)
        sum_sound_from_multiple_frequencies(freq_lst, weighting, weights_outpath, results_dir, results_label, sum_temp, freq_fill, truncate = truncate)
        my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Sum energy across all frequencies completed")

    # 2: summarize over frequencies for each point separately
    if summarize == 2 or summarize == 4:
        # Loop through each point  
        # Loop through points
        for point_id in point_lst:
            main_outpath = setup_outpath(path_lst_input, outpath, summarize)
            weights_outpath = main_outpath + "weightings/"
            weights_temp = weights_outpath + "temp/"
            sum_temp = main_outpath + "sum_temp%s/" % point_id
            make(weights_temp)
            make(sum_temp)

            # Re-weight value for each frequency (e.g., A-weighting)      
            #weighting = "A" # "A" for humans; "F" for flat, "C" for C-weighting. "W" is for wildlife-specific weightings
            weight_sound_from_multiple_frequencies(freq_lst, weighting, main_outpath, weights_outpath, weights_temp, freq_fill, point_fill, point_id)
            my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "%s-weighting completed for all frequencies" % weighting[0])

            # Sum energy across all frequencies to get overall energy level (should this require all octave bands to be represented?)
            sum_sound_from_multiple_frequencies(freq_lst, weighting, weights_outpath, results_dir, results_label, sum_temp, freq_fill, point_fill, point_id, truncate = truncate)
            my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Sum energy across all frequencies completed")

    # Add summary option to get a maximum across points instead of a sum.
    # No summary across frequencies is made as this option is intended mainly
    # for applications using 1/3 octave results
    if summarize == 5:
        path_lst = setup_path_lst(path_lst_input, point_lst)
        main_outpath = setup_outpath(path_lst_input, outpath, summarize)
        points_outpath = main_outpath + "point_summaries_max/"
        points_temp = points_outpath + "temp/"
        make(points_temp)
        
        # Get maximum from multiple points
        #e.g. for a single noise source traveling a line or area
        max_sound_from_multiple_points(freq_lst, point_lst, path_lst, point_fill, freq_fill, points_outpath, points_temp)
        my_times, my_time_labels = get_time(timelog, my_times, my_time_labels, "Max Noise from multiple points completed")
        

    if keep_intermediates == 0:
        # Weights_outpath and sum_temp are not used for the 5th summary option or for no summary
        if summarize != 5 and summarize != 0:
            shutil.rmtree(weights_outpath, ignore_errors = True) # Delete entire frequency-specific weights directory - do not need those files!
            shutil.rmtree(sum_temp, ignore_errors = True)
        # This directory is only created for summarize 1, 3, or 5 and will throw an error for other summarize options as points_temp will not have been defined.
        if summarize == 1 or summarize == 3 or summarize == 5:
            shutil.rmtree(points_temp, ignore_errors = True)

    return(my_times, my_time_labels)

# ---------------------------------------------------------------------------
# Modified from D_SumPropagationforMultiplePoints.py
# Author: Sarah E. Reed
# Created on: 27 August 2010
# Revised: 18 April 2016 by A.C. Keyel and Jessica Sushinsky
# Description: Sums acoustic energy from multiple simultaneous point 
#              sound sources
# ---------------------------------------------------------------------------
def sum_sound_from_multiple_points(freq_lst, point_lst, path_lst,
                                   point_fill, freq_fill, points_outpath, points_temp):

    arcpy.AddMessage("Summing noise propagation values for %s point sources and %s frequencies..." % (len(point_lst), len(freq_lst)))
   
    for freq_s in freq_lst:        
        # Int. drops any fractions from file names. This assumes that frequency bands are at least 1 Hz apart        
        freq_lbl = str(int(freq_s))
        freq_lbl = freq_lbl.zfill(freq_fill) # Add leading zeros to make file names sort nicely       
        # Convert each input raster from SPL (dB) to acoustic energy (z) and sum them    

        for i in range(len(point_lst)):
            point = str(point_lst[i])
            point = point.zfill(point_fill) # Add leading zeros to make file names sort nicely
            path = path_lst[i]
            this_raster = path + "pt%s_pr%s.tif" % (point, freq_lbl) 
            prev_total = points_temp + "Etot%s%s.tif" % (freq_lbl, i - 1)
            bells_raster = points_temp + "bells%s%s.tif" % (freq_lbl, i)
            acoustic_energy_raster = points_temp + "E%s%s.tif" % (freq_lbl, i)
            # Changes identity at each time step, saving over the previous file crashed ArcGIS/Python            
            total_energy_raster = points_temp + "Etot%s%s.tif" % (freq_lbl, i)
            
            # Convert from dB to acoustic energy
            convert_dB_to_acoustic_energy(this_raster, bells_raster, acoustic_energy_raster)            
            # Add acoustic energy to running total (skip first iteration, no prev raster)
            if i == 0:
                new_total = arcpy.Raster(acoustic_energy_raster)
                new_total.save(total_energy_raster)
            if i != 0:
                new_total = arcpy.Raster(prev_total) + arcpy.Raster(acoustic_energy_raster)
                new_total.save(total_energy_raster)
                
        # Convert acoustic energy (z) back to SPL (dB)
        db_raster = points_outpath + "sumenergy_%s.tif" % freq_lbl
        convert_acoustic_energy_to_dB(total_energy_raster, db_raster)        

# ---------------------------------------------------------------------------
# Modified from D_SumPropagationforMultiplePoints.py and the derivative 
# sum_sound_from_multiple_points
# Authors: Sarah E. Reed, A.C. Keyel, and Jessica Sushinsky
# Created on: 27 August 2010; Last revised 25 April 2017
# Description: Get maximum observed acoustic energy from multiple point 
#              sound sources
# ---------------------------------------------------------------------------
def max_sound_from_multiple_points(freq_lst, point_lst, path_lst,
                                   point_fill, freq_fill, points_outpath, points_temp):

    arcpy.AddMessage("Getting maximum noise propagation values for %s point sources and %s frequencies..." % (len(point_lst), len(freq_lst)))
   
    for freq_s in freq_lst:        
        # Int. drops any fractions from file names. This assumes that frequency bands are at least 1 Hz apart        
        freq_lbl = str(int(freq_s))
        freq_lbl = freq_lbl.zfill(freq_fill) # Add leading zeros to make file names sort nicely       
        # Convert each input raster from SPL (dB) to acoustic energy (z) and sum them    

        for i in range(len(point_lst)):
            point = str(point_lst[i])
            point = point.zfill(point_fill) # Add leading zeros to make file names sort nicely
            path = path_lst[i]
            this_raster = path + "pt%s_pr%s.tif" % (point, freq_lbl) 
            prev_max = points_temp + "Emax%s%s.tif" % (freq_lbl, i - 1)
            bells_raster = points_temp + "bells%s%s.tif" % (freq_lbl, i)
            acoustic_energy_raster = points_temp + "E%s%s.tif" % (freq_lbl, i)
            # Changes identity at each time step, saving over the previous file crashed ArcGIS/Python            
            max_energy_raster = points_temp + "Emax%s%s.tif" % (freq_lbl, i)
            
            # Convert from dB to acoustic energy
            convert_dB_to_acoustic_energy(this_raster, bells_raster, acoustic_energy_raster)            
            # Add acoustic energy to running total (skip first iteration, no prev raster)
            if i == 0:
                new_max = arcpy.Raster(acoustic_energy_raster)
                new_max.save(max_energy_raster)
            if i != 0:
                prev_max_raster = arcpy.Raster(prev_max)
                cur_max = arcpy.Raster(acoustic_energy_raster)
                new_max = arcpy.sa.Con(cur_max > prev_max_raster, cur_max, prev_max_raster) 
                new_max.save(max_energy_raster)
                
        # Convert acoustic energy (z) back to SPL (dB)
        db_raster = points_outpath + "maxenergy_%s.tif" % freq_lbl
        convert_acoustic_energy_to_dB(max_energy_raster, db_raster)        


# Need an input path for each point, but input may be either a single directory of a list of directories.
def setup_path_lst(path_lst_input, point_lst):
    # Check if model_dir is a list. If not, turn it into one
    if type(path_lst_input) is not list:
        path_lst = [path_lst_input] * len(point_lst)
    else:
        path_lst = path_lst_input

    return path_lst

# Simple function to handle diverse inputs & give a correct output path
def setup_outpath(path_lst_input, outpath, summarize):

    if summarize == 1 or summarize == 3 or summarize == 5:
        if outpath == "none":
            main_outpath = path_lst_input # This assumes path_lst_input is a single directory. Otherwise outpath must not be "none"
        else:
            main_outpath = outpath

    if summarize == 2 or summarize == 4:
        if type(path_lst_input) is str:
            main_outpath = path_lst_input
        else:
            main_outpath = outpath

    return main_outpath

# Weight noise from multiple frequencies
def weight_sound_from_multiple_frequencies(freq_lst, weighting, in_dir, weights_outpath, weights_temp, freq_fill, point_fill = "", point_id = "", outpath = "none", is_ambient = 0):

    arcpy.AddMessage("Weighting frequencies to perform %s-weighting" % weighting[0])
  
    #import arcpy    
    point_label = ""
    if point_id != "":
        point_id = point_id.zfill(point_fill) #update local point id to include appropriate number of digits
        point_label = "_pt%s" % point_id.zfill(point_fill)
    
    w_index = 0
    # Loop through each frequency
    for freq_s in freq_lst:
        freq_s = float(freq_s)
        freq_lbl = str(int(freq_s)).zfill(freq_fill)
        if point_id == "":
            summary_raster = in_dir + "sumenergy_%s.tif" % freq_lbl
        else:
            summary_raster = in_dir + "pt%s_pr%s.tif" % (point_id, freq_lbl)
        if is_ambient == 1:
            summary_raster = in_dir + "ambient%s" % freq_lbl
        weighted_raster = weights_outpath + "pr%s%s%s.tif" % (weighting[0], freq_lbl, point_label)
        
        if weighting[0] == "A":
            aweight_raster(freq_s, freq_lbl, summary_raster, weighted_raster, weights_temp)
        # If flat weighting, no weighting needs to be applied, each raster is equal (even though octave bands contain different numbers of frequencies, the energy level is for the entire octave band).
        if weighting[0] == "F" or weighting[0] == "Z":
            # No weighting is applied
            arcpy.CopyRaster_management(summary_raster, weighted_raster)

        if weighting[0] == "C":
            cweight_raster(freq_s, freq_lbl, summary_raster, weighted_raster, weights_temp)
            
        if weighting[0] == "W":
            w_index += 1  # Index immediately increments to 1 to skip the weighting code that is the first entry.
            weight = weighting[w_index]
            wweight_raster(weight, summary_raster, weighted_raster)

# Create a weighting function for a list of values instead of an entire raster
def scalar_weighting(freq_lst, db_lst, weighting):
    if weighting[0] == "A":
        weighted_db = AWeight(freq_lst, db_lst)
        
    # If flat weighting, no weighting needs to be applied, each raster is equal (even though octave bands contain different numbers of frequencies, the energy level is for the entire octave band).
    if weighting[0] == "F" or weighting[0] == "Z":
        # No weighting is applied, sum energy levels from each frequency
        weighted_db = FWeight(db_lst)
        
    if weighting[0] == "C":
        weighted_db = CWeight(freq_lst, db_lst)
        
    if weighting[0] == "W":
        weighted_db = WWeight(weighting, db_lst)
    
    return(weighted_db)
  
  
# Sum sound levels across multiple frequencies
def sum_sound_from_multiple_frequencies(freq_lst, weighting, in_dir, results_dir, results_label, temp_path,
                                        freq_fill, point_fill = "", point_id = "", truncate = 0):
    #import arcpy
    arcpy.AddMessage("Summing %s-weighted sound levels over all frequencies" % weighting[0])
    
    point_label = ""
    if point_id != "":
        point_id = point_id.zfill(point_fill)
        point_label = "_pt%s" % point_id.zfill(point_fill)
    
    # Loop through frequencies
    for i in range(len(freq_lst)):
        freq_lbl = str(int(freq_lst[i])).zfill(freq_fill)

        # Convert from dB to acoustic energy
        this_raster = in_dir + "pr%s%s%s.tif" % (weighting[0], freq_lbl, point_label)
        bells_raster = temp_path + "bells%s" % freq_lbl
        acoustic_energy_raster = temp_path + "E%s" % freq_lbl
        # Changes identity at each time step, I think I was having trouble saving over it            
        total_energy_raster = temp_path + "EFt%s_%s%s.tif" % (freq_lbl, i, point_id)

        convert_dB_to_acoustic_energy(this_raster, bells_raster, acoustic_energy_raster)            

        if i == 0:
            new_total = arcpy.Raster(acoustic_energy_raster)
            new_total.save(total_energy_raster)
        # Sum acoustic energies
        if i != 0:
            prev_freq = int(freq_lst[i-1])
            prev_freq_lbl = str(prev_freq).zfill(freq_fill)
            prev_raster = temp_path + "EFt%s_%s%s.tif" % (prev_freq_lbl, (i - 1), point_id)
            new_total = arcpy.Raster(prev_raster) + arcpy.Raster(acoustic_energy_raster)
            new_total.save(total_energy_raster)

    # Convert acoustic energy (z) back to SPL (dB)
    db_raster_raw = temp_path + "%s_%s%s_raw.tif" % (results_label, weighting[0], point_label)
    convert_acoustic_energy_to_dB(total_energy_raster, db_raster_raw)        

    # Truncate values to some specified minimum (default is 0 dB)
    db_raster = results_dir + "%s_%s%s.tif" % (results_label, weighting[0], point_label)
    truncated_raster = arcpy.sa.Con(arcpy.Raster(db_raster_raw) > truncate, db_raster_raw, truncate)
    truncated_raster.save(db_raster)

# Helper function to convert from dB to acoustic energy
def convert_dB_to_acoustic_energy(dB_raster, bells_raster, acoustic_energy_raster):
    # Convert from dB to Bells
    times_result = arcpy.sa.Times(dB_raster, 0.1)
    times_result.save(bells_raster)
    #sum_energy = arcpy.sa.Raster(Workspace + "/results/sum_energy")
    acoustic_energy = arcpy.sa.Power(10, bells_raster)
    acoustic_energy.save(acoustic_energy_raster)
    
# Helper function to convert from acoustic energy to dB
def convert_acoustic_energy_to_dB(acoustic_energy_raster, SPL_raster):
    calc_result = 10 * arcpy.sa.Log10(acoustic_energy_raster)
    calc_result.save(SPL_raster)   

            
# Create a-weighting
# For now based on AWeight function from Bruce Ikelheimer
def aweight_raster(freq_s, freq_lbl, summary_raster, weighted_raster, temp_w_dir):    
    #  This function, with the input of frequency, F, in Hz, computes
    #  the A-weighting frequency characteristic, Wa in decibels,
    #  according to ANSI S1.4-1983 Appendix C.
    #**# Does this differ from ANSI S1.42-2001 R2016?
    #   The following are factors used in the calculation
        
    k1 = 2.242881e16
    k3 = 1.562339

    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    f4 = 12194.22

    # Calculate Wc term
    term1 = k1 * math.pow(freq_s, 4)
    term2 = math.pow((math.pow(freq_s, 2) + math.pow(f1, 2)), 2)
    term3 = math.pow((math.pow(freq_s, 2) + math.pow(f4, 2)), 2)
    Wc = term1 / (term2 * term3)

    # Then calculate the A-weight
    term4 = math.pow(freq_s, 4)
    term5 = math.pow(freq_s, 2) + math.pow(f2, 2)
    term6 = math.pow(freq_s, 2) + math.pow(f3, 2)
    Wa = k3 * term4 / (term5 * term6)
    Wa = 10 * math.log10(Wa * Wc)

    #Define temporary rasters
    term7_raster_file =  temp_w_dir + "term7_%s" % freq_lbl
    #weighted_energy_raster = model_dir + "temp/E_A_%s" % freq_lbl

    # Delete any existing term7 raster (needed because GUI is not honoring overwrite = TRUE)
    arcpy.Delete_management(term7_raster_file)

    '''
    # Original, but I realized that we are not adding energy here as in Bruce's code
    # So the conversion to energy and back that was introducing the log10(0) problem due to rounding could be bypassed.

    # Compute weighted level
    term7_raster = (Wa + arcpy.Raster(summary_raster)) * 10**-1
    term7_raster.save(term7_raster_file)
    energy = arcpy.sa.Power(10, term7_raster)
    #energy = energy + math.pow(10, (Wa + level[i]) * 10**-1)
    energy.save(weighted_energy_raster)

    dBA = 10 * arcpy.sa.Log10(energy) # This will cancel the arcpy.sa.Power(10, term7_raster) step when there is no addition.
    dBA.save(weighted_raster)
    '''

    term7_raster = (Wa + arcpy.Raster(summary_raster)) * 10**-1
    term7_raster.save(term7_raster_file)
    #energy = arcpy.sa.Power(10, term7_raster)
    #energy = energy + math.pow(10, (Wa + level[i]) * 10**-1)
    #energy.save(weighted_energy_raster)

    dBA = 10 * term7_raster
    dBA.save(weighted_raster)


# Compute C-weighting for a raster
def cweight_raster(freq_s, freq_lbl, summary_raster, weighted_raster, temp_sub_dir):
    
    #  This function, with the input of frequency, F, in Hz, computes
    #  the C-weighting frequency characteristic, Wa in decibels,
    #  according to ANSI S1.4-1983 Appendix C.

    #   Input:
    #   freq - an array of frequencies of interst
    #   level - an array of the input signal level in dB

    #   The following are factors used in the calculation
    k1 = 2.242881e16

    f1 = 20.598997
    f4 = 12194.22

     #Calculate the C-weight
    Wc = k1 * math.pow(freq_s, 4) / (
        math.pow((math.pow(freq_s, 2) + math.pow(f1, 2)), 2) * math.pow((math.pow(freq_s, 2) + math.pow(f4, 2)),
                                                                         2))
    Wc = 10 * math.log(Wc)
    term1_raster_file = temp_sub_dir + "term7_%s" % freq_lbl
    energy_raster = temp_sub_dir + "E_C_%s" % freq_lbl
    term1_raster = (Wc + arcpy.Raster(summary_raster))* 10**-1
    term1_raster.save(term1_raster_file)
    energy = arcpy.sa.Power(10, term1_raster)
    energy.save(energy_raster)

    dBC = 10 * arcpy.sa.Log10(energy)
    dBC.save(weighted_raster)
    return dBC

# Really simple - subtract the weighting from the raster value
# (this could also be done with A & C weighting, based on tables with the weighting changes)
def wweight_raster(weight, summary_raster, weighted_raster):
    dBW = arcpy.Raster(summary_raster) - weight
    dBW.save(weighted_raster)
    return dBW

# Simple version - just one audiogram, set up for analysis, extract everything
def get_w_weights_simple(audiogram):
    
    w_weight = ["W"]   
    weights = {}
    with open(audiogram, 'r') as a:
        count = 0
        for line in a.readlines():
            #Skip header row
            if count == 0:
                count += 1
                continue
        
            line = line[:-1]
            fields = line.split(',')
            # Float ensures the sorting will be correct.
            frequency = float(fields[0])
            db = float(fields[1])
            weights[frequency] = db
            
    for key in sorted(weights.iterkeys()):
        # Add the db values in order of increasing frequency
        w_weight.append(weights[key])

    return w_weight

# Convert an audiogram to a custom weighting approach
def get_w_weights(audiogram, species_lst, target_frequencies, interpolate_missing_frequencies):
    '''
    audiogram = "C:/smt/csu/drafts/auditory_sensitivity_ms/Auditory_sensitivity.csv"
    species_lst = ['Odocoileus virginianus']
    #**# Need to document audiogram format! But maybe use keywords to assign columns?
    '''
          
    # Read in species data #**# Need to extract all frequencies to give option of interpolation.
    db_dict, n_dict = read_species(audiogram, species_lst, target_frequencies)       

    # convert dictionary to a list of frequency values (#**# does sample size matter?)
    db_lst = dict_to_lst(db_dict, target_frequencies, interpolate_missing_frequencies)
    
    # Compute weighting across frequencies
    weighting = make_weighting(db_lst)
    
    return(weighting)

# Convert a drop-down menu to an actual species-specific weighting
def convert_to_weighting(drop_down_weight):
    
    if drop_down_weight == "A":
        weighting = ["A"]
    elif drop_down_weight == "C":
        weighting = ["C"]
    elif drop_down_weight == "Z-flat":
        weighting = ["Z"]
    else:
        raise ValueError("%s weight is not supported at this time" % drop_down_weight)

    return weighting


# Convert a drop-down menu to an actual species-specific weighting
def convert_to_weighting_old(drop_down_weight):
    
    if drop_down_weight == "A-human":
        weighting = ["A"]
    if drop_down_weight == "C-human":
        weighting = ["C"]
    if drop_down_weight == "Z-flat":
        weighting = ["Z"]
    if drop_down_weight == "W-deer":
        raise ValueError("This still needs to be extracted from the Piceance Basin example")
        weights = []
        weighting = ["W", weights]
    #if drop_down_weight == "W-bird":
        raise ValueError("We're sorry, a bird-specific weighting is not scripted yet")

    return weighting

# Read in species data
def read_species(audiogram, species_lst, target_frequencies):
    # Create data dictionaries
    db_dict = {}
    n_dict = {}

    # Read in file
    with open(audiogram, 'r') as a:
        count = 0
        for line in a.readlines():
            line = line[:-1] # strip carriage return
            count += 1
            fields = line.split(",") #**# My use of .csv may cause problems for European colleagues            

            # check header
            if count == 1:
                index = -1
                for field in fields:
                    index += 1
                    if field == "Frequency":
                        freq_index = index
                    if str.lower(field) == "db spl":
                        db_index = index
                    if field == "Species":
                        sci_index = index
                    #if field == "Species (common name)":
                    #    com_index = index
                    if field == "N":
                        n_index = index
                
                continue #Skip rest of loop
                
                
            # Extract data into data dictionary    
            freq = float(fields[freq_index])
            db = float(fields[db_index])
            n = int(fields[n_index])
            
            # Sci name required, eventually accept either scientific or common name, allow switching between the two
            species = fields[sci_index]
            
            sp_ok = 0
            for sp in species_lst:
                if sp == species:
                    sp_ok = 1
                    
            freq_ok = 0
            for f in target_frequencies:
                if f == freq:
                    freq_ok = 1
            
            if sp_ok == 1 and freq_ok == 1:
                # check if frequency is already in dictionary
                if freq in db_dict:
                    # if so, average it
                    old_db = db_dict[freq]
                    old_n = n_dict[freq]
                    new_n = old_n + n
                    new_db = (old_db * old_n + db * n) * new_n**-1
                    db_dict[freq] = new_db
                    n_dict[freq] = new_n
                else:
                    # if not, add it!
                    db_dict[freq] = db
                    n_dict[freq] = n
                    
    return(db_dict, n_dict) 

# Convert my dictionary to an ordered list
def dict_to_lst(db_dict, target_frequencies, interpolate_missing_frequencies):
    
    db_lst = []
    for i in xrange(len(target_frequencies)):
        freq = target_frequencies[i]
        # If frequency exists, add it to the list in order
        if freq in db_dict:
            db_lst.append(float(db_dict[freq]))
        else:
            if interpolate_missing_frequencies == 1:
                raise ValueError("This hasn't been scripted yet")
            else:
                db_lst.append("NA")
    
    return db_lst

# Create a weighting based on an input list of db values
def make_weighting(db_lst):
    #**# This function is not yet ready for public use and requires further consideration
    # Simple solution - subtract minimum value from each value. Then these are the weights to substract from your dB reading.
    min_val = min(db_lst) # min will ignore text when assessing the weighting. it will also ignore numbers stored as text.
    
    db_lst_out = []
    for db in db_lst:
        db_lst_out.append(db - min_val)
    
    return(db_lst_out)

# Define a function for computing A weights
# Function written by B. Ikelheimer
def AWeight(freq, level):
    
    #  This function, with the input of frequency, F, in Hz, computes
    #  the A-weighting frequency characteristic, Wa in decibels,
    #  according to ANSI S1.4-1983 Appendix C.

    #   Input:
    #   freq - an array of frequencies of interst
    #   level - an array of the input signal level in dB

    #   The following are factors used in the calculation
    k1 = 2.242881e16
    k3 = 1.562339

    f1 = 20.598997
    f2 = 107.65265
    f3 = 737.86223
    f4 = 12194.22

    energy = 0

    for i in range(0, len(freq)):
        this_freq = float(freq[i])
        # First calculate the C-weight
        Wc = k1 * math.pow(this_freq, 4) / (
            math.pow((math.pow(this_freq, 2) + math.pow(f1, 2)), 2) * math.pow((math.pow(this_freq, 2) + math.pow(f4, 2)),
                                                                             2))
        # Then calculate the A-weight
        Wa = k3 * math.pow(this_freq, 4) / (
            (math.pow(this_freq, 2) + math.pow(f2, 2)) * (math.pow(this_freq, 2) + math.pow(f3, 2)))
        Wa = 10 * math.log10(Wa * Wc)
        energy = energy + math.pow(10, (Wa + level[i]) / 10)

    dBA = 10 * math.log10(energy)
    return dBA

# Define a function for computing C-weighting
# Function written by B. Ikelheimer
def CWeight(freq, level):
    
    #  This function, with the input of frequency, F, in Hz, computes
    #  the C-weighting frequency characteristic, Wa in decibels,
    #  according to ANSI S1.4-1983 Appendix C.

    #   Input:
    #   freq - an array of frequencies of interst
    #   level - an array of the input signal level in dB

    #   The following are factors used in the calculation
    k1 = 2.242881e16

    f1 = 20.598997
    f4 = 12194.22

    energy = 0

    for i in range(0, len(freq)):
        #Calculate the C-weight
        Wc = k1 * math.pow(freq[i], 4) / (
            math.pow((math.pow(freq[i], 2) + math.pow(f1, 2)), 2) * math.pow((math.pow(freq[i], 2) + math.pow(f4, 2)),
                                                                             2))
        Wc = 10 * math.log(Wc)
        energy = energy + math.pow(10, (Wc + level[i]) / 10)

    dBC = 10 * math.log10(energy)
    return dBC

# Define a function for producing flat-weighted estimates
# Function written by B. Ikelheimer
def FWeight(level):
    
    # This calculates the Flat Weight of a provided signal
    #   Input:
    #   level - an array of the input signal level in dB

    energy = 0

    for i in range(0, len(level)):
        energy = energy + math.pow(10, (level[i]) / 10)

    dBF = 10 * math.log10(energy)
    return dBF

# Do a (wildlife-specific) custom weighting
def WWeight(weighting, db_lst):
    
    
    energy = 0
    for i in range(0, len(db_lst)):
        w_index = i + 1 # first element of weighting is the weighting type, not the weighting to be applied
        weighted_db_level = db_lst[i] - weighting[w_index]
        energy = energy + math.pow(10, weighted_db_level / 10)

    dBW = 10 * math.log10(energy)
    return dBW


# Code to handle frequencies in different inputs
def extract_frequencies(in_freq, freq_only = 0):
    
    value_weighting = "NA"
    min_source_height = "NA"    
    
    # Frequencies may come in as a single frequency, a .src from NMSim, or as
    # a table with frequency, sound level & measurement distance
    if in_freq[-4:] == ".src":
        #extract from NMSim source
        freq_lst, sound_level_lst, measurement_distance_lst, min_source_height = get_nmsim_source(in_freq)
    elif in_freq[-4:] == ".csv":
        out_info = read_table(in_freq, ".csv", freq_only)
        freq_lst = out_info[0]
        sound_level_lst = out_info[1]
        measurement_distance_lst = out_info[2] #[in_Measurement_Distance] * len(freq_lst)
        value_weighting = out_info[4]
    elif in_freq[-4:] == ".txt":
        out_info = read_table(in_freq, ".txt", freq_only)        
        freq_lst = out_info[0]
        sound_level_lst = out_info[1]
        measurement_distance_lst = out_info[2] #[in_Measurement_Distance] * len(freq_lst)
        value_weighting  = out_info[4]
    elif in_freq[-4:] == ".dbf":
        freq_lst, sound_level_lst, in_Measurement_Distance = read_dbf(in_freq)
    else:        
        parts = in_freq.split(';')
        in_freq = parts[0]
        in_Sound_Level = parts[1]
        in_Measurement_Distance = parts[2]
        freq_lst = [in_freq]
        sound_level_lst = [in_Sound_Level]
        measurement_distance_lst = [in_Measurement_Distance]
    
    # If extracting all fields
    if freq_only == 0:
        #Check that lists are all same length
        if len(freq_lst) != len(sound_level_lst):
            raise ValueError("Number of frequencies (%s) must match number of sound levels (%s)" % (len(freq_lst), len(sound_level_lst)))
        
        if len(freq_lst) != len(measurement_distance_lst):
            raise ValueError("Number of frequencies (%s) must match number of measurement distances (%s)" % (len(freq_lst), len(measurement_distance_lst)))
        
    return(freq_lst, sound_level_lst, measurement_distance_lst, min_source_height, value_weighting)        

# Function to extract ambient values from a table
def extract_ambient(in_table, freq_s, wind, snow):
    
    #import arcpy    
    
    # Convert to required string format
    wind = str(wind)
    snow = str(snow)    
    
    if str.lower(in_table[-4:]) == ".csv":
        sep = ","
    elif str.lower(in_table[-4:]) == ".txt":
        sep = "\t"
    else:
        raise ValueError("We're sorry, support for file extension %s has not yet been scripted" % in_table[-4:])
        
    # Create a data dictionary to hold results
    ambient_dict = {}        
        
    # Read table into Python
    table = open(in_table, 'r')

    count = 0
    for line in table.readlines():
        count += 1
        line = line[:-1] #Strip carriage return
        fields = line.split(sep)
        if count == 1:
            # Extract headers
            windmin_index = ambient_fields_setup(wind, "windmin", fields)
            windmax_index = ambient_fields_setup(wind, "windmax", fields)
            snow_index = ambient_fields_setup(snow, "snow", fields)
            continue # skip rest of loop for header row        
        
        freq = fields[0]
        cover = fields[1]
        value = fields[2]

        # Only define optional fields if values are requested.
        windmin = windmax = "#"        
        if wind != "#":
            # If the user set wind to something other than #, but did not include the required columns, give an informative error message
            if windmin_index == "#":
                raise ValueError("Ambient dataset did not include a recognized 'windmin' column. Please add the desired column or set wind = #")                
            if windmax_index == "#":
                raise ValueError("Ambient dataset did not include a recognized 'windmax' column. Please add the desired column or set wind = #")
                
            windmin = fields[windmin_index]
            windmax = fields[windmax_index]
        snowval = "#" # Set default of no value
        if snow != "#":
            # If user forgot to include a snow column, or did not set snow to '#' but meant to, give an informative error message.
            if fields[snow_index] == "#":
                raise ValueError("Ambient dataset did not include a recognized 'snow' column. Please add the desired column or set snow = #")
            snowval = fields[snow_index]
        
        # if criteria match, assign value to dictionary
        if str(freq) == str(freq_s) and snow == snowval:
            if windmin == "#":
                ambient_dict[cover] = value 
            if windmin != "#":
                if float(wind) >= float(windmin) and float(wind) < float(windmax):
                    ambient_dict[cover] = value

    # Check that something was extracted, if not, raise error
    if len(ambient_dict) == 0:
        raise ValueError("Something went wrong with table extraction. Perhaps"
                         "your desired frequency is not in the table?\n"
                         "Frequency = %s, wind = %s, snow = %s" % (freq_s, wind, snow))
        
    # Fails on windspeed of 16. Getting string & int objects!  Same error with snow == 0. Probably need a more informative error message
    # Extract values from dictionary
    BAR_sl = ambient_dict_extract("BAR", ambient_dict)
    CON_sl = ambient_dict_extract("CON", ambient_dict)
    HEB_sl = ambient_dict_extract("HEB", ambient_dict)
    HWD_sl = ambient_dict_extract("HWD", ambient_dict)
    SHB_sl = ambient_dict_extract("SHB", ambient_dict)
    URB_sl = ambient_dict_extract("URB", ambient_dict)
    WAT_sl = ambient_dict_extract("WAT", ambient_dict)

    if len(ambient_dict) < 7:
        raise ValueError("Something went wrong, not all land cover types received" +
                        "a background value:\n BAR %s\nCON: %s\n" % (BAR_sl, CON_sl) +
                        "HEB: %s\nHWD: %s\nSHB: %s\nURB: %s\nWAT: %s" % (HEB_sl, HWD_sl, SHB_sl, URB_sl, WAT_sl))

    if len(ambient_dict) > 7:
        raise ValueError("Something went wrong, somehow you obtained more valus than there are landcover types.")
        
    return(BAR_sl, CON_sl, HEB_sl, HWD_sl, SHB_sl, URB_sl, WAT_sl)

# Optional fields setup
def ambient_fields_setup(variable, field_name, fields):
    # Check that variable is actually included as an input    
    index = "#"    
    if variable != "#":
        # Loop through fields
        for i in range(0,len(fields)):
            this_field = fields[i]
            # If the field matches the field name, create an index for field retrieval. This will support columns in any order and missing columns
            if str.lower(this_field) == str.lower(field_name):
                index = i
    return(index)

# Extract values and if not present, give a no-data value
def ambient_dict_extract(code, ambient_dict):
    if code in ambient_dict:
        out_val = ambient_dict[code]
    else:
        out_val = -9999        
        
    return out_val
    

# Assess area where points would be audible using a user-defined audibility threshold
# Started with a copy of AssessImpacts function and modified it for dprime calculation
def AssessAreaAudible(output_file, point_source_file, source_id_field, n_points, dprime_threshold, focal_area, dprime_path, results_label):
    # Set up a temporary directory for intermediate files
    temp_dir = dprime_path + "calcs/" # Variable name is a code artifact
    make(temp_dir)    
    
    # set up inputs for multipoints
    point_lst = get_points(point_source_file, source_id_field, n_points)

    # Add constants for adding leading 0's, but not more than 3 to avoid problems with file names being too long: 13 char max
    point_fill = len(str(len(point_lst))) # Get number of digits. len of str is number of characters. len of pointlist gets number of elements. So it's the number of digits in len(point_lst)     
    point_fill = min(point_fill, 3)
       
    area_above_threshold_lst = []

    output_vec = ["AREA_ABOV"]
    
    # Loop through points
    for point_id in point_lst:
        #print point_id
        # Find results files that correspond to this point
        point_label = "_pt%s" % point_id.zfill(point_fill)
        dprime_raster = dprime_path + "%s%s.tif" % (results_label, point_label)

        ## Compute impact on focal area from this point    
        # Compute average sound impact
        focal_area_field = "FID"
        
        aat_out = compute_area_above_threshold(focal_area, focal_area_field, dprime_raster, dprime_threshold, temp_dir)        
        area_above_threshold_lst.append(aat_out)
        
    # Copy source file
    arcpy.CopyFeatures_management(point_source_file, output_file)

    # Add fields to hold summary info    
    for item in output_vec:
        arcpy.AddField_management (output_file, item, "FLOAT", 16, 4) #16 is precision, 4 is scale.

    # use update cursor to update point file fields
    id_check = []
    with arcpy.da.UpdateCursor(output_file, [source_id_field, "AREA_ABOV"]) as rows:
        count = 0
        for row in rows:

            # Check that rows are updated in the correct order. I would assume so, but might need a more explicit check for that!
            id_check.append(str(row[0]))
            row[1] = area_above_threshold_lst[count]
            rows.updateRow(row)            
            count += 1
            
            # Stop once the number of specified points has been reached
            if n_points != "all":            
                if count == n_points:
                    break

    if id_check != point_lst:
        raise ValueError("Rows were not updated in the correct order, your updated",
                         "file is corrupted")

    # Remove directory with temporary calculations
    shutil.rmtree(temp_dir)

# Assess length of a line where a noise source or sources would be audible
def AssessLengthAudible(line, line_final, dprime_raster, temp_dir, threshold = 7.3):
    '''
    line A line shapefile where an audibility analysis is desired
    line_final The final shape file containing the length and percent of the shapefile where noise would be audible
    dprime_raster A raster layer containing the dprime statistic
    temp_dir a directory to contain intermediate files
    threshold The threshold of dprime to use for determining audibility
    '''
    
    
    # Convert the dprime raster to one that can have an attribute table
    #dprime_32bit = temp_dir + "dprime32bit.tif"
    #arcpy.CopyRaster_management(dprime_raster, dprime_32bit, "", "", "","", "", "32_BIT_FLOAT")    
        
    #**# May need to calculate the attribute table. But I think the VALUE field ought to exist        
        
    # Reclassify dprime raster according to threshold
    rc_file = temp_dir + "rc_dprime.tif"
    dprime_raster = arcpy.Raster(dprime_raster)
    rc_dprime = arcpy.sa.Con(dprime_raster >= threshold, 1, 0)
    rc_dprime.save(rc_file)
    #mapping = arcpy.sa.RemapRange([])
    #reclass = arcpy.sa.Reclassify(dprime_32bit, "VALUE", expression, "DATA")
    
    # Convert dprime raster to polygon
    dprime_poly = temp_dir + "rc_poly.shp"
    arcpy.RasterToPolygon_conversion(rc_dprime, dprime_poly, "NO_SIMPLIFY")        
    
    # Intersect polygon and line & associate attributes of the polygon with the line segments
    line2 = temp_dir + "line2.shp"
    arcpy.Intersect_analysis([dprime_poly, line], line2, "ALL", 0, "LINE")    
    
    # Calculate length (km) of each line segment
    arcpy.AddField_management(line2, "LENGTH", "DOUBLE")    
    arcpy.CalculateField_management(line2, "LENGTH", "!shape.length@kilometers!", "PYTHON")

    # Calculate length (km) audible of each line segment
    arcpy.AddField_management(line2, "L_AUDIBLE", "DOUBLE")
    arcpy.CalculateField_management(line2, "L_AUDIBLE", '!LENGTH! * !GRIDCODE!', "PYTHON")    
    
    # Sum up the total audible and total length of the line
    arcpy.Dissolve_management(line2, line_final,"", [["LENGTH", "SUM"], ["L_AUDIBLE", "SUM"]])
    
    # Calculate the percent of the line where the noise would be audible
    arcpy.AddField_management(line_final, "P_AUDIBLE", "DOUBLE")
    arcpy.CalculateField_management(line_final, "P_AUDIBLE", '!SUM_L_AUDI! * !SUM_LENGTH!**-1', "PYTHON")    
    

# Assess impact of individual source points on an area of interest
# Uses data summed across all frequencies for the analysis
def AssessImpact(output_file, point_source_file, source_id_field, n_points, focal_area,
                 results_dir, results_label_base, threshold):
    '''
    output_file is the name and path of the resulting point file with the metrics
    point_source_file contains the point data,
    source_id_field is a unique point id,
    n_points is the number of points to run the analysis for,
    focal_area is the area of interest,
    results_dir is the directory containing model analysis outputs,
    results_label_base is the pattern for results file, omitting the point numbers
    threshold is the threshold for computing area impacted above a set sound level threshold.
    '''
    
    # Set up a temporary directory for intermediate files
    analysis_dir = results_dir + "calcs/" # Variable name is a code artifact
    make(analysis_dir)    
    
    # set up inputs for multipoints
    point_lst = get_points(point_source_file, source_id_field, n_points)

    # Add constants for adding leading 0's, but not more than 3 to avoid problems with file names being too long: 13 char max
    point_fill = len(str(len(point_lst))) # Get number of digits. len of str is number of characters. len of pointlist gets number of elements. So it's the number of digits in len(point_lst)     
    point_fill = min(point_fill, 3)
       
    mean_lst = []
    max_lst = []
    min_lst = []    
    area_above_threshold_lst = []
    percent_above_threshold_lst = []
    area_lst = []

    output_vec = ["MEAN", "MAX", "MIN", "AREA_ABOV", "AREA_PERC"]
    
    # Loop through points
    for point_id in point_lst:
        print point_id
        # Find results files that correspond to this point
        point_label = "_pt%s" % point_id.zfill(point_fill)
        db_levels = results_dir + "%s%s.tif" % (results_label_base, point_label)

        ## Compute impact on focal area from this point    
        # Compute average sound impact
        focal_area_field = "FID" #**# The .tif doesn't have an attribute table, so this step might crash
        mean_out, max_out, min_out, area_out = compute_average_sound_contribution(focal_area, focal_area_field, db_levels, analysis_dir)
        mean_lst.append(mean_out)
        max_lst.append(max_out)
        min_lst.append(min_out)
        area_lst.append(area_out)

        # Compute area impacted above a threshold # Default is 55 dB # max above threshold is not different than max (assuming max > threshold)
        aat_out = compute_area_above_threshold(focal_area, focal_area_field, db_levels, threshold, analysis_dir)
        area_above_threshold_lst.append(aat_out)
        percent_above_threshold = aat_out * area_out**-1 * 100 # Needs to be multiplied by 100 to be a percent
        percent_above_threshold_lst.append(percent_above_threshold) 


    # Copy source file
    # Get just the source name    
    #source_parts = point_source_file.split("/") 
    #if len(source_parts) == 1:
    #    source_parts = point_source_file.split("\\")
    #source_base_name = source_parts[-1] # Last element is the file name
    #updated_source_file = results_dir + source_base_name[:-4] + "_updated" + source_base_name[-4:] # keep the same name, insert updated into it before the file ending
    arcpy.CopyFeatures_management(point_source_file, output_file)

    # Add fields to hold summary info    
    for item in output_vec:
        arcpy.AddField_management (output_file, item, "FLOAT", 16, 4) #16 is precision, 4 is scale.

    # use update cursor to update point file fields
    id_check = []
    with arcpy.da.UpdateCursor(output_file, [source_id_field, "MEAN", "MAX", "MIN", "AREA_ABOV", "AREA_PERC"]) as rows:
        count = 0
        for row in rows:

            # Check that rows are updated in the correct order. I would assume so, but might need a more explicit check for that!
            id_check.append(str(row[0]))
            row[1] = mean_lst[count]
            row[2] = max_lst[count]
            row[3] = min_lst[count]
            row[4] = area_above_threshold_lst[count]
            row[5] = percent_above_threshold_lst[count]
            rows.updateRow(row)            
            count += 1
            
            # Stop once the number of specified points has been reached
            if n_points != "all":            
                if count == n_points:
                    break

    if id_check != point_lst:
        raise ValueError("Rows were not updated in the correct order, your updated",
                         "file is corrupted")

    # Remove directory with temporary calculations
    shutil.rmtree(analysis_dir)
   
def compute_average_sound_contribution(focal_area, focal_area_field, db_levels, analysis_dir):

    zst_table = analysis_dir + "Zonal_Statistics.dbf" 
    arcpy.sa.ZonalStatisticsAsTable (focal_area, focal_area_field, db_levels, zst_table, "DATA", "ALL")
    
    # Use search cursor to extract mean & max #**# OLDER SEARCH CURSOR SYNTAX MAY LEAVE PERSISTENT SCHEMA LOCKS
    srows = arcpy.SearchCursor(zst_table)
    count = 0
    for row in srows:
        if count == 1:
            raise ValueError("Only a single focal area is supported. Please dissolve your focal area such that there is only one value in the focal_area_field varaible")

        # Extract statisticsaverage
        mean_out = row.getValue("MEAN") # Arithmetic mean, not reflective of average energy!
        max_out = row.getValue("MAX")
        min_out = row.getValue("MIN")
        total_area = row.getValue("AREA")
        #**# any value to the standard deviation STD?

        count += 1

    # Delete search cursor objects
    del srows, row

    return(mean_out, max_out, min_out, total_area)
   
   
# Compute area impacted above a threshold
def compute_area_above_threshold(focal_area, focal_area_field, db_levels, threshold, analysis_dir):

    # Reclassify all values below threshold to NoData (nodata will be ignored in the summary statistics)
    rc_file = analysis_dir + "rc_db_levels"
    #my_remap = arcpy.sa.RemapRange([[0,threshold,"NODATA"]]) #**# if the existing doesn't work, substitute this!
    rc_db_levels = arcpy.sa.Reclassify(db_levels, "VALUE", "0 %s NODATA" % threshold, "DATA")
    rc_db_levels.save(rc_file)

    # Convert reclassification to polygon
    rc_poly = analysis_dir + "rc_poly.shp"
    arcpy.RasterToPolygon_conversion(rc_file, rc_poly, "NO_SIMPLIFY")

    # Dissolve polygons into a single polygon
    dissolve_poly = analysis_dir + "dissolve_poly.shp"
    arcpy.Dissolve_management(rc_poly, dissolve_poly)

    # Intersect polygon with buffer area
    intersect_poly = analysis_dir + "intersect_poly.shp"
    arcpy.Intersect_analysis([dissolve_poly, focal_area], intersect_poly)    

    # Calculate area of intersection to determine area affected
    arcpy.AddField_management(intersect_poly, "AREA", "DOUBLE")
    arcpy.CalculateField_management(intersect_poly, "AREA", "!shape.area@squaremeters!", "PYTHON") 

    #zst_table = analysis_dir + "Zonal_Statistics_aat.dbf"
    #arcpy.sa.ZonalStatisticsAsTable (focal_area, focal_area_field, rc_db_levels, zst_table, "DATA", "ALL")
        
    # Use search cursor to extract area
    with arcpy.da.SearchCursor(intersect_poly, "AREA") as srows:
        count = 0
        for row in srows:
            if count == 1:
                raise ValueError("Only a single focal area is supported. Please dissolve your focal area such that there is only one value in the focal_area_field varaible")
        
            # Extract statisticsaverage
            area_out = row[0]
    
            count += 1
    
     # Check if there are any rows in the table, or if everything is NoData
    if count == 0:
        # Assign an area of 0 if all is NoData
        area_out = 0
        
    return(area_out)

# Compute average excess above ambient sound impact
#**# IN PROGRESS
def compute_average_excess_contribution(focal_area, ambient, excess_levels, excess_adjustment):
    #**# Requires revisiting the excess calculations
    
    # Subtract ambient conditions from propagation values to calculate excess noise propagation
    ex_file = results_dir + "ex" + freq_s #**# Fix
    ambient = arcpy.sa.Raster(ambient) #**# FIx
    trueValue = pr - ambient 
    ex = arcpy.sa.Con((pr > (ambient + excess_adjustment)), trueValue, 0)
    ex.save(ex_file)
    
    zst_table = analysis_dir + "Zonal_Statistics_ex.dbf" #**# Where to put this?
    arcpy.sa.ZonalStatisticsAsTable (focal_area, focal_area_field, ex, zst_table, "DATA", "ALL")

    # Use search cursor to extract mean & max
    srows = arcpy.SearchCursor(zst_table)
    count = 0
    for row in srows:
        if count == 1:
            raise ValueError("Only a single focal area is supported. Please dissolve your focal area such that there is only one value in the focal_area_field varaible")

        # Extract statisticsaverage
        average_out = row.getValue("MEAN") # Arithmetic mean, not reflective of average energy!
        max_out = row.getValue("MAX")
        min_out = row.getValue("MIN")
        total_area = row.getValue("AREA")
        #**# any value to the standard deviation STD?

        count += 1

    # Delete search cursor objects
    del srows, row

    return(average_out, max_out, min_out, total_area)


# Run Sensitivity analysis to get an estimate of range of sound propagation possibilities
def Sensitivity_Analysis(model, point_source_file, source_id_field, Model_Extent, ambient_path,
                         elevation, landcover, d_max, d_increment, temp_extremes, rh_extremes, 
                         temp_increment, rh_increment, in_wind_sp, seas_cond_lst,
                         base_output_dir, results_dir, in_freq, receiver_offset, n_points, tbx_root,
                         timelog, delete_existing_timelog, keep_intermediates,
                         weighting, truncate):

    ''' for trouble shooting in piceance_basin
    from soundprophlpr import *
    point_source_file = source1
    source_id_field = source1_id
    in_freq = drill_src
    source_info = drill_source_info
    source_offset = drill_offset
    in_Sound_Level = Sound_Level
    in_Measurement_Distance = Measurement_Distance
    base_output_dir = opath
    seas_cond = "none"
    '''
    
    ''' for troubleshooting based on designed landscape
    from soundprophlpr import *
    model = run_name
    point_source_file = sources
    Model_Extent = lc_extent
    in_wind_sp = max_wind_sp
    base_output_dir = output_dir
    in_freq = srcfile
    in_Sound_Level = Sound_Level
    in_Measurement_Distance = Measurement_Distance
    '''    

    raise ValueError("Regretably, the sensitivity analysis option does not work at this time")
    
    #import time    
   
    # Summarize is set to 1. I just want one file to deal with!    
    summarize = 1

    # Set up NMSim dll (doesn't actually matter for SPreAD-GIS)
    nmsim_dll = tbx_root + "dlls/NMSim_Libraries.dll"

    # Get list of frequencies    
    out_lst = extract_frequencies(in_freq)
    freq_lst = out_lst[0]
    spl_lst = out_lst[1]
        
    if model == "nmsimgis":
        # Get list of source levels
        from nmsimhlpr import NMSim_Source
        #**# This will cause problems if the sensitivity analysis is for more than one point.
        #**# And realistically, people will want sensitivity for more than one point.
        source_info, source_offset = read_SoundSource(in_freq, point_source_file)
        head, roll, vel, pitch, engpow, srcfile = source_info
        # Source levels shouldn't depend on distance to source, these only matter for the distance output, which is not use here
        xyzsrc = [0,0,0] # Arbitrary
        xyzrec = [1,1,1] # Arbitrary
        # [2] only keeps the dbband_out output
        dbband_out = NMSim_Source(nmsim_dll, xyzsrc, xyzrec, source_offset,
                                           receiver_offset, head, roll, vel, pitch, engpow,
                                           in_freq, len(freq_lst))[2]
        for i in xrange(len(dbband_out)):
            value = dbband_out[i]
            spl_lst[i] = value

    wind_dir_lst = [0, 45, 90, 135, 180, 225, 270, 315]    
    min_dir = results_dir + "min/"
    max_dir = results_dir + "max/"
    make(min_dir)
    make(max_dir)

    # Screen atmospheric absorption values, find maximum & minimum atmospheric absorption
    temp_dir = base_output_dir + "temporary/"
    make(temp_dir) # create a temporary directory, because ArcGIS can't do some processes in memory like an efficient software
    weather_lst, minmax_lst = atmo_screen_v2(model, freq_lst, spl_lst, temp_extremes, rh_extremes, d_max, d_increment, elevation, point_source_file,
                                           temp_dir, temp_increment, rh_increment,
                                           weighting, nmsim_dll, tbx_root)
                
    #time.sleep(5) # Add a 1 second delay to allow the schema-lock to clear. I hate schema locks! Did not fix problem
    #try:
    shutil.rmtree(temp_dir) # remove temporary directory
    #**# This does not work - there is a schema lock that I can't get rid of! I tried deleting every object that could possibly be causing the schema lock, but it is still there!\
    #**# I think this was fixed - check?
    #except:
    #    arcpy.AddMessage("unable to delete temp_dir")
    for w in xrange(len(weather_lst)): # min & max weather_lsts should have same length
        weather = weather_lst[w]
        minmax = minmax_lst[w] # get id for whether this is a min or max        
        temp = weather.split(';')[0]
        rh = weather.split(';')[1]
    
        # Run for each set of desired seasonal conditions
        loopcount = 0
        for seas_cond in seas_cond_lst:
            loopcount += 1
            # Only run once for nmsimgis (basically bypass the seas_cond_lst loop, as seas_cond is not a nmsimgis input)
            if loopcount > 1 and model == "nmsimgis":
                break
            
            seas_output_dir = base_output_dir + "%s_%s/" % (seas_cond, w)
            make(seas_output_dir)        
            
            # Run model for calm for min & max
            wind_sp = 0
            wind_dir = 0
            
            # Run for max conditions
            output_dir = seas_output_dir + "%s_calm/" % minmax
            make(output_dir)
            #max_dir_lst.append(output_dir)
            results_label = "%s_calm" % minmax

            results_dir = min_dir
            if minmax == "max":
                results_dir = max_dir
            
            SoundPropagation(model, point_source_file, source_id_field, Model_Extent, ambient_path,
                                 elevation, landcover, temp, rh, wind_dir, wind_sp,
                                 seas_cond, output_dir, results_dir, results_label, in_freq, receiver_offset,
                                 n_points, tbx_root, timelog,
                                 delete_existing_timelog, summarize, keep_intermediates, weighting, truncate = truncate)        
                                                             
            # nmsimgis currently does not include wind direction
            if model == "spreadgis":
                # Run model for each wind direction #wind from 8 directions for min & max
                wind_sp = in_wind_sp
                for wind_dir in wind_dir_lst:
                    # Run for max conditions
                    output_dir = seas_output_dir + "%s_wind_dir_%s/" % (minmax, wind_dir)
                    #max_dir_lst.append(output_dir)
                    make(output_dir)
                    results_label = "%s_wind_dir_%s" % (minmax, wind_dir)
    
                    SoundPropagation(model, point_source_file, source_id_field, Model_Extent, ambient_path,
                                         elevation, landcover, temp, rh, wind_dir, wind_sp,
                                         seas_cond, output_dir, results_dir, results_label, in_freq, receiver_offset,
                                         n_points, tbx_root, timelog,
                                         delete_existing_timelog, summarize, keep_intermediates, weighting, truncate = truncate)        
                                         
            
    # Create map of composite maximum propagation conditions (but note that it is a composite of minimums from wind analyses)
    out_raster = base_output_dir + "max_%s" % model
    make_composite_map(model, weighting, out_raster, max_dir, "MAXIMUM")
    
    # Create map of composite minimum propagation conditions (but note that it is a composite of maximum from wind analyses)
    out_raster = base_output_dir + "min_%s" % model
    make_composite_map(model, weighting, out_raster, min_dir, "MINIMUM")

    # Return conditions for maximum and minimum atmospheric absorption
    return(weather_lst, minmax_lst)

def make_composite_map(model, weighting, out_raster, in_dir, min_or_max):
    
    # Get file names using dir_lst
    raster_lst = []
    dir_lst = os.listdir(in_dir)
    for raster in dir_lst:
        # Just run for files with a .tif extension (avoid .tfw and .tif.aux.xml files that are also in the directory)
        if raster[-4:] == ".tif": 
            raster_dir = in_dir + raster
            raster_lst.append(raster_dir)

    # Compute minimum or maximum value
    oraster = arcpy.sa.CellStatistics(raster_lst, min_or_max)
    oraster.save(out_raster)

# Run model for different atmospheric absorption conditions & return min & max
#**# Think about substituting pressure for elevation
#**# For now, there is not a good way to incorporate either - I can screen for
# elevation effects, but I can't change the mean elevation used in the model
# Assumes homogeneous elevation, temperature, and relative humidity over study area
def atmo_screen_v2(model, freq_lst, spl_lst, temp_extremes, rh_extremes, d_max, d_increment, dem,
                   sound_source, temp_dir, temp_increment, rh_increment, weighting, nmsim_dll, tbx_root):

    #import subprocess
    # Get mean elevation for this landscape
    elev_m = atmo_screen_get_mean_elevation(temp_dir, dem, sound_source)
    
    # create a distance vector at which to evaluate effects of distance
    # d_max must be in meters
    dist_vec = seq(0, d_max, d_increment)

    # Create lists to hold outputs
    min_weather = [] * len(dist_vec)
    max_weather = [] * len(dist_vec)
    
    # Have distance as the outer loop, then do the inner loop for each distance
    for i in xrange(len(dist_vec)):
        distance = dist_vec[i]
        wc = atmo_screen_core(distance, freq_lst, spl_lst, weighting, elev_m, temp_extremes, temp_increment, rh_extremes, rh_increment, model, nmsim_dll)
        min_weather.append(wc[0])
        max_weather.append(wc[1])
        
    # Reduce list to minimum number of temp & rh combinations
    min_weather = list(set(min_weather))
    max_weather = list(set(max_weather))

    # compile into a single list for iterating over, with a second list to indicate minimum or maximum
    weather_lst = min_weather + max_weather # + here combines the two lists
    minmax_lst = ["min"] * len(min_weather) + ["max"] * len(max_weather)

    return (weather_lst, minmax_lst)

# Helper function to get mean elevation for atmo_screen_v2
def atmo_screen_get_mean_elevation(temp_dir, dem, sound_source):
    # Get mean elevation for this landscape (can't use code from spreadgishlpr - that is only the mean elevation of the source point)
    mean_raster_file = temp_dir + "mean_elevation.tif"
    mean_raster_table = temp_dir + "mean_elevation.dbf"
    mean_raster = arcpy.sa.CellStatistics([dem], "MEAN", "DATA")
    mean_raster.save(mean_raster_file)

    '''
    hlpr_script = tbx_root + "scripts/sensitivity_schema_lock_bypass.py"
    command = "%s" % hlpr_script # %s %s %s" % (hlpr_script, mean_raster, sound_source, mean_raster_table)
    print command
    elev_m = subprocess.check_output([command])
    '''

    arcpy.sa.Sample(mean_raster, sound_source, mean_raster_table, "NEAREST") # Doesn't need to be the sound source, but that will at least be in the landscape.
    table_fields = arcpy.ListFields(mean_raster_table)
    field_name = "NA"    
    for field in table_fields:
        if str.lower(str(field.name[0:4])) == "mean":
            field_name = field.name
    if field_name == "NA":
        raise ValueError("Could not find the appropriate value field for calculating mean elevation")
    
    with arcpy.da.SearchCursor(mean_raster_table, field_name) as cursor:
        for row in cursor:
            elev_m = row[0]
    
    return elev_m

# Do the actual screen for each distance & SPL
def atmo_screen_core(distance, freq_lst, spl_lst, weighting, elev_m, temp_extremes, temp_increment, rh_extremes, rh_increment, model, nmsim_dll):
    # Set a counter to avoid never ending while loops (e.g., if min & max are switched!)
    #**# value-error here throws away all the work done! However, the values are high enough they should not be triggered

    min_dB = 1000 # Arbitrary high value
    max_dB = -1000 # Arbitrary low value
    # Get maximum temperature, rh   
    temp_max = temp_extremes[1]
    rh_max = rh_extremes[1]
    temp = temp_extremes[0]
    tcount = 0
    while temp < temp_max:
        rh = rh_extremes[0] # Reset humidity at start of each temperature loop
        tcount += 1        
        if tcount == 5000:
            raise ValueError("Maximum temperature count reached %s iterations" % tcount)

        rh_count = 0
        while rh < rh_max:
            rh_count += 1
            if rh_count == 5000:
                raise ValueError("Maximum rh count reached %s iterations" % rh_count)
                
            db_lst = ["NA"] * len(freq_lst) # Get a list to hold absorptions by frequency
            if model == "spreadgis":
                from spreadgishlpr import atmospheric_absorption_loss_core as aal_core

                # Run calculation
                for i in xrange(len(freq_lst)):
                    freq = float(freq_lst[i])
                    temp_k = temp + 273.15
                    absorption = aal_core(elev_m, rh, temp_k, freq)
                    db_lst[i] = float(spl_lst[i]) - float(absorption) * float(distance) # Distance needs to be in m
                    #absorption_lst[i] = absorption * distance 

            if model == "nmsimgis":
                from nmsimhlpr import NMSim_atmospheric_absorption
                temp = float(temp)
                rh = float(rh)
                freq = map(float, freq_lst)
                nfreq = len(freq)
                absorption_lst = NMSim_atmospheric_absorption(nmsim_dll, temp, rh, freq,
                                                                    nfreq, distance)
                for i in xrange(len(absorption_lst)):
                    absorption = absorption_lst[i]
                    db_lst[i] = float(spl_lst[i]) - float(absorption) # Already multiplied by distance
            
            # Summarize absorptions across frequencies with weighting
            dB_level = scalar_weighting(freq_lst, db_lst, weighting)
            # Maximum dB level corresponds to maximum propagation
            if dB_level > max_dB:
                max_dB = dB_level
                temp_at_max = temp                    
                rh_at_max = rh
                
            # Minimum absorption corresonds to maximum propagation
            if dB_level < min_dB:
                min_dB = dB_level
                temp_at_min = temp
                rh_at_min = rh
                    
            rh = rh + rh_increment
        temp = temp + temp_increment

    min_weather = "%s;%s" % (temp_at_min, rh_at_min)
    max_weather = "%s;%s" % (temp_at_max, rh_at_max)
    dBs = [min_dB, max_dB]
    #max_vals = [max_dB, temp_at_max, rh_at_max, elev_m]
    #min_vals = [min_dB, temp_at_min, rh_at_min, elev_m]
     
    return(min_weather, max_weather, dBs, elev_m)           
    #return(max_vals, min_vals)

    
# Extract nmsim source info for SPreAD-GIS 
def use_nmsimgis_source(csv_out, source_info, source_offset, tbx_root, allowed_frequencies_only):
    #import arcpy   
    import os
    import nmsimhlpr
    nmsim_dll = tbx_root + "dlls/NMSim_Libraries.dll"
    head, roll, pitch, vel, engpow, srcfile = source_info
    freqs, sound_levels, measurement_distances, min_source_height, value_weighting = extract_frequencies(srcfile)
    nfreq = len(freqs)    

    # Run once to get refdist, then re-run at that reference distance to get levels
    xyzsrc = [0,0,0]
    xyzrec = [15,0,0] # Take sound measurement at 15 m
    zsrc = source_offset
    zrec = 1 # Assign an arbitrary receiver height of 1 m.
    
    # Run once to get refdist, then re-run at that reference distance to get levels
    refdist_out, distance_out, dbband_out = nmsimhlpr.NMSim_Source(nmsim_dll, xyzsrc, xyzrec, zsrc, zrec, head, roll, vel, pitch,
                 engpow, srcfile, nfreq)
    xyzrec = [refdist_out, 0, 0]
    refdist_out, distance_out, dbband_out = nmsimhlpr.NMSim_Source(nmsim_dll, xyzsrc, xyzrec, zsrc, zrec, head, roll, vel, pitch,
                 engpow, srcfile, nfreq)

    # write outputs to file for reading in to SPreAD-GIS
    #need to delete an existing file
    if os.path.isfile(csv_out):
        os.remove(csv_out)
    with open(csv_out, 'a') as csv:
        # write header
        csv.write("Frequency,dB,Measurement_Distance\n")

        # Extract & write data
        for i in xrange(len(dbband_out)):
            freq = int(freqs[i])
            db = dbband_out[i]
            if allowed_frequencies_only == 1:
                write_csv = 0
                if freq == 125 or freq == 160 or freq == 200 or freq == 200:
                    write_csv = 1
                if freq == 250 or freq == 315 or freq == 400 or freq == 500:
                    write_csv = 1
                if freq == 630 or freq == 800 or freq == 1000 or freq == 1250:
                    write_csv = 1
                if freq == 1600 or freq == 2000:
                    write_csv = 1
                if write_csv == 1:
                    csv.write("%s,%s,%s\n" % (freq, db, refdist_out))
            else:
                csv.write("%s,%s,%s\n" % (freq, db, refdist_out))

# Add the SPREADTYPE field for a known landcover type (currently only NLCD)
def Add_SPREADTYPE(landcover, landcover_type):

    # Check if SPREADTYPE Field exists, if so, just use that
    exists = check_for_field(landcover, "SPREADTYPE")
    if exists == 1:
        arcpy.AddWarning("SPREADTYPE field already exists, and the existing field will be used.")
    else:
  
        if landcover_type == "autodetect":
            raise ValueError("This was a great idea, but it has not been scripted yet!") #**#
    
            # Extract values from raster
            raster_values = "etwas"
            # compare values against those used in known raster systems
            nlcd_value_lst = [11, 12, 21, 22, 23, 24, 31, 41, 42, 43,51, 52, 71, 72, 73, 74, 81, 82, 90, 95]        
            
            #**# Need a way of matching subsets
            if raster_values == nlcd_value_lst:
                landcover_type = "nlcd" #**# Recode landcover type
            #elif raster_value == othertype:
            else:
                raise ValueError("Could not autodetect land cover values for this land cover input")
            
            
        if str.lower(landcover_type) == "nlcd":
            # Wetlands are now classified as water
            codeblock = '''def calcNLCD(x):
                newval = ""            
                if x == 42 or x == 43:
                    newval = "CON"
                if x == 41:
                    newval = "HWD"
                if x == 51 or x == 52:
                    newval = "SHB"
                if x == 71 or x == 72 or x == 73 or x == 74 or x == 81 or x == 82:
                    newval = "HEB"
                if x == 31:
                    newval = "BAR"
                if x == 11 or x >= 90:
                    newval = "WAT"
                if x == 21 or x == 22 or x == 23 or x == 24:
                    newval = "URB"
                return(newval)'''
            
            arcpy.AddField_management(landcover, "SPREADTYPE", "TEXT")
            arcpy.CalculateField_management(landcover, "SPREADTYPE", "calcNLCD(!VALUE!)", "PYTHON", codeblock) # calculate it to be equal to value, then reclassify the value field in the next step

        if str.lower(landcover_type) == "landfire":
            arcpy.AddWarning("CLASSIFICATION UNDER DEVELOPMENT, only codes 3001 through 3968 are currently supported by this tool, and please check your classified raster for mis-classifications!")
            codeblock = '''def calcLANDFIRE(x):
                newval = ""            
                if x == 3001 or x == 3002 or x == 3006 or x == 3153 or x == 3218 or x == 3219 or x == 3220:
                    newval = "BAR"
                if x == 3221 or x == 3222 or x == 3294 or x == 3295:
                    newval = "BAR"
                if x == 3011 or x == 3016 or x == 3017 or x == 3019 or x == 3027 or x == 3028 or x == 3030:
                    newval = "CON"
                if x == 3031 or x == 3032 or x == 3033 or x == 3034 or x == 3043 or x == 3044 or x == 3049:
                    newval = "CON"
                if x == 3050 or x == 3051 or x == 3052 or x == 3054 or x == 3055 or x == 3057 or x == 3058:
                    newval = "CON"
                if x == 3059 or x == 3061 or x == 3062 or x == 3159 or x == 3208 or x == 3229 or x == 3231:
                    newval = "CON"
                if x == 3901 or x == 3902 or x == 3911 or x == 3921 or x == 3941:
                    newval = "CON"
                if x == 3029 or x == 3154 or x == 3160 or x == 3180 or x == 3900 or x == 3910 or x == 3920:
                    newval = "HWD"
                if x == 3940:
                    newval = "HWD"
                if x >= 3064 and x <= 3127:
                    newval = "SHB"
                if x >= 3210 and x <=3217:
                    newval = "SHB"
                if x == 3153 or x == 3220 or x == 3251 or x == 3252 or x == 3255 or x == 3259 or x == 3904:
                    newval = "SHB"
                if x == 3914 or x == 3923 or x == 3943 or x == 3960: # Orchard is treated as shrub? Or should it be hardwood?
                    newval = "SHB"
                if x >= 3135 and x <= 3146:
                    newval = "HEB"
                if x == 3164 or x == 3903 or x == 3913 or x == 3924 or x == 3944 or x == 3961 or x == 3963: # Vineyard as HEB?
                    newval = "HEB"
                if x == 3964 or x == 3965 or x == 3966 or x == 3967 or x == 3968:
                    newval = "HEB"
                if x >= 3181 and x <= 3182:
                    newval = "HEB"
                if x == 3292 or x == 3293: #**# 3293 is snow/ice, which is misclassified!
                    newval = "WAT"
                if x == 3296 or x == 3297 or x == 3298 or x == 3299:
                    newval = "URB"
                return(newval)'''

            arcpy.AddField_management(landcover, "SPREADTYPE", "TEXT")
            arcpy.CalculateField_management(landcover, "SPREADTYPE", "calcLANDFIRE(!VALUE!)", "PYTHON", codeblock) # calculate it to be equal to value, then reclassify the value field in the next step
            

'''
landcover = "C:/smt/csu/TWS_Case_Studies/Travel_Management_Plan/nlcd_p_cpy"
landcover_type = "nlcd"
Add_SPREADTYPE(landcover, landcover_type)
'''

# Function to extract landscover bounding rectangle
def get_extent(Extent_File):
    #import arcpy

    all_info = arcpy.Describe(Extent_File)
    extent_info = all_info.Extent
    # Clip tool needs: X-Minimum, Y-Minimum, X-Maximum, Y-Maximum
    xmin = extent_info.XMin
    ymin = extent_info.YMin
    xmax = extent_info.XMax
    ymax = extent_info.YMax
    
    out_extent = "%s %s %s %s" % (xmin, ymin, xmax, ymax)
    
    return(out_extent)

# Add time information to points file (create copy of points file? Probably the best idea. But then you could update this directly)
#**# This code intended to try to do a more explicit consideration of time, but has
# not been used and still needs improvement
def add_time(points_file, order_id, route_file, network, vehicle_speed, speed_units, analysis_dir, start_time):
    arcpy.env.overwriteOutput = True #Turn Overwrite on
    arcpy.CheckOutExtension("Network") # Check out network analyst

    '''
    points_file = "C:/smt/csu/TWS_Case_Studies/Travel_Management_Plan/temp_points.shp"
    route_file = "C:/smt/csu/TWS_Case_Studies/Travel_Management_Plan/temp_route_clip.shp"
    analysis_dir = "C:/smt/csu/TWS_Case_Studies/Travel_Management_Plan/temp/"
    network = analysis_dir + "split_route_ND.nd"
    start_time = 3600 * 9 # 9 am. 3600 seconds in an hour * 9 hours #**# For now, assume it's all the same day
    vehicle_speed = 20
    speed_units = "mph"
    order_id = "FID"
    '''    
    # Set up new routes
    new_points_file = analysis_dir + "temp_points2.shp" #**# Fix naming here
    new_routes = analysis_dir + "temp_routes2.shp"
    
    # Copy points file & add time field
    arcpy.CopyFeatures_management(points_file, new_points_file)
    arcpy.AddField_management (new_points_file, "TIME_S", "LONG") 
    arcpy.AddField_management (new_points_file, "HOUR", "SHORT")
    arcpy.AddField_management (new_points_file, "MIN", "SHORT")
    arcpy.AddField_management (new_points_file, "SEC", "SHORT")

    # Order of route travel is determined by order_id field

    #**# I think this is completely irrelevant now
    # Split routes at sample points
    arcpy.SplitLineAtPoint_management (route_file, new_points_file, new_routes, 1) #1 is search radius

    # Convert points to individual shape files #**# consider directory for this! This will be a lot of clutter
    point_lst, point_fill = points_to_shps(new_points_file, order_id, analysis_dir)

    # Calculate time for each point
    # Use time since start of day in s for analysis, but also add time fields that are human readable for checks.
    network_analyst_helper(new_points_file, point_lst, order_id, network, analysis_dir, vehicle_speed, speed_units, start_time)
       

# Use network analyst tools to calculate time schedule for points
def network_analyst_helper(new_points_file, point_lst, order_id, network, analysis_dir, vehicle_speed, speed_units, start_time):

    #Set up parameters for clinics and addresses:
    dest_id = order_id #unique ID field.
    inc_id = order_id #unique ID field
    searchdist = 5 #Maximum distance allowed for connecting a point to the network.
        
    #Set up network parameters
    #Set path to network information & set up network parameters
    #arcpy.env.workspace = "S:\\GIS_Center_Projects\\Simeonova\\Denmark_Health\\map_addresses\\map_addresses\\KDV_Data\\ALR_DATA.gdb"
    odlayer = "MyTest" #The facility layer to create
    distance = "Length" #Set the impedance layer to be based on distance
    #travel = "TRAVEL_TO"
    cutoff = 20000 #Cut point, in meters
    no_to_find = 1 #number of facilities to find
    dist_ref = ["Length"] #Not sure if this is needed, or not.  I think this outputs the accumulated distances
    uturn = "ALLOW_UTURNS"
    restrictions = []
    hierarchy = "NO_HIERARCHY"
    pathshape = "STRAIGHT_LINES" #This doesn't set it to straight line distance - this only has to do with what ArcGIS displays as a output shapefile (the values are network distance).  But if you calculate the length of these, you will get the straight line distances

    #Run network analysis

    # Loop through points
    n_reps = len(point_lst) - 1 # Looking at segments between points, so one fewer segments than points, assuming a straight line. Huh, what if you don't assume a straight line? That breaks this code plan in a bad, bad way.
    #**# Well, code for now to get network analysis working, then think about how to deal with route junctions
    #**# A vehicle can't split in two. Route *MUST* be continuous as an input. For loops, the user needs to break them into pieces.
    val_vec = [] # vector to hold values
    for i in xrange(n_reps): #Remove [0:1] when running for all destinations
        origin = analysis_dir + "point_%s.shp" % point_lst[i]
        destination = analysis_dir + "point_%s.shp" % point_lst[i + 1]
        
        #Calculate network distances    
        #Set up the cost distance framework
        arcpy.MakeODCostMatrixLayer_na (network, odlayer, distance, cutoff, no_to_find, dist_ref, uturn, restrictions, hierarchy, "", pathshape)

        #Add destination layer
        dlyr = destinations(odlayer,destination, dest_id,searchdist)

        #Add address (incident) layer
        ilyr = incidents(odlayer, origin, inc_id, no_to_find, cutoff, searchdist)

        #Solve analysis
        arcpy.Solve_na (odlayer, "SKIP", "TERMINATE")

        #Output network information
        this_val = netout(odlayer)
        val_vec.append(this_val)

    # convert speed to m/s
    meters_per_second = convert_speed(vehicle_speed, speed_units)

    # Convert lengths (m) to times (s)
    time_vec = [start_time]
    cum_val = start_time
    for val in val_vec:
        time = val * meters_per_second**-1 + cum_val
        time_vec.append(time)
        cum_val = time # update cumulative time
    
    # update new_points_file #**# Watch out for things being out of order. Is there a way to sort? No, but they should be in order when you create point_lst - control order there!
    id_check = []
    with arcpy.da.UpdateCursor(new_points_file, [order_id, "TIME_S", "HOUR", "MIN", "SEC"]) as rows:
        count = 0
        for row in rows:
            id_check.append(str(row[0]))
            this_time = time_vec[count]
            hour = this_time // 3600 # Floor division, get base number of hours 
            remainder = this_time % 3600
            minutes = remainder // 60 # Floor division get number of minutes
            sec = remainder % 60 # Remainder gives seconds
            row[1] = this_time
            row[2] = hour
            row[3] = minutes
            row[4] = sec
            rows.updateRow(row)            
            count += 1

    if id_check != point_lst:
        raise ValueError("Rows were not updated in the correct order, your updated",
                         "file is corrupted")

#Set up destination layer.  Requires the odlayer created by MakeODCostMatrixLayer_na.
def destinations(odlayer,destinations,uniqueid,search_dist = 500):
    '''Network layer, destinations, tolerance for locating destinations on the network (optional)'''
    #Search Distance = distance to search to assign an incident location
    #For incidents map the SourceID, SourceOID, PosAlong and SideOfEdge properties
    sub_layerf = "Destinations"
    field_mappingsf = "Name %s #; CurbApproach # Either side of vehicle; Attr_Minutes # 0; Attr_Minutes_Scales # 0; Attr_Length # 0" % uniqueid
    arcpy.AddLocations_na (odlayer, sub_layerf, destinations, field_mappingsf, search_dist)
    return sub_layerf

#Set up starting layer.  Requires the odlayer created by MakeODCostMatrixLayer_na
def incidents(odlayer,incidents,uniqueid,no_to_find,cutoff = 20000, search_dist = 500):
    '''Network layer, incidents, search tolerance for incidents (optional)'''
    #For incidents map the SourceID, SourceOID, PosAlong and SideOfEdge properties
    sub_layer = "Origins"
    search_dist = 500 # Distance to search to assign an incident location
    field_mappings = "Name %s #;TargetFacilityCount # %s; CurbApproach # Either side of vehicle; Attr_Minutes # 0; Attr_Minutes_Scales # 0; Attr_Length # 0;Cutoff_Length # %s" % (uniqueid,no_to_find,cutoff)
    arcpy.AddLocations_na (odlayer, sub_layer, incidents, field_mappings, search_dist)
    return sub_layer

def netout(odlayer):

    rlayer = "%s/Lines" % odlayer
    srows = arcpy.SearchCursor(rlayer)

    # Extract value from network analysis    
    for row in srows:
        value = row.getValue("Total_Length")

    return value


# Convert to meters per second
def convert_speed(vehicle_speed, speed_units):

    meters_per_hour = "NA"

    if str.lower(speed_units) == "m/s":
        return(vehicle_speed)

    if str.lower(speed_units) == "kph" or str.lower(speed_units) == "k/h":
        meters_per_hour = vehicle_speed * 1000
    
    if str.lower(speed_units) == "mph" or str.lower(speed_units) == "m/h":
        meters_per_hour = vehicle_speed * 1609.34 # Google says 1609.34 meters in 1 mile
    
    if meters_per_hour == "NA":
        raise ValueError("%s is not a recognized unit. Use kph for kilometers per hour, mph for miles per hour, or m/s for meters per second" % speed_units)
    meters_per_second = meters_per_hour * 3600**-1 # divide by number of seconds in an hour to get meters / sec

    return(meters_per_second)

# Split a multipoint shapefile into individual points
def points_to_shps(point_source_file, source_id_field, output_dir):
    point_lst = get_points(point_source_file, source_id_field, "all")
    point_lst = sorted(point_lst)
    point_fill = int(len(point_lst) / 10 )     
    point_fill = min(point_fill, 4)

    # Sort point_lst to ensure that points are run in order
    for point_id in point_lst:
        if len(point_id) > 4 or point_fill > 4:
            raise ValueError("This model uses the ArcGIS GRD Format. It has a maximum number of characters of 13 in file names. This leaves a maximum of 4 characters for the point ID. Please recode your point ID to be 4 or fewer characters Note that this should allow for >1.5 million unique point IDs using numbers and letters.")
        this_point = output_dir + "point_%s.shp" % point_id.zfill(point_fill)
        arcpy.Select_analysis(point_source_file, this_point, source_id_field + " = " + point_id)
        
    return(point_lst, point_fill)


## extract frequency information from a NMSim source file (.src)
# Description of NMSim Source file from Bruce Ikelheimer
    # The first line is a description.
    # The next three lines describe the icon that is used in NMSim.
    # The fifth line tells the number of frequencies used
    # The sixth line lists the number of source data files
    # The following lines provide the root name for the source data files.
    # After the source data files, the frequencies used in those files is listed.
    # The final line is the minimum height for a given source in meters.
def get_nmsim_source(in_src):

    # Define variables used later in script at arbitrary values that won't crash code
    lower_bound = 6
    upper_bound = 6

    freq_lst = []    
    
    line_number = 1 # DOES NOT FOLLOW the typical Python convention of starting with 0 instead of 1
    with open(in_src, "r") as nmsim_source:

        for line in iter(nmsim_source):
            
            # For troubleshooting purposes
            #print(line[:-1])  #Need to strip the carriage return otherwise there is an extra blank line after each line on the print output
            
            # If line number is other than specified, only increment the line number!
            if line_number == 5: # 5th line
                n_freq = int(line[:-1]) # -1 is to strip the carriage return
            if line_number == 6:
                n_sources = int(line[:-1])
                lower_bound = 6 + n_sources # Skip
                upper_bound = 6 + n_sources + n_freq
            if line_number > lower_bound and line_number <= upper_bound: 
                # if on a line containing a frequency, append it to freq_lst
                freq_lst.append(float(line[:-1]))
            if line_number == upper_bound + 1:
                minimum_source_height = float(line[:-1]) # get the minimum specified receiver height for a source
        
            line_number += 1
        
    # check that freq_lst was properly assigned
    if len(freq_lst) != n_freq:
        raise ValueError("Number of frequencies extracted from source file (%s) does not match number of frequencies specified in source file (%s)" % (len(freq_lst), n_freq))
    
    sound_level_lst = measurement_distance_lst = ["NA"] * n_freq
    
    # Convert minimum source height from feet to m
    minimum_source_height = minimum_source_height / 3.28084    
    
    return(freq_lst, sound_level_lst, measurement_distance_lst, minimum_source_height)

# Read in frequencies from a table
def read_table(in_freq, file_ending, freq_only):
    my_table = open(in_freq, "r")
    
    sep = "NA"    
    
    if file_ending == ".csv":
        sep = ","
    if file_ending == ".txt":
        sep = "\t"
        
    if sep == "NA":
        raise ValueError("We're sorry, the specified file ending %s"
                         "is not yet supported.\nPlease convert to .csv or"
                         ".txt.")
    # Define lists to contain outputs
    freq_lst= []
    sound_level_lst = []
    measurement_distance_lst = []
              
    count = 0
    # Go through each line of the file and extract the data            
    for line in my_table.readlines():
        line = line[:-1] # Strip carriage return
        fields = line.split(sep)

        # Check header for accuracy
        if count == 0:
            value_weighting = read_table_check_header(fields, freq_only)
            count += 1
            continue # Skip rest of processing for header
                    
        freq = fields[0]
        freq_lst.append(freq)

        # If all fields are required, update other lists as well        
        if freq_only == 0:
            dB = fields[1]
            measurement_distance = fields[2]
            sound_level_lst.append(dB)
            measurement_distance_lst.append(measurement_distance)
        
        count += 1        
      
    # "" is because reading the nmsim .src file returns an argument in that list place
    # and this is meant to maintain consistency in downstream code.
    return(freq_lst, sound_level_lst, measurement_distance_lst,"", value_weighting)

def read_table_check_header(fields, freq_only):
    value_weighting = "NA"    
    if str.lower(fields[0]) != "frequency":
        raise ValueError("The first column must be labeled 'Frequency'"
                         "\n In the given table, the first column was"
                         "%s" % fields[0])

    # If all fields are required
    if freq_only == 0:
             
        if str.lower(fields[1]) != "db" and str.lower(fields[1]) != "dba":
            raise ValueError("The second column must be either dB or dBA"
                             "In the given table, the second column was"
                             "%s" % fields[1])
                             
        if str.lower(fields[2]) != "measurement_distance":
            raise ValueError("The third column must be 'measurement_distance'"
                             "In the given table, the third column was"
                             "%s" % fields[2])
                             
        if str.lower(fields[1]) == "db":
            value_weighting = "Z"
        if str.lower(fields[1]) == "dba":
            value_weighting = "A"
    return value_weighting

def read_dbf(in_freq):
    raise ValueError("This function has not yet been scripted")
    # Read in a .dbf file (using a search cursor)


# Helper function to set up directories
def dir_setup(base_dir, results_dir, model, summarize):

    # Check if base_dir exists, if not throw an error
    if not os.path.exists(base_dir):
        raise ValueError("Base directory %s does not exist." % base_dir + "Please make sure you have specified the correct "
        "directory and that it exists. Then try re-running the model. If you are using Python, you can create the "
        "directory with os.makedirs('%s')" % base_dir)

    # Set up directories
    model_dir = base_dir + "%s/" % model    # Set up main model directory
    idir = model_dir + "temp/"          # Create a sub-folder for intermediates
    pdir = model_dir + "points/"          # Create a sub-folder for point-specific results
    freq_dir = model_dir + "frequency_propagation/" # Create a folder for frequency-specific results
    # Note: there is also a point_dir defined later to contain point-specific results

    dir_lst = [results_dir, idir, pdir, freq_dir]

    if summarize < 3:
        # Delete any pre-existing model directory (but only if doing a new model run!)
        if os.path.exists(model_dir):    
            shutil.rmtree(model_dir)

    # Make required directories (makedirs will recursively make model_dir)
    for my_dir in dir_lst:
        # This line is here, because sometimes directories don't delete properly
        make(my_dir)
        
    return(model_dir, idir, pdir, freq_dir)

# Helper function to get list of points for running multiple points
def get_points(point_source_file, source_id_field, n_points):
    #import arcpy

    # Check if source_id_field is in shapefile, if not, throw a meaningful error
    fields_lst = arcpy.ListFields(point_source_file)
    field_ok = 0
    for test_field in fields_lst:
        if test_field.name == source_id_field:
            field_ok = 1

    if field_ok == 0:
        raise ValueError("source id, %s, is not a field in the source file, %s" % (source_id_field, point_source_file))


    with arcpy.da.SearchCursor(point_source_file, source_id_field) as srows:    

        point_lst = []
    
        count = 0
        for row in srows:
            
            # Allow an early break if only a few points should be run.        
            if count == n_points:
                break
            count += 1
    
            point = row[0]
            #point = 1        
            point_id = str(point)
            point_lst.append(point_id)
    
    if count == 0:
        raise ValueError("There was a problem extracting point information from your source file. Does the source file exist? Does it contain data? Please fix your source file and try again.")
    return(point_lst)  

# Function write values to table
def extract_values(model, sources, source_id_field, receivers, in_freq, results_table, results_path, n_points):
    #import arcpy, arcpy.sa, os
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True


    # Loop through sources    
    srows = arcpy.SearchCursor(sources)

    count = 0
    for row in srows:

        # Add option to run for a subset of points        
        if count == n_points:
            break
        count += 1
        
        point = row.getValue(source_id_field)

        freq_lst = extract_frequencies(in_freq)[0]
        
        for freq_s in freq_lst:
        
            raster = results_path + "pt%s_pr%s.tif" % (point, freq_s)
            temp_table = arcpy.sa.Sample([raster], receivers, results_path + "temp_table.dbf", "NEAREST", "ID")
    
            # Extract results from temp_table into a Python-readable format
            results_dict = {}
            cur2 = arcpy.SearchCursor(temp_table)
            row2 = cur2.next()
            while row2 <> None:
                raster_name = str.split(raster, "/")[-1]
                field = raster_name[0:10]
                rc = row2.getValue("receivers") # Get number of receiver
                rc_value = row2.getValue(field) # Get value at that receiver
                results_dict[rc] = rc_value
                
                row2 = cur2.next()
            
            del cur2, row2
        
            # Write results from this source to file (#should I extract to dictionary and then do all writing at once?)
            # Desired output table format: Model Frequency Source Receiver Value
            existing_file = 1
            if not os.path.exists(results_table):
                existing_file = 0
        
            rt = open(results_table, 'a')
    
            # If file is being created, print header row (indicater is needed, because os.path.exists after the open call will always return TRUE)
            if existing_file == 0:        
                print >> rt, "Model,Frequency,Source,Receiver,Value"
                
            for key in sorted(results_dict.iterkeys()):
                print >> rt, "%s,%s,%s,%s,%s" % (model, freq_s, point, key, results_dict[key])            
        
    # Delete search cursor objects to remove any associated schema locks
    del srows, row

# Set up for timing results
def setup_timing(timelog):
    my_times = "NA" # Create a string for when timing is not occuring
    my_time_labels = "NA"
    if timelog != "none":
        import time
        this_time = time.time()
        my_times = [this_time]
        my_time_labels = ["start"] # Create a list of labels for the time for the log

    return(my_times, my_time_labels)
    
#**# Duplicated in spreadgishlpr to avoid having to re-import this function
# from that script
# Create a timer function for timing code
def get_time(timelog, my_times, my_time_labels, label):
    import time    
    if timelog != "none":
        this_time = time.time()
        my_times.append(this_time)
        my_time_labels.append(label)

    return(my_times, my_time_labels)

# Function to output final timing
def output_timing(timelog, my_times, my_time_labels, delete_existing_timelog):
    '''timelog is the path and file for the timing results, my_times give the
    time since epoch in seconds, my_time_labels label each time step, delete_existing
    indicates whether or not to delete the existing time file.'''
    import os
        
    if timelog != "none":
        # Check that times and labels are correct
        if len(my_times) != len(my_time_labels):
            raise ValueError("There must be one label for every timed event")
        
        # Delete any pre-existing timelog
        if os.path.exists(timelog) and delete_existing_timelog == 1:
            os.remove(timelog)

        # Calculate start & end times
        start_time = my_times[0] # Just the first timestep. This will always be substracted
        end_time = my_times[-1] - start_time  # Need the subtraction to convert from time since last epoch.
    
        t = open(timelog, 'a')
        
        # Create header row
        print >> t, "Task, Process Time, Elapsed Time, Percent of Total"

        
        #Start at 1, no sense calculating elapsed for the start time!
        for i in range(1,len(my_times)):
            # Get time since code started
            elapsed_time = my_times[i] - start_time
            # Get time for this specific item
            item_time = my_times[i] - my_times[i - 1]
            perc_time = (item_time * end_time**-1) * 100 # Division is screwy in Python, so instead multiply by the reciprocal
            print >> t, "%s,%s,%s,%.1f" % (my_time_labels[i], item_time, elapsed_time, perc_time)
        
        # Stop writing to file upon completion
        t.close()

# Function to give a more informative error message when the source point
# is not within the specified extent (and also check that the source_id_field field
# is actually unique and only one point is selected)
def check_selection(Sound_Source):
    #**# These two lines can be deleted once the updated cursor code is tested and works
    #cur = arcpy.SearchCursor(Sound_Source)
    #row = cur.next()
    
    with arcpy.da.SearchCursor(Sound_Source, "FID") as srows:   
        #row = cur.next()

        count = 0
        #while row:
        for row in srows:
            count += 1
            #row = cur.next()

    if count == 0:
        raise ValueError("No points were included in the selection table. Please check that the point is within the model extent")
    
    if count > 1:
        raise ValueError("source_id_field field was not unique. You must select a unique ID field")

# Make the directory recursively, if it does not already exist. I hate the 3 lines, and want a 1-line function!
def make(in_path):
    if not os.path.exists(in_path):
        os.makedirs(in_path)

# Make a NMSim-style source file
def make_srcfile(csv_lst, file_type, analysis_type, outpath, source_name, source_description, avg_name, source_type, AGL_m, ref_dist_m):

    ''' For testing:
    csv_lst = drilling_src
    file_type = drill_file_type
    analysis_type = drill_analysis_type
    source_name = drill_source_name
    source_description = drill_source_desc
    avg_name = drill_avg_name
    AGL_m = drill_AGL_m
    ref_dist_m = drill_ref_dist_m
    '''

    print("WARNING: INPUTS FOR THIS FUNCTION MAY CHANGE TO ACCOMODATE SOURCE SPEED AND/OR SOURCE DIRECTIONALITY")

    ft_per_m = 3.28084

    speed_vec = [100] # Does it matter? If there is only one speed, then it shouldn't matter. Think about this.    
    ref_dist_ft_vec = [ref_dist_m * ft_per_m]
    src_values_lst = ["NA"] * len(speed_vec)
    #**# Think about best time to correct for A-weighting - probably want to do that after reading in files but before QC. Alternatively, can do so with Damon's program after QC.    

    # check number of characters in source_name, avg_name
    if len(source_name) > 12:
        message = "source_name must be less than or equal to 12 characters (input was %s characters)" % len(source_name)
        raise ValueError(message)
    if len(avg_name) > 6: 
        message = "avg_name must be less than or equal to 6 characters (input was %s characters)" % len(avg_name)
        raise ValueError(message)
        

    for k in xrange(len(csv_lst)):
        in_csv = csv_lst[k]

        # Begin to read in .csv & identify whether it is a Larson Davis-type file, an NVSPL file, or a simple entry (should I autodetect this??? Or have the user specify it?)
        all_values, freq_lst = read_raw_src(in_csv, file_type)
        
        # Check all_values for obvious outliers
        #all_values = check_raw_values(all_values, STUFF)
        
        # Decide on handling of file
        if analysis_type == "steady":
            # Compare median value to mean
            src_values = compare_mean_median(all_values, freq_lst)
        else:
            src_values = all_values
        #**# Think about adding a source_type for a cyclical source (or just take a mean value??)
    
        #**# Think about adding support for a vehicle type source (approach, peak, & decline, only use time around peak)

        src_values_lst[k] = src_values

    # if no error raised above:
    # Increasing the number of speed variants could be very useful for cyclic sources - you could specify a particular part of the cycle.
    remove_a_weighting = 0 # Hard-coded to not remove A-weighting - this can be done above before the file checks
    # Convert m to ft
    #AGL_ft = AGL_m * ft_per_m

    # Extract lower & upper frequency from freq_lst
    #**# Reenable this code if one starts getting errors    
    #for m in freq_lst:
    #    freq_lst[m] = float(freq_lst[m])
    lower_frequency = min(freq_lst)
    upper_frequency = max(freq_lst)

    srcfile = write_nmsim_src(outpath, src_values_lst, source_name, source_description,
                    avg_name, source_type, freq_lst, lower_frequency,
                    upper_frequency, remove_a_weighting, AGL_m,
                    ref_dist_ft_vec, speed_vec)
    
    return srcfile    
    
# Read in raw src values
def read_raw_src(in_csv, file_type):

    if file_type == "Larson Davis" or file_type == "LD":
        src_values, freq_lst = extract_larson_davis(in_csv)

    if file_type == "NVSPL":
        raise ValueError("NVSPL format is not yet supported. Sorry for the inconvenience. You may manually compute the final summary, and choose 'simple' to extract the values")
        
    if file_type == "simple":
        output = extract_frequencies(in_csv)
        freq_lst = output[0]
        src_values = output[1]
        #**# YOU COULD ALSO GET MEASUREMENT DISTANCE THIS WAY.

    return(src_values, freq_lst)

# Extract data in converted larson davis csv format
# Data can be converted using the slm utility
# http://www.larsondavis.com/support/softwareproductssupport/SoftwareUtilities.aspx
#**# NOTE: THERE IS A G4 SLM utility, perhaps we should use that instead
def extract_larson_davis(in_csv):    

    # Read in csv file
    with open(in_csv, 'r') as csv:
        src_values = {} # initialize a dictionary to hold the data
        one_third_octave_bands = [8, 10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000]
        one_third_octave_indices = []
        old_band = 0

        #**# Damon suggested there was a readlines function in R that one could use grep to find the correct line. Then can skip a couple of lines and read the whole thing in as a table. Might be an option to speed this up.            
        begin = 0
        count = -1
        for line in csv.readlines():
            count += 1
            # remove carriage return & split line into fields
            line = line[:-1]
            fields = line.split(",")
            
            if fields[0] == "Time History":
                begin = 1
                continue # Don't do any more processing

            # Set up header row
            if begin == 1:
                begin = 2
                
                # Loop through fields & extract data
                # For now, making the simplifying assumption that Emma's data are in a standard Larson Davis 831 output format
                for i in xrange(len(fields)):
                    field = fields[i]
                    if field == "Record #":
                        rec_index = i
                    if field == "Date":
                        date_index = i
                    if field == "Time":
                        time_index = i
                    one_third_octave_indices, old_band = extract_octave_indices(field, one_third_octave_bands, one_third_octave_indices, i, old_band)
                continue # Do not continue with rest of code with header data

            # First line does not contain data, so don't do anything with it. Update the begin index
            if begin == 2:
                begin = 3
                continue

            # Stop when you reach the bottom of the data table
            if fields[4] == "Key": # There is a "Key" above this point but the begin == 2: continue code bypasses. Careful if moving either this or the above block of code
                break
            
            if begin == 3:
                # Begin extracting data into a dictionary (#**# is there a better data format?)
                src_values = extract_src_data(src_values, fields, rec_index, date_index, time_index, one_third_octave_indices)


    # check that all 1/3 octave bands were extracted
    if len(one_third_octave_bands) != len(one_third_octave_indices):
        raise ValueError("not all 1/3 octave band values were found in the data set")

    # or I could just rename one_third_octave_bands as freq_lst & update all the above code
    freq_lst = one_third_octave_bands

    return(src_values, freq_lst)
                    
# function to extract field indices for 1/3 octave bands
def extract_octave_indices(field, one_third_octave_bands, one_third_octave_indices, i, old_band):
    for j in xrange(len(one_third_octave_bands)):
        band = one_third_octave_bands[j]
        # Fields that have text will throw an error with the float command
        # Try to float field, otherwise skip it.
        try:
            field = float(field)            
        except:
            continue

        #band > old_band checks that the bands are increasing order
        if band == field and band > old_band:
            one_third_octave_indices.append(i)

            # update the old_band designation
            old_band = band

    return(one_third_octave_indices, old_band)                                

def extract_src_data(src_values, fields, rec_index, date_index, time_index, one_third_octave_indices):
    record = fields[rec_index]
    date = fields[date_index].split("/")
    month = date[0]
    day = date[1]
    year = date[2]
    time = fields[time_index] #.split(":")
    #hour = time[0]
    #minute = time[1]
    #second = time[2]
    
    key = "%s %s%s%s %s" % (record, year, month.zfill(2), day.zfill(2), time)

    # get the values for each 1/3 octave band    
    one_third_octave_values = []
    for index in one_third_octave_indices:
        value = float(fields[index])
        one_third_octave_values.append(value)
    
    # put the 1/3 octave band values into a dictionary keyed by record, date, and time
    src_values[key] = one_third_octave_values
    
    # Return the dictionary for the next loop
    return src_values

def check_raw_values(all_values, STUFF):
    raise ValueError("Not yet scripted")
    # Check whether there are obvious problems with the data set

def compare_mean_median(all_values, freq_lst, median_tolerance = 5, Q10_tolerance = 10, Q90_tolerance = 10):
    
    #**# Discuss w/ Kurt & Damon what reasonable tolerance values are.    
    print "Please discuss soundprophlpr compare_mean_median with Kurt & Damon"    

    # NEED TO LOOP THIS BY FREQUENCY    
    src_values = ["NA"] * len(freq_lst)
    raise_value_error = 0
    for i in xrange(len(freq_lst)):
        
        freq = freq_lst[i]

        # Calculate mean, median, L10, L90 from all_values data table
        mean_value, median_value, Q10_value, Q90_value = calculate_src_statistics(all_values, i)    
        # Set the value for this frequency to the mean value
        src_values[i] = mean_value
    
        # Check if median, L10, & L90 are similar to the mean
        error_vec = [1,1,1] #Default: All three errors: median, Q10, Q90
        if mean_value <= (median_value + median_tolerance) and mean_value >= (median_value - median_tolerance):
            error_vec[0] = 0
            
        if mean_value <= (Q10_value + Q10_tolerance) and mean_value >= (Q10_value - Q10_tolerance):
            error_vec[1] = 0
            
        if mean_value <= (Q90_value + Q90_tolerance) and mean_value >= (Q90_value - Q90_tolerance):
            error_vec[2] = 0
    
        # If not, figure out what to do
        if max(error_vec) > 0:
            print "Frequency %s: median: %s, L10: %s, L90: %s." % (freq, error_vec[0], error_vec[1], error_vec[2])
            "Values of 1 indicate that the calculated value relative to the mean is not within the specified thresholds."
            "Thresholds are %s,%s,%s for median, L10, and L90 respectively" % (median_tolerance, Q10_tolerance, Q90_tolerance)
            "while mean, median, L10, L90 are %s, %s, %s, %s" % (mean_value, median_value, Q10_value, Q90_value)
            raise_value_error = 1

    if raise_value_error == 1:
        raise ValueError("Steady source is not constant enough for automated processing. Please change tolerance thresholds or manually assess the data & enter the desired values as entry type 'simple'")

    return(src_values)

#Calculate mean, median, 10% & 90% quantiles
def calculate_src_statistics(all_values, i):
    # Somehow calculate this from the dictionary

    #**# THINK ABOUT A BURN-IN TIME option & BURN-OUT - how many rows from the beginning and the end to exclude    

    mean_value = "NA"
    values_lst = ["NA"]
    N = 0
    
    for key in all_values:
        value = all_values[key][i] # i just gets the frequency of interest
        if mean_value == "NA":
            mean_value = value
        else:
            mean_value = sp_update_mean(mean_value, value, N)
        if values_lst[0] == "NA":
            values_lst[0] = value
        else:
            values_lst.append(value)             
        # update N after processing the record
        N += 1

    # Convert list of values for each function to a quantile
    median_value = compute_quantile(values_lst, 50)
    Q10_value = compute_quantile(values_lst, 10)
    Q90_value = compute_quantile(values_lst, 90)
    
    return [mean_value, median_value, Q10_value, Q90_value]

# Formula for updating a mean   #MODIFIED FROM PICEANCEHLPR.PY - IT SHOULD PROBABLY JUST REFERENCE THIS VERSION, FOR CONSISTENCY             
def sp_update_mean(mvl_value, value, N):
    # take the average of the old mean weighted by the previous sample size, plus the new value, then divide by the new sample size
    new_mean = (mvl_value * N + value) * (N + 1)**-1
    return(new_mean)

# Function to convert a list of values by frequency to a single quantile value per frequency
def compute_quantile(values_lst, quantile):
    from numpy import percentile

    quantile_value = percentile(values_lst, quantile)
   
    return(quantile_value)

# Create a nmsim source for an input data set
# A Python command line implementation of NMSIM Source File Editor was Public Domain, & was written by Damon Joyce at the National Park Service Natural Sounds, Night Skies Division
# NMSIM Source File Editor was Public Domain, & was written by Damon Joyce at the National Park Service Natural Sounds, Night Skies Division
def write_nmsim_src(outpath, src_values_lst, source_name, source_description,
                    avg_name, source_type, freq_lst, lower_frequency, upper_frequency,
                    remove_a_weighting, AGL_m, ref_dist_ft_vec,
                    speed_vec):

    # WRITE NMSIM SRC FILE
    # Keep outpath flexible, but have it try to use tbx_root & put into the Sources directory
    srcfile = write_src_file(outpath, source_name, source_type, source_description, freq_lst, avg_name, speed_vec, AGL_m)

    # WRITE AVG FILES
    # Loop through to allow variation with speed #**# Not actually functional, but adding the infrastructure now to make it easy to add
    for i in xrange(len(speed_vec)):
        speed = speed_vec[i]
        ref_dist_ft = ref_dist_ft_vec[i]
        src_values = src_values_lst[i]
        write_avg_file(outpath, avg_name, speed, freq_lst, src_values, ref_dist_ft)
                   
    return srcfile

def write_src_file(outpath, source_name, source_type, source_description, freq_lst, avg_name, speed_vec, AGL_m): 
    outfile = outpath + source_name + ".src"
    
    # Delete any existing src file with the same name
    #if os.path.exists(outfile):
    #    os.remove(outfile)

    with open(outfile, 'w') as of:
        # Line 1: source type & description
        of.write("%s-%s\n" % (source_type, source_description))
        # Line 2: # identifier for the icon that is displayed in NMSim
        thing = "808VCZNZRVR080" 
        of.write("s %s\n" % thing)
        # Line 3 & 4: The dimensions of the icon when displayed in NMSim
        of.write("xsize: 1.5\n")
        of.write("ysize: 1.5\n")
        # Line 5 - number of frequencies
        of.write("%s\n" % len(freq_lst))
        # Line 6: the number of avg files
        of.write("%s\n" % len(speed_vec))
        # Line 7+: a line for each avg file
        for speed in speed_vec:
            buffer_size = 9 - len(str(speed))
            this_buffer = " " * buffer_size # Add spaces based on the buffer size
            speed_plus_buffer = "%s%s." % (this_buffer, speed) # length of speed_plus_buffer needs to be 10 & needs to end with a period.
            of.write("%s%s.avg%s\n" % (avg_name, speed, speed_plus_buffer))
        # Next lines: one for each frequency band
        for freq in freq_lst:
            of.write("%s\n" % float(freq))
        # Last line: minimum source height (m)
        of.write("%s\n" % (float(AGL_m) * 3.28084)) # The .src file entry is actually in feet.

    return outfile

# Write avg_files for NMSim
def write_avg_file(outpath, avg_name, speed, freq_lst, src_values, ref_dist_ft):
    
    
    avg_file = "%s%s%s.avg" % (outpath, avg_name, speed)
    # Delete any existing avg file with the same name
    #if os.path.exists(avg_file):
    #    os.remove(avg_file)
    with open(avg_file, 'w') as af:
        #Line 1 = header
        header = "   Theta    PHI       "
        for freq in freq_lst:
            if freq == freq_lst[-1]:
                header = "%s%s" % (header, freq)
            else:
                freq_buffer_size = 5 - len(str(freq))
                freq_buffer = " " * freq_buffer_size
                header = "%s%s%s" % (header, freq, freq_buffer)
        af.write("%s\n" % header)
        
        #Line 2 -blank (but has a single space?)
        af.write(" \n")
        # Remaining Lines: one for each Theta value. For now keeping Phi constant at 0, and making sources omni-directional.
        #**# Consider modifying this to get full source functionality.
        theta_values = seq(5,180,5)
        theta_lst = [".0"]
        for value in theta_values:
            theta_lst.append(float(value))
        for theta_value in theta_lst:
            theta_buffer_size = 8 - len(str(theta_value))
            theta_buffer = " " * theta_buffer_size
            #**# Talk w/overflights people about what they need in this respect
            phi = "0.0" #**# Consider making this an input? Esp. for the aircraft overflights people. Also consider making it flexible
            phi_buffer_size = 8 - len(str(phi))
            phi_buffer = " " * phi_buffer_size
            first_buffer = " " * 4
            
            source_level_entries = ""
            count = 0
            for source_level in src_values:
                count += 1
                source_level_cB = float(source_level) * 10 # Convert to centibells, the format used by NMSim
                # Reference distance correction from Damon: _distCorr = 20 * Math.Log10(_distCorr / 1000.0);
                # NMSim uses a reference distance of 1000 ft. Need to correct the SPL values based on the reference distance
                spherical_spreading_loss_correction = 20 * math.log10(ref_dist_ft * 1000**-1)
                source_level_cB = source_level_cB + spherical_spreading_loss_correction * 10 # ssl correction must also be in centibels!
                source_level_cB = int(round(source_level_cB, 0))
                entry_buffer_size = 5 - len(str(source_level_cB))
                entry_buffer = " " * entry_buffer_size
                source_level_entries = "%s%s%s" % (source_level_entries, entry_buffer, source_level_cB)
                if count <5 and theta_value < 2:
                    print source_level_entries
                    print "%s%s%s" % ("b", entry_buffer, "B")
            af.write("%s%s%s%s%s%s\n" % (theta_buffer, theta_value, phi_buffer, phi, first_buffer, source_level_entries))
    
# Create a nmsim source for an input data set
def write_nmsim_src_NMSIM_SOURCE_FILE_EDITOR_VERSION(outpath, src_values_lst, nmsim_source_file_editor_path, source_name,
                    source_description, avg_name, source_type, freq_lst,
                    remove_a_weighting, AGL_m, left_column_vec, ref_dist_m_vec, speed_variants, speed):
    raise ValueError("NMSIM_SOURCE_FILE_EDITOR does not take command line arguments")

    import subprocess

    # Convert m to ft
    ft_per_m = 3.28084 # Conversion from Google, 2016-07-11
    AGL_ft = AGL_m * ft_per_m

    # Extract lower & upper frequency from freq_lst
    lower_frequency = min(float(freq_lst))
    upper_frequency = max(float(freq_lst))

    command_pt1 = "%s %s %s %s %s" % (nmsim_source_file_editor_path, source_name, source_description, avg_name, source_type)
    command_pt2 = "%s %s %s %s %s" % (lower_frequency, upper_frequency, AGL_ft, remove_a_weighting, outpath)
    command = "%s %s" % (command_pt1, command_pt2)
                        
    # Loop through to allow variation with speed #**# Not actually functional, but adding the infrastructure now to make it easy to add
    for i in xrange(speed_variants):
        ref_dist_m = ref_dist_m_vec[i]
        ref_dist_ft = ref_dist_m * ft_per_m
    
        left_column = left_column_vec[i]
        src_values = src_values_lst[i]
        src_command_part = ""
        for value in src_values:
            src_command_part = "%s %s" % (src_command_part, value)
        
        command = "%s %s %s" % (left_column, src_command_part, ref_dist_ft)

    # Run NMSIM Source File Editor
    #  NMSIM Source File Editor was Public Domain, & was written by Damon Joyce at the National Park Service Natural Sounds, Night Skies Division
    subprocess.call(command) # Modify to accept command line arguments & to allow people to choose the directory in which to save. Think about how to deal with multiple speed variants

# Check whether all data are projected, and if all data are in the same projection
def check_projections(file_lst):
    
    proj_error = 0    
    mismatch_error = 0
    e_messages = ""
    
    # 11/8/2016 Changed: Multiple names correspond to the same projection and codes correspond to multiple UTMs!
    # Now using the spatial reference properties, and checking that they are the same.
    proj_names = []    
    linear_units = []
    #angular_units = [] # Got blank, when it should have been degrees. skipping this for now
    false_easting = []
    false_northing = []
    central_meridian = []
    scale_factor = []
    latitude_of_origin = []
    #Datum = [] # Got blank or 0. So not sure what's going on here either. Skipping for now
    
    # Loop through files
    for my_file in file_lst:
        describe_object = arcpy.Describe(my_file)
        if describe_object.SpatialReference.type != "Projected":
            e_messages = e_messages + "Input %s is not projected.\n"
            proj_error = 1

        # Changed from projection code to name because codes are not unique between UTM zones!
        SR = describe_object.SpatialReference
        proj_names.append(SR.name)            
        linear_units.append(SR.linearunitname)
        #angular_units.append(SR.angularunitname)
        false_easting.append(SR.falseeasting)
        false_northing.append(SR.falsenorthing)
        central_meridian.append(SR.centralmeridian)
        scale_factor.append(SR.scalefactor)
        latitude_of_origin.append(SR.latitudeoforigin)
        #Datum.append(SR.datum)

    # Check if all projections are the same
    #n_projections = len(list(set(proj_codes)))  # Set == unique from R, except it is a special object type and needs to be converted back to a list
    #  Get list of different values
    names = list(set(proj_names))    
    units = list(set(linear_units)) 
    eastings = list(set(false_easting))
    northings = list(set(false_northing))
    meridians = list(set(central_meridian))
    factors = list(set(scale_factor))
    latitudes = list(set(latitude_of_origin))

    # Get number of different values
    n_units = len(units)
    n_eastings = len(eastings)
    n_northings = len(northings)
    n_meridians = len(meridians)
    n_factors = len(factors)
    n_latitudes = len(latitudes)
    # Names is deliberately excluded - NAD83_UTM_zone_11N is the same as NAD_1983_UTM_Zone_11N but the names won't match! I don't know why multiple duplicate names even exist!
    n_projections = max(n_units, n_eastings, n_northings, n_meridians, n_factors, n_latitudes)
    if n_projections != 1:
        e_messages = e_messages + "Not all files are in the same projection. %s different projections were used\n" % n_projections
        e_messages = e_messages + "Projection names were: %s\n" % names
        e_messages = e_messages + "Linear Units were: %s\n" % units
        e_messages = e_messages + "False eastings were: %s\n" % eastings
        e_messages = e_messages + "False northings were: %s\n" % northings
        e_messages = e_messages + "Central Meridians were: %s\n" % meridians
        e_messages = e_messages + "Scale factors were: %s\n" % factors
        e_messages = e_messages + "Latitudes of origin were: %s\n" % latitudes
        mismatch_error = 1

    if proj_error == 1 or mismatch_error == 1:      
        raise ValueError("%sAll files must be projected and in a common (UTM) coordinate system" % e_messages)

    # No return, just designed to throw an error message if there is a problem.

def setup_message_display(point_counter):
    display_messages = 0            
    if point_counter < 5:
        display_messages = 1
    if (point_counter % 5) == 0:
        display_messages = 1
    
    '''      
    # No point in any of these as long as the Initializing and Start Euclidian... messages can't be turned off!
    if point_counter < 30 and (point_counter % 5) == 0:
        display_messages = 1
    if point_counter < 100 and (point_counter % 10) == 0:
        display_messages = 1
    if point_counter < 1000 and (point_counter % 100) == 0:
        display_messages = 1
    if point_counter < 10000 and (point_counter % 500) == 0:
        display_messages = 1
    '''
    
    return(display_messages)

# Try to make an R-like sequence function (see xrange function)
def seq(start, stop, by):    

    #**# For now, this only works for ascending sequences    
    # Test if stop is greater than or less than start -     
    if stop < start:
        raise ValueError("'stop' must be greater than 'start'")
    # Check if by is positive
    if by <= 0:
        raise ValueError("'by' must be positive")
    
    value = start
    my_lst = []
    
    while value <= stop:
        my_lst.append(value)
        value = value + by
    
    return(my_lst)

# Speed of sound calculation from Bruce Ikelheimer
def calculate_speed_of_sound(temp_C):
    # Convert T to K
    temp_k = float(temp_C) + 273.15
    # Speed of sound (m/s)
    c = math.sqrt(1.4 * 286.9 * temp_k) 
    return c

# Add a function to add a raster estimating propagation time based on the speed of sound
def create_propagation_time_raster(source_point, temp, point_dir, model_dir, point_id, point_fill):
    # Set up a point label
    point_lbl = "pt%s" % point_id.zfill(point_fill)

    # Get speed of sound
    c = calculate_speed_of_sound(temp)
    # Convert T to K
    #temp_k = float(temp) + 273.15
    ## Speed of sound (m/s)
    #c = math.sqrt(1.4 * 286.9 * temp_k) 
    time_per_m = c**-1 # inverse is seconds per meter

    # Create a cost raster where the cost is the speed of sound
    speed_raster = arcpy.sa.CreateConstantRaster(time_per_m, "FLOAT")
    speed_raster.save(point_dir + "speed_of_sound.tif")
    time_raster = arcpy.sa.PathDistance(source_point, speed_raster)
    time_raster.save(model_dir + "propagation_time_%s.tif" % point_lbl)
    
# Convert a time-course to an animated gif for visualization purposes
# Convert a series of images into a movie
def results_to_movie(path):
    raise ValueError("This still is not yet scripted")

    # List contents of the path

    # Sort contents of path

    # Call a function to compile the images into a movie

# Read NMSim inputs from point_source_file
def read_SoundSource(in_freq, SoundSource, allow_defaults):
    
    # Get list of fields to know what is present and what is omitted
    fields = arcpy.ListFields(SoundSource)
    field_indices = [0,0,0,0,0,0]
    order_index = []
    my_fields = []
    # Go through fields and see which are present
    for field in fields:
        field = str(field.name) # Change field to the actual field name & convert to string
        if str.lower(field) == "head" or str.lower(field) == "heading" or str.lower(field) == "course":
            field_indices[0] = 1
            my_fields.append(field)
            order_index.append(0)
        if str.lower(field) == "roll":
            field_indices[1] = 1
            my_fields.append(field)
            order_index.append(1)
        if str.lower(field) == "pitch":
            field_indices[2] = 1
            my_fields.append(field)
            order_index.append(2)
        if str.lower(field) == "speed" or str.lower(field) == "vel" or str.lower(field) == "velocity":
            field_indices[3] = 1
            my_fields.append(field)
            order_index.append(3)
        if str.lower(field) == "engpow" or str.lower(field) == "engine_power" or str.lower(field) == "engine":
            field_indices[4] = 1
            my_fields.append(field)
            order_index.append(4)
        if str.lower(field) == "s_offset" or str.lower(field) == "source_offset" or str.lower(field) == "soffset" or str.lower(field) == "msl" or str.lower(field) == "msl_m":
            field_indices[5] = 1
            my_fields.append(field)
            order_index.append(5)

    # Check if any default values are required, and check whether default values
    # are allowed (allow y or true or t as alternatives for yes)
    if str.lower(allow_defaults) != "yes" and str.lower(allow_defaults) != "y" and str.lower(allow_defaults) != "t" and str.lower(allow_defaults) != "true":
        for i in field_indices: 
            if i == 0:
                m1 = "Please make sure the following fields are present in your shape file, "
                m2 = "or set allow_defaults to 'YES': "
                m3 = "HEAD, PITCH, SPEED, ENGPOW, S_OFFSET. "
                m4 = "See User's Manual or help for more information."
                raise ValueError(m1 + m2 + m3 + m4)

    # An empty fields object will crash the search cursor
    if len(my_fields) > 0:

        # Use a search cursor to read parameters from the shapefile
        with arcpy.da.SearchCursor(SoundSource, my_fields) as cursor:
            for ss in cursor:
        
                # Check that all fields have been assigned, if not, assign a default
                if field_indices[0] == 0:
                    head = 0
                else:
                    head_index = get_right_field(order_index, 0)
                    head = ss[head_index]
                if field_indices[1] == 0:
                    roll = 0
                else:
                    roll_index = get_right_field(order_index, 1)
                    roll = ss[roll_index]
                if field_indices[2] == 0:
                    pitch = 0
                else:
                    pitch_index = get_right_field(order_index, 2)
                    pitch = ss[pitch_index]
                if field_indices[3] == 0:
                    vel = 0               
                else:
                    vel_index = get_right_field(order_index, 3)
                    vel = ss[vel_index]
                if field_indices[4] == 0:
                    engpow = 0
                else:
                    engpow_index = get_right_field(order_index, 4)
                    engpow = ss[engpow_index]
                if field_indices[5] == 0:
                    source_offset = 1
                else:
                    so_index = get_right_field(order_index, 5)
                    source_offset = ss[so_index]
            
            source_info = [head, roll, pitch, vel, engpow, in_freq]

    if len(my_fields) == 0:
        source_info = [0, 0, 0, 0, 0, in_freq]
        source_offset = 1
        
    return(source_info, source_offset)
    
# Function to return the correct field when the fields come in in any order
def get_right_field(order_index, target_value):
    for i in xrange(len(order_index)):
        index = order_index[i]
        if index == target_value:
            out_index = i

    return out_index

# Create a function to compare propagated values to the ambient conditions
def compute_excess(freq_lst, point_lst, ambient_path, ambient_raster, weighting, results_dir, results_label, freq_fill, point_fill, temp_dir, summarize):
    
    weights_outpath = temp_dir + "weighted_ambients/"
    make(weights_outpath)
    weights_temp = weights_outpath + "temp/"
    make(weights_temp)
    temp_path = temp_dir + "temp_ambient_sum/"
    make(temp_path)
    
    # If a single file is not provided, assume it is a path and sum all existing ambient conditions
    overall_background = ambient_raster
    if ambient_path != "NA":
        # Create summary background layer
        # Weight sound from multiple frequencies
        bg_results_label = "background_%s" % results_label
        weight_sound_from_multiple_frequencies(freq_lst, weighting, ambient_path, weights_outpath, weights_temp, freq_fill, point_fill = "", point_id = "", outpath = "none", is_ambient = 1)    
        
        # Sum sound from multiple frequencies
        sum_sound_from_multiple_frequencies(freq_lst, weighting, weights_outpath, results_dir, bg_results_label, temp_path,
                                                freq_fill, point_fill = "", point_id = "")
    
        overall_background = results_dir + "%s_%s.tif" % (bg_results_label, weighting[0])
        
    # Check that overall_background exists
    if not os.path.exists(overall_background):
        raise ValueError("Background raster, %s, was not properly created" % overall_background)
    # Convert the background raster into a raster object for later steps
    overall_background = arcpy.Raster(overall_background)      

    #Overall excess
    if summarize == 1 or summarize == 3:
        overall_results = arcpy.Raster(results_dir + "%s_%s.tif" % (results_label, weighting[0]))
        overall_excess_file = results_dir + "excess_%s_%s.tif" % (results_label, weighting[0])
        trueValue = overall_results - overall_background
        # Make minimum excess 0
        overall_excess = arcpy.sa.Con((overall_results > overall_background), trueValue, 0)
        overall_excess.save(overall_excess_file)
    
    # Compute for each individual point if summarize == 2 or 4
    if summarize == 2 or summarize == 4:
        for point in point_lst:
            point_label = "_pt%s" % point.zfill(point_fill)
            overall_results = arcpy.Raster(results_dir + "%s_%s%s.tif" % (results_label, weighting[0], point_label))
            #Overall excess
            overall_excess_file = results_dir + "ex_%s_%s%s.tif" % (results_label, weighting[0], point_label)
            trueValue = overall_results - overall_background
            # Minimum excess at 0 dB
            overall_excess = arcpy.sa.Con((overall_results > overall_background), trueValue, 0)
            overall_excess.save(overall_excess_file)
                    
# Create a function to compare 
def compute_excess_by_freqs(point_dir, excess_freq_dir, freq_lst, freq_fill, ambient_path, model, point_id):

    # Check that individual frequency ambients exist
    ambient_path, ambient_raster = prep_ambient(ambient_path, freq_lst, freq_fill)

    # Ambient path may be modified to be NA, if only a global raster is given
    # Only compute ambient if ambient levels for each frequency are provided
    if ambient_path != "NA":
        
        # Drop the file-specific part of ambient path and add a trailing slash
        #ambient_path = os.path.split(ambient_path)[0] + "/"
    
        for i in range(len(freq_lst)):
            freq_s = freq_lst[i]
            freq_lbl = str(int(freq_s)).zfill(freq_fill) # create label with appropriate number of leading 0's
            background = ambient_path + "ambient%s" % freq_lbl
    
            # Subtract ambient conditions from propagation values to calculate excess noise propagation
            #**# NEED TO DO A BETTER JOB OF PASSING THESE FILE NAMES
            pr_file = point_dir + "pt%s_pr%s.tif" % (point_id, freq_lbl)
            if model == "spreadgis":
                pr_file = point_dir + "results/pr%s.tif" % freq_s #**# eventually this should be updated to freq_lbl            
            ex_file = excess_freq_dir + "ex%s.tif" % freq_lbl
            background = arcpy.sa.Raster(background)
            trueValue = pr_file - background 
            #**# Do we want to truncate to 0 here?
            ex = arcpy.sa.Con((pr_file > background), trueValue, 0)
            ex.save(ex_file)

# Function to convert from string summary inputs to 0 - 4 codes
def unpack_summary(in_text):
    
    if in_text == 'no summary':
        summarize = 0
    if in_text == 'points & frequencies':
        summarize = 1
    if in_text == 'frequencies only':
        summarize = 2
    if in_text == "points & frequencies, do not run model":
        summarize = 3
    if in_text == "frequencies only, do not run model":
        summarize = 4
    if in_text == "maximum across frequencies":
        summarize = 5

    return summarize
    
# Function to convert from YES & NO to 1 & 0
def unpack_ki(in_text):
    
    if str.upper(in_text) == "YES":
        outval = 1
    if str.upper(in_text) == "NO":
        outval = 0
    
    return outval
    
# Helper function to allow ambient inputs to be processed correctly
def prep_ambient(ambient_path, freq_lst, freq_fill):
    
    # Check if ambient is NA, if so, set ambient_raster to NA
    if ambient_path == "NA" or ambient_path == "" or ambient_path == "#":
        ambient_path = "NA"
        ambient_raster = "NA"
    else:
        
        # Check if there is a file ending
        file_ending_length = 0 # Default for ESRI GRID format
        if ambient_path[-4] == ".":
            file_ending_length = 4
            
        if ambient_path[-5] == ".":
            file_ending_length = 5
        
        is_path = 0
        # Check if the ending of the ambient input matches one of the frequencies.
        # If so, extract the entire path and assume there are other frequencies of interest
        # If not, assume it is an overall background raster
        for freq in freq_lst:    
            freq_lbl = str(int(freq)).zfill(freq_fill)
            len_lbl = len(freq_lbl) + file_ending_length
            # If ambient path is shorter than the label, it cannot be a match, and might crash the code
            if len(ambient_path) >= len_lbl:
                start = -len_lbl
                # 2 tests one if file_ending_length == 0, other for != 0. Need different syntax since -0 does not run to the end!
                if file_ending_length == 0:
                    if ambient_path[start:] == freq_lbl:
                        is_path = 1
                else:
                    end = -file_ending_length # Exclude the file ending from the check
                    if ambient_path[start:end] == freq_lbl:
                        is_path = 1
                    
        if is_path == 1:
            # Convert the path + file name to just the path (need to add a trailing slash)
            ambient_path = os.path.split(ambient_path)[0] + "/"
            ambient_raster = "NA"
        else:
            ambient_raster = ambient_path
            ambient_path = "NA"
        
    return(ambient_path, ambient_raster)
    
# Add a function to reclassify landcover to the correct format for the code
def reclass_landcover(landcover, landcover_type, model):

    if model == "spreadgis":
        # Add SPREADTYPE FIELD
        Add_SPREADTYPE(landcover, landcover_type)
    
    if model == "nmsimgis":
        Add_NMSIMTYPE(landcover, landcover_type)


# This directly modifies the landcover - think about the implications of this.
def Add_NMSIMTYPE(landcover, landcover_type):
    
    # Check if NMSIMTYPE field already exists, if so, use it.
    exists = check_for_field(landcover, "NMSIMTYPE")
    
    if exists == 1:
        arcpy.AddWarning("NMSIMTYPE FIELD ALREADY EXISTS, THE EXISTING FIELD WILL BE USED")
    else:
        
        if str.lower(landcover_type) == "nlcd":
            codeblock = '''def calcNLCD(x):
                newval = "" 
                if x == 1:
                    newval = "AK"
                if x == 11:
                    newval = "WATER"
                if x == 12:
                    newval = "SNOW"
                if x == 21:
                    newval = "URBAN1"
                if x == 22:
                    newval = "URBAN2"
                if x == 23:
                    newval = "URBAN3"
                if x == 24:
                    newval = "URBAN4"
                if x == 31:
                    newval = "BARREN"
                if x == 32:
                    newval = "SHORE"
                if x == 41 or x == 42 or x == 43:
                    newval = "FOREST"
                if x == 51 or x == 52:
                    newval = "SHRUB"
                if x == 71 or x == 72 or x == 81 or x == 82:
                    newval = "GRASS_HERB"
                if x == 73:
                    newval = "LICHEN"
                if x == 74:
                    newval = "MOSS"
                if x >= 90:
                    newval = "WETLAND"
                return(newval)'''
            
            arcpy.AddField_management(landcover, "NMSIMTYPE", "TEXT")
            arcpy.CalculateField_management(landcover, "NMSIMTYPE", "calcNLCD(!VALUE!)", "PYTHON", codeblock)

        # If landfire, apply the landfire classification (may need modification)
        elif str.lower(landcover_type) == "landfire":
            arcpy.AddWarning("CLASSIFICATION UNDER DEVELOPMENT, only codes 3001 through 3968 are currently supported by this tool, and please check your classified raster for mis-classifications!")
            codeblock = '''def calcLANDFIRE(x):
                newval = ""            
                if x == 3001 or x == 3002 or x == 3006 or x == 3153 or x == 3218 or x == 3219 or x == 3221:
                    newval = "BARREN"
                if x == 3222 or x == 3294 or x == 3295:
                    newval = "BARREN"
                if x >= 3011 and x <= 3062:
                    newval = "FOREST"
                if x == 3114 or x == 3152 or x == 3154 or x == 3159 or x == 3160 or x == 3180 or x == 3208:
                    newval = "FOREST"
                if x == 3229 or x == 3231 or x == 3900 or x == 3901 or x == 3902 or x == 3910 or x == 3911:
                    newval = "FOREST"
                if x == 3920 or x == 3921 or x == 3940 or x == 3941:
                    newval = "FOREST"
                if x >= 3064 and x <= 3113
                    newval = "SHRUB"
                if x >= 3115 and x <= 3127: #**# should 3117 be forest or shrub? Pine savannah. Not sure this is consistently classified. 3119 too.
                    newval = "SHRUB"
                if x >= 3211 and x <=3217:
                    newval = "SHRUB"
                if x == 3153 or x == 3220 or x == 3251 or x == 3252 or x == 3255 or x == 3259 or x == 3904:
                    newval = "SHRUB"
                if x == 3914 or x == 3923 or x == 3943 or x == 3960 or x == 3961:
                    newval = "SHRUB"
                if x >= 3135 and x <= 3146:
                    newval = "GRASS_HERB"
                if x >= 3181 and x <= 3182:
                    newval = "GRASS_HERB"
                if x == 3903 or x == 3913 or x == 3924 or x == 3944 or x == 3963 or x == 3964 or x == 3965:
                    newval = "GRASS_HERB"
                if x == 3966 or x == 3967 or x == 3968:
                    newval = "GRASS_HERB"
                if x == 3164:
                    newval = "WETLAND"
                if x == 3294:
                    newval = "BARREN"
                if x == 3292:
                    newval = "WATER"
                if x == 3296:
                    newval = "URBAN2"                    
                if x == 3297:
                    newval = "URBAN3"
                if x == 3298:
                    newval = "URBAN4"
                if x == 3299:
                    newval = "URBAN4" #**# Is this a reasonable classification of roads?
                if x == 3293:
                    newval = "SNOW"
                return(newval)'''

            arcpy.AddField_management(landcover, "SPREADTYPE", "TEXT")
            arcpy.CalculateField_management(landcover, "SPREADTYPE", "calcLANDFIRE(!VALUE!)", "PYTHON", codeblock) # calculate it to be equal to value, then reclassify the value field in the next step

    
        else:
            raise ValueError("Your chosen landcover type, %s, is not supported." % landcover_type +
            "Supported types are the custom format (landcover_type = 'default')"
            "or NLCD (landcover_type = 'nlcd')\n"
            "The default type uses the following case-sensitive landcover categories:\n"
            "WATER, SNOW, SHORE, WETLAND, FOREST, SHRUB, GRASS, BARREN,"
            "URBAN1, URBAN2, URBAN3, URBAN4, AK, LICHEN, MOSS.\n These must be in a column labeled 'NMSIMTYPE'.")
    
        #**# I'm happy to add other landcovers, as people can explain them 
        # to me! Some of these are pretty confusing, and I don't have the
        # time to figure out the details. I assumed they'd be straight-forward
        # like NLCD!
        '''        
        if str.lower(landcover_type) == "regap":
            landcover = reclass_regap_nmsimgis(landcover)
        
        if str.lower(landcover_type) == "geocover":
            landcover = reclass_geocover_nmsimgis(landcover)
        
        if str.lower(landcover_type) == "globecover":
            landcover = reclass_globecover_nmsimgis(landcover)
        '''

# Check for an existing field in an ArcGIS attribute table
def check_for_field(in_file, in_field):
    exists = 0
    fields = arcpy.ListFields(in_file)
    for field in fields:
        if str.upper(str(field.name)) == in_field:
            exists = 1
    
    return exists


# Summarize over the input layers
# Idea is to have a list of layers to be summarized, and a list of weights
# to be applied to each layer
#**# NEEDS TESTING
def custom_summary(list_of_layers, list_of_weights, results_path, results_name):

    # For a road input, first, sum together all the points in the road when
    # doing the model run. Then apply this tool with the weights for each vehicle
    # type.

    #check if list of layers is same length as list of weights
    if len(list_of_layers) != len(list_of_weights):
        raise ValueError("List of layers (%s rasters) must match list of weights (%s weigths)" % (len(list_of_layers), len(list_of_weights)))

    # Create a directory to hold temporary files
    # Consider using the scratch workspace here
    temppath = results_path + "temp/" 
    make(temppath)
 
    for i in xrange(len(list_of_layers)):
        layer = list_of_layers[i]
        traffic_weight = list_of_weights[i]
        # Convert layer from SEL (dB) to energy basis
        energy_layer = temppath + "energy%s.tif" % i
        bells_layer = temppath + "bells%s.tif" % i
        temp_total = temppath + "partial%s.tif" % i
        convert_dB_to_acoustic_energy(layer, bells_layer, energy_layer)
        
        # Multiply layer by traffic weight #**# Check that this is sufficient - it should be, but I thought Bruce sent me a formula, and I'm not sure it matched this.
        weighted_energy_layer = arcpy.sa.Times(energy_layer, traffic_weight)

        # Add to running layer total
        if i == 0:
            #new_total = arcpy.Raster(weighted_energy_layer)
            #new_total.save(temp_total)
            weighted_energy_layer.save(temp_total)
        if i != 0:
            prev_total = temppath + "partial%s.tif" % (i - 1)
            #new_total = arcpy.Raster(prev_total) + arcpy.Raster(weighted_energy_layer)
            new_total = arcpy.Raster(prev_total) + weighted_energy_layer
            new_total.save(temp_total)

    # Convert back to dB
    # Use the last temp_total raster, which following the loop contains the weighted sum
    weighted_db_layer = results_path + results_name
    arcpy.AddMessage(weighted_db_layer)
    convert_acoustic_energy_to_dB(temp_total, weighted_db_layer)        
    
    # Delete temporary files
    #shutil.rmtree(temppath) # Does not work - ArcGIS is retaining a schema lock on prev_total variable that prevents this.

    # Would it be worth returning weighted_db_layer?

# Convert from SPL to Leq
# Units don't matter, but must match between sound_source_duration and Leq_duration
#**# NEEDS TESTING
def convert_to_Leq_batch(list_of_rasters, sound_source_duration, Leq_duration, results_path):
    
    temppath = results_path + "temp/" 
    make(temppath)
            
    # For each raster:
    for i in xrange(len(list_of_rasters)):
        raster = list_of_rasters[i]  
        # Get raster base name
        Leq_file = results_path + os.path.basename(raster)
        # if base name has a file extension, drop it
        if Leq_file[-4] == ".":
            Leq_file = Leq_file[:-4]
        if Leq_file[-5] == ".":
            Leq_file = Leq_file[:-5]
        Leq_file = Leq_file + "Leq%s.tif" % Leq_duration
        
        # Convert the raster to an energy basis
        energy_layer = temppath + "energy%s.tif" % i
        bells_layer = temppath + "bells%s.tif" % i
        wel_file = temppath + "weighted_energy%s.tif" % i
        convert_dB_to_acoustic_energy(raster, bells_layer, energy_layer)
        
        # Multiply by duration
        weighted_energy_layer = energy_layer * sound_source_duration
        weighted_energy_layer.save(wel_file)
        
        # Divide by total time window
        Leq = weighted_energy_layer / float(Leq_duration)
        Leq.save(Leq_file)
        
# Convert a single, fluctuating noise source to a continuous Leq level
# List of rasters to allows for fluctating conditions
#**# NEEDS TESTING
def convert_to_Leq_fluctuating(list_of_rasters, list_of_sound_source_durations, Leq_duration, results_path, results_name):
    
    # check that list of rasters & list of durations match
    if len(list_of_rasters) != len(list_of_sound_source_durations):
        raise ValueError("Number of input rasters (%s) must have an equal number of durations (%s)" % (len(list_of_rasters), len(list_of_sound_source_durations)))
    
    temppath = results_path + "temp/" 
    make(temppath)
    
    # Check that sound source duration is equal to or shorter than sound source duration
    tot_duration = 0
    for i in list_of_sound_source_durations:
        tot_duration = tot_duration + list_of_sound_source_durations[i]
    
    # Check that the total duration is less than the Leq duration, and warn user if it is not
    if tot_duration > Leq_duration:
        warning_pt1 = "The sum of individual durations %s " % (tot_duration)
        warning_pt2 = "exceeds the total duration time %s. " % (Leq_duration)
        warning_pt3 = "If this was not intentional, (e.g., to compute an Leq8h) "
        warning_pt4 = "please ensure that the total duration of the source layers"
        warning_pt5 = " is less than the Leq duration time"
        arcpy.AddWarning("%s%s%s%s%s" % warning_pt1, warning_pt2, warning_pt3, warning_pt4, warning_pt5)
    
    # For each raster:
    for i in xrange(len(list_of_rasters)):
        raster = list_of_rasters[i]
        duration = list_of_sound_source_durations[i]
        
        # Convert the raster to an energy basis
        energy_layer = temppath + "energy%s.tif" % i
        bells_layer = temppath + "bells%s.tif" % i
        weighted_energy_layer = temppath + "weighted_energy%s.tif" % i
        partial_total = temppath + "partialtotal%s.tif" % i
        convert_dB_to_acoustic_energy(raster, bells_layer, energy_layer)
        
        # Multiply by duration
        weighted_energy_layer = energy_layer * duration
        weighted_energy_layer.save(weighted_energy_layer)
        
        if i == 0:
            new_total = arcpy.Raster(weighted_energy_layer)
            new_total.save(partial_total)
        if i != 0:
            prev_total = temppath + "partialtotal%s.tif" % (i - 1)
            new_total = arcpy.Raster(prev_total) + arcpy.Raster(weighted_energy_layer)
            new_total.save(partial_total)
  
   # Divide final total by total time window
    Leq = new_total / float(Leq_duration)
    Leq.save(results_name)


# What about calculating an L10, L50, L90? Isn't that the sound level present 10% of the time, 50% of the time and 90% of the time?
# But that would require information on the temporal distribution of the sound sources.
# And that requires something far more sophisticated.

# Add a function to give a warning if the model is expected to take a long time
def estimate_model_run_time(in_raster, CellSize, model):
    
    Extent = arcpy.Describe(in_raster).Extent
    
    # Estimate number of rows
    xmax = Extent.XMax   
    xmin = Extent.XMin
    rows = (xmax - xmin) / float(CellSize)
    
    # Estimate number of columns
    ymax = Extent.YMax
    ymin = Extent.YMin
    cols = (ymax - ymin) / float(CellSize)    
    
    # Estimate number of cells.    
    cells = rows * cols

    # Create messages to warn the user
    w1a = "YOUR MODEL RUN MAY TAKE > 30 MIN TO COMPLETE"
    w1b = "YOUR MODEL RUN MAY TAKE >24 HOURS TO COMPLETE"
    w2 = "You have %s cells in your landscape due to an extent of:" % cells
    w3 = "(xmin: %s; xmax: %s; ymin: %s; ymax: %s;) " % (xmin, xmax, ymin, ymax)
    w4 = "and cell size %s." % (CellSize)
    w5 = " Consider aborting the run and trying a smaller extent or a larger cell size."
    
    # Give estimate if model is SPreAD-GIS & warn user if # cells exceeds a threshold
    if model == "spreadgis":
        if cells > 10**9:
            arcpy.AddWarning("%s%s%s%s%s" % (w1a,w2,w3,w4,w5))
    
    # Give estimate if model is nmsimgis & warn user if # cells exceeds a threshold
    if model == "nmsimgis":
        if cells > 4 * 10**6:
            arcpy.AddWarning("%s%s%s%s%s" % (w1b,w2,w3,w4,w5))
        if cells > 300000 and cells <= 4 * 10**6:
            arcpy.AddWarning("%s%s%s%s%s" % (w1a,w2,w3,w4,w5))
           
# Check path length against path lengths used in the code to give the user an error
#**# Check how many characters are actually added in the current verison of the script
def check_path_length(base_dir):
    
    char_limit = 128
    added_char = 69 #**# May be an overestimate
    len_limit = char_limit - added_char
    
    if len(base_dir) > len_limit:
        e1 = "Unfortunately, this script uses files in ESRI GRD format, "
        e2 = "and files in GRD format can have a maximum path length of %s characters. " % char_limit
        e3 = "The code adds %s characters to the path length during processing, " % added_char
        e4 = "restricting the longest path length you can input here to %s characters. " % (len_limit)
        e5 = "your chosen path has %s characters." % (len(base_dir))
        raise ValueError("%s%s%s%s%s" % (e1,e2,e3,e4,e5))


# Function to facilitate the creation of many square extents. This function
# only makes one extent, but placed within a loop, it can create many exents
# Note the dependency on R!
def make_extent(R_exe, tbx_root, source_point, point_dir, extent_size, extent_file, utm_zone):
    R_script = tbx_root + "scripts/make_extents.R"
    
    # THIS ASSUMES THERE IS ONLY ONE SOURCE POINT IN THE SOURCE POINT FILE
    buffer_size = int(math.ceil(extent_size / 2.0)) # Divide extent size in half
    with arcpy.da.SearchCursor(source_point, ["SHAPE@XY"]) as cursor:
        count = 0
        for row in cursor:
            X, Y = row[0]   
            count += 1
            
    if count != 1:
        raise ValueError("The source point shapefile (%s) must have only one point in it" % source_point )

    args = "%s %s %s %s %s" % (X, Y, buffer_size, extent_file, utm_zone)
    command = "%s %s %s" % (R_exe, R_script, args)    
    subprocess.call(command)

# END OF FILE
