# -*- coding: utf-8 -*-
'''
Description: Run the core algorithms of the SPreAD-GIS model. Calculates
             decline in sound level due to spherical spreading,
             atmospheric absorption, vegetation, and wind loss and
             topographic effects; summarizes noise propagation patterns
             
Dependencies: ArcGIS >=10.3, Python 2.7, NumPy, soundprophlpr.py

@author Sarah E. Reed, Alexander "Sasha" Keyel <skeyel@gmail.com>

Copyright Information For spreadgishlpr.py:
Copyright (C) 2010 - 2017 Sarah E. Reed, A.C. Keyel <skeyel@gmail.com>
Upgraded to 10.X by Jessica Sushinski

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

# Import arcpy package
import arcpy, os, math, numpy, shutil, arcpy.sa

from soundprophlpr import seq

# Call the SPreAD-GIS tool
def spreadgis(run_directory, results_path, Sound_Source, point_id, n_points,Model_Extent,
              freq_lst, sound_level_lst, measurement_distance_lst, dem_ft, landcover,
              idir, CellSize, source_offset, receiver_offset,
              temp_s, hum_s, wind_dir, wind_sp, seas_cond, tbx_root,
              point_fill, freq_fill, my_times, my_time_labels, point_counter=1, keep_intermediates = 0, use_old_barrier = 1):

    # point_counter comes from soundprophlpr. Here 1 prints messages, all other values do not.

    # Set up for function call
    # Local variables
    eucdist_ft = idir + "eucdist_ft"
    eucdist = idir + "eucdist"
    eucdir = idir + "eucdir"
    sound_src = idir + "sound_src"

    # CALCULATE EUCLIDEAN DISTANCE & DIRECTION
    euclidean_dist_dir(Sound_Source, sound_src, eucdist, eucdist_ft, eucdir, CellSize)

    # CALCULATE MEAN ELEVATION
    elev = calculate_mean_elevation(dem_ft, sound_src, idir)
    #my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CALCULATE MEAN ELEVATION %s Hz" % freq_s)

    # If using the new barrier condition, set up barrier info and vegetation loss
    # once outside of frequency loop.
    BPD = "NA"
    max_veg_loss = "NA"
    if use_old_barrier == 0:
        # CALCULATE BARRIER PATH DISTANCE & MAX VEG COST
        BPD, max_veg_loss = calculate_barrier_path_distance_and_vegmax(Sound_Source, dem_ft, landcover, eucdist_ft, source_offset, receiver_offset, idir)

    # Loop through frequencies
    for i in range(len(freq_lst)):
        freq_s = freq_lst[i]
        Sound_Level = sound_level_lst[i]
        Measurement_Distance = measurement_distance_lst[i]
        
        # Add leading zeros for better sorting
        point_lbl = str(point_id).zfill(point_fill)        
        freq_lbl = str(int(freq_s)).zfill(freq_fill)
            
        my_times, my_time_labels = NoisePropagationOnePoint(run_directory, Sound_Source, sound_src, Model_Extent, freq_s,
                                 Sound_Level, Measurement_Distance, dem_ft, landcover,
                                 elev, eucdist, eucdist_ft, eucdir, 
                                 temp_s, hum_s, wind_dir, wind_sp, seas_cond,
                                 tbx_root, keep_intermediates, my_times, my_time_labels, point_counter, use_old_barrier, BPD, max_veg_loss)
                                 
                                 
        Raster1a = run_directory + "results/pr" + freq_s + ".tif"
        Raster1b = results_path + "pt" + point_lbl + "_pr" + freq_lbl + ".tif"
        arcpy.CopyRaster_management(Raster1a, Raster1b)

    return(my_times, my_time_labels)

# Main function for noise propagation for a single point source
def NoisePropagationOnePoint(base_directory, Sound_Source, sound_src, Model_Extent, freq_s, Sound_Level,
                                Measurement_Distance, dem_ft, landcover,
                                elev, eucdist, eucdist_ft, eucdir, 
                                temp_s, hum_s,
                                wind_dir_correct, wind_sp, seas_cond,
                                tbx_root, keep_intermediates, my_times, my_time_labels, point_counter, use_old_barrier, BPD, max_veg_loss): # timelog, delete_existing_timelog):
      
    # For no timing, set timelog to "none" and delete_existing_timelog to 1
    
    if use_old_barrier == 1:
        arcpy.AddWarning("Use of the old barrier approach for SPreAD-GIS has been deprecated. It has been kept as the default to maintain backwards compatiblity with earlier code. However, use of the new barrier calculations is recommended.")

    # determine whether timing is desired
    do_timing = 1
    if my_time_labels == "NA":
        do_timing = 0

    # compute set up time
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "Initialization %s Hz" % freq_s)


    # CONVERT FROM METRIC TO OLD ENGLISH FOR HISTORICAL REASONS
    temp_s = float(temp_s) * 1.8 + 32 # convert to F    
    wind_sp = float(wind_sp) * 0.621371 #Conversion values from Google convert to mph
    Measurement_Distance = float(Measurement_Distance) * 3.28084  # conversion value from Google

    # Convert wind direction from direction wind is blowing FROM to direction wind is blowing TO
    wind_dir = convert_wind_dir(wind_dir_correct)
        
    # Check out the spatial analyst extension
    arcpy.CheckOutExtension("Spatial")
    
    # Check that inputs are valid
    input_check(freq_s)    
    
    intermediates_dir = base_directory + "intermediate_%s/" % freq_s    
    
    # Delete and re-create required directories
    shutil.rmtree(intermediates_dir, ignore_errors=True)
    try:
        make_dirs(intermediates_dir) # requires os package
    except:
        raise ValueError("There was a problem deleting and re-creating the model directories.\n"
                         "Please try again, or if running in python, try opening a new Python console.\n"
                         "Please also look for any open files that may be preventing old model runs from"
                         " being deleted.")

    results_dir = base_directory + "results/" 
    prelim_dir = results_dir + "prelim_%s/" % freq_s
    if not os.path.exists(prelim_dir):
        os.makedirs(prelim_dir) #Note that this is different from os.mkdir, whick does not work recursively!

    shutil.copy(tbx_root + "tables/Table_13.dbf", intermediates_dir + "wind/Table_13.dbf")
    
    # Delete intermediate datasets, if they exist
    arcpy.env.overwriteOutput = True
    
    # Issue warning that data from prior model runs will be deleted
    if point_counter == 1:
        arcpy.AddWarning("Warning: Intermediate data from prior runs of SPreAD-GIS were deleted and results for the " + freq_s + " Hz frequency band will be overwritten ...")

    '''
    # Clip dem to desired extent "NO_MAINTAIN_EXTENT" prevents resampling the rasters & may result in a slightly larger final extent
    # Convert input from meters to feet
    dem_clip = intermediates_dir + "dem_clip"
    #arcpy.Clip_management(dem, Model_Extent, dem_ft, "", "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #ClippingGeometry option causes it to crash when only a bounding rectangle is input.
    #arcpy.Clip_management(dem, '295530.76 4103795.76 295743.08 4104008.08', dem_ft, "", "", "", "NO_MAINTAIN_EXTENT")
    arcpy.Clip_management(dem, Model_Extent, dem_clip, "", "", "", "NO_MAINTAIN_EXTENT")
    '''    
    
    # Convert dem to feet #**# MOVED TO SOUNDPROPHLPR.PY & DIRECTORY CHANGED
    #dem_ft = intermediates_dir + "dem_ft"
    #this_dem = arcpy.Raster(dem_clip) * 3.28084 # Conversion values taken from Google on 2016-03-16
    #this_dem.save(dem_ft)    

    # RasterExtent & CellSize are used below, but ideally, these parameters should
    # be inherited from soundprophlpr, where they are first defined.
    # Set general parameters for rasters
    dscRD = arcpy.Describe(dem_ft)
    RasterExtent = dscRD.Extent
    CellSize = dscRD.MeanCellHeight

    # Set Environment settings # Now set here for all tools instead of in each tool.
    arcpy.env.extent = RasterExtent
    arcpy.env.cellSize = CellSize
    arcpy.env.snapRaster = arcpy.Raster(dem_ft)

    
    # Copy elevation dataset to intermediate folder
    #dem_ft = intermediates_dir + "dem_ft"
    #arcpy.CopyRaster_management(dem, dem_ft)
    
    # compute set up time
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "SETUP %s Hz" % freq_s)

    #my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CALCULATE EUCLIDEAN DISTANCE & DIRECTION %s Hz" % freq_s)
                
    # CREATE A RASTER OF THE SOUND SOURCE LEVEL at specified frequency AT ALL POINTS (i.e. with no sound loss)
    source = baseline(Sound_Level, freq_s, intermediates_dir, CellSize, RasterExtent)    
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CREATE CONSTANT SOURCE RASTER %s Hz" % freq_s)
    
    # CALCULATE SPHERICAL SPREADING LOSS #
    print_message("Calculating spherical spreading loss ...", point_counter)    
    ssl = spherical_spreading_loss(eucdist_ft, Measurement_Distance, Sound_Level, freq_s, intermediates_dir)
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CALCULATE SPHERICAL SPREADING LOSS %s Hz" % freq_s)
    
    # CALCULATE ATMOSPHERIC ABSORPTION LOSS #
    aal = atmospheric_absorption_loss(elev, eucdist_ft, hum_s, temp_s, freq_s, intermediates_dir, CellSize, RasterExtent, point_counter)
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CALCULATE ATMOSPHERIC ABSORPTION LOSS %s Hz" % freq_s)
    
    if use_old_barrier == 1:
        # CALCULATE FOLIAGE & GROUND COVER LOSS #
        print_message("Calculating foliage and ground cover loss ...", point_counter)    
        veg = foliage_groundcover_loss(eucdist_ft, landcover, Sound_Source, freq_s, intermediates_dir, RasterExtent, CellSize)
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CALCULATE FOLIAGE & GROUND COVER LOSS %s Hz" % freq_s)
    
    # CALCULATE UPWIND AND DOWNWIND LOSS #
    print_message("Calculating upwind and downwind loss for wind blowing from %s degrees and wind speed of %0.1f k/hr on a %s ..." % (wind_dir_correct, wind_sp * 1.60934, seas_cond), point_counter)    
    wind = windloss(wind_dir, wind_sp, seas_cond, eucdist_ft, eucdir, freq_s, intermediates_dir, CellSize, RasterExtent, sound_src, dem_ft)
    my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "CALCULATE UPWIND AND DOWNWIND LOSS %s Hz" % freq_s)

    # Option to use original SPreAD-GIS barrier calculation
    if use_old_barrier == 1:
    
        # Delineate barrier for calculation of barrier loss & topographic zones
        print_message("Calculating decline in sound levels due to barrier loss ...", point_counter)
        barrier, ground = delineate_barrier(sound_src, dem_ft, intermediates_dir)    
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "DELINIATE BARRIER & COMPUTE TOPO ZONES %s Hz" % freq_s)
        
        # CALCULATE TOPOGRAPHIC EFFECTS AND BARRIER LOSS #
        bar = topographic_barrier_effects(barrier, elev, eucdist_ft, eucdir, dem_ft, freq_s, intermediates_dir, CellSize, RasterExtent)
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "TOPOGRAPHIC EFFECTS AND BARRIER LOSS %s Hz" % freq_s)
    
        # CALCULATE LOCATIONS WHERE GROUND, ATMOSPHERIC, AND BARRIER EFFECTS DOMINATE
        print_message("Identifying areas where ground, barrier, and atmospheric effects dominate ...", point_counter)
        topo_zones = calculate_topozones(Sound_Source, dem_ft, ground, intermediates_dir, sound_src)
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "LOCATIONS WHERE GROUND & ATMOSPHERIC & BARRIER EFFECTS DOMINATE %s Hz" % freq_s)

        # CALCULATE SUMMARY NOISE PROPAGATION PATTERNS # 
        print_message("Calculating final noise propagation patterns ...", point_counter)
        compute_noise_propagation(freq_s, source, ssl, aal, veg, wind, bar, topo_zones, results_dir, prelim_dir, sound_src, keep_intermediates, intermediates_dir)                           
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "SUMMARY NOISE PROPAGATION PATTERNS %s Hz" % freq_s)


    if use_old_barrier == 0:
        # Calculate barrier effect
        bar = barrier_effects_v2(BPD, freq_s, intermediates_dir)

        # Check if barrier + wind exceeds 25 dB, if so, cap at 25 dB
        barwind = check_barrier_wind(bar, wind, freq_s, intermediates_dir)

        # CALCULATE SUMMARY NOISE PROPAGATION PATTERNS # 
        print_message("Calculating final noise propagation patterns ...", point_counter)
        compute_noise_propagation_v2(freq_s, source, ssl, aal, max_veg_loss, barwind, results_dir, prelim_dir, sound_src, keep_intermediates, intermediates_dir)                           
        my_times, my_time_labels = sg_get_time(do_timing, my_times, my_time_labels, "SUMMARY NOISE PROPAGATION PATTERNS %s Hz" % freq_s)
        

    return(my_times, my_time_labels)

# convert wind dir #**# the reverse conversion is done elsewhere - this would be something that could be made more efficient!
def convert_wind_dir(wind_dir):
    
    wind_dir = float(wind_dir)    
    
    if wind_dir < 180 and wind_dir >= 0:
        wind_dir_out = wind_dir + 180
    elif wind_dir >= 180 and wind_dir <=360:
        wind_dir_out = wind_dir - 180
    else:
        raise ValueError("%s wind value is not valid. Wind direction must be in range 0 - 360 degrees")
        
    return(wind_dir_out)


# Function to assist in rapid generation of directories
def make_dirs(base_directory):
    '''makes the base directory and aal, ssl, topo, veg, wind subfolders'''
    os.makedirs(base_directory) # makedirs is recursive
    os.mkdir(base_directory + "aal")
    os.mkdir(base_directory + "ssl")
    os.mkdir(base_directory + "topo")
    os.mkdir(base_directory + "veg")
    os.mkdir(base_directory + "wind")
   

# Function to check that inputs are valid and within reasonable ranges
# Note that single-point & multipoint in the toolbox have constraints, but
# multi-frequency tables do not. Now only values in the Toolbox drop-down are
# allowed
def input_check(freq_s):
    '''Input frequency'''
    freq = float(freq_s)
    
    #if freq < 400:
    #    raise ValueError("SPreAD-GIS currently does not support frequencies lower than 400")
    #if freq > 2000:
    #    raise ValueError("SPreAD-GIS currently does not support frequencies greater than 2000")
    is_ok = 0

    if freq == 125 or freq == 160 or freq == 200 or freq == 250 or freq == 315:
        is_ok = 1
    
    if freq == 400 or freq == 500 or freq == 630 or freq == 800 or freq == 1000:
        is_ok = 1
        
    if freq == 1250 or freq == 1600 or freq == 2000:
        is_ok = 1
    
    if is_ok != 1:
        raise ValueError("SPreAD-GIS currently does not support frequency %s" % freq)
    
#**# Moved to soundprophlpr to try to avoid schema locks
'''
# CALCULATE EUCLIDEAN DISTANCES & DIRECTIONS #
def euclidean_dist_dir(Sound_Source, sound_src, eucdist, eucdist_ft, eucdir, CellSize):

    # Check if "ONE" field already exists, if so, delete it
    fields = arcpy.ListFields(Sound_Source)
    for field in fields:
        if field.name == "ONE":
            arcpy.DeleteField_management(Sound_Source, "ONE")

    # Convert sound source to raster dataset
    arcpy.AddField_management(Sound_Source, "ONE", "SHORT")
    arcpy.CalculateField_management(Sound_Source, "ONE", "1", "PYTHON")
    #tempEnvironment0 = arcpy.Extent()
    #arcpy.env.extent = Model_Extent
    arcpy.PointToRaster_conversion(Sound_Source, "ONE", sound_src, "", "", CellSize)
    #arcpy.env.extent = tempEnvironment0
    
    # Calculate Euclidean allocation, distance, and direction around sound source
    eucDirection = arcpy.sa.EucDirection(sound_src, "", CellSize, eucdist)
    eucDirection.save(eucdir)
    
    # Convert Euclidean distance from meters to feet
    times_result = arcpy.sa.Times(eucdist, 3.2808)
    times_result.save(eucdist_ft)
'''

# CALCULATE SOUND PROPAGATION IF THERE WAS NO LOSS #
def baseline(Sound_Level, freq_s, intermediates_dir, CellSize, RasterExtent):
    NoReduction = intermediates_dir + "source" + freq_s
    source = arcpy.sa.CreateConstantRaster(Sound_Level, "FLOAT", CellSize, RasterExtent)
    source.save(NoReduction)
    return(source)

# CALCULATE SPHERICAL SPREADING LOSS #
def spherical_spreading_loss(eucdist_ft, Measurement_Distance, Sound_Level, freq_s, intermediates_dir):
    ''' Documentation goes here '''

    # Add notification message

    # Divide Euclidean distance by measurement distance (X/y)
    xy_grid = intermediates_dir + "ssl/xy_grid"
    #Measurement_Distance_i = int(Measurement_Distance)
    #divide_result = arcpy.sa.Divide(eucdist_ft, Measurement_Distance_i)
    divide_result = arcpy.sa.Divide(eucdist_ft, Measurement_Distance)
    divide_result.save(xy_grid)
    
    # Reclassify 0 cell at center to be Spherical Spreading Loss at 1 m # Reclassify doesn't work, using Con instead
    ssl_patch_file = intermediates_dir + "ssl/ssl_patch.tif"    
    #reclass_value = 3.28084 / Measurement_Distance # conversion from meters into ft from Google, float division because the numerator is a float.
    # Using 1 m (arbitrarily chosen) This may be more or less than NMSIMGIS, depending on the distance between the source point and the cell of origin's center
    reclass_value = 3.28084 / Measurement_Distance # 1 m to ft, conversion from meters into ft from Google, float division because the numerator is a float.
    ssl_patch = arcpy.sa.Con((divide_result > 0), divide_result, reclass_value) # replace 0 with 1 meter distance
    ssl_patch.save(ssl_patch_file)
    
    # Calculate spherical spreading loss
    ssl_reclass = intermediates_dir + "ssl/ssl_reclass"
    #ssl = 8.64 * arcpy.sa.Ln(ssl_patch) + 0.149
    ssl = 20 * arcpy.sa.Log10(ssl_patch) # Change to equation used by NMSIMGIS & ISO 9613-2
    ssl.save(ssl_reclass)
    
    return(ssl)
    

# CALCULATE ATMOSPHERIC ABSORPTION LOSS #
def atmospheric_absorption_loss(elev, eucdist_ft, hum_s, temp_s, freq_s, intermediates_dir, CellSize, RasterExtent, point_counter):

    # Convert input parameters to numbers
    rh = float(hum_s)
    temp_f = float(temp_s)
    freq = float(freq_s)

    # Conversions for calcualtions 
    elev_m = elev / 3.28084
    temp_c = (temp_f - 32) * 5 / 9 # Convert to Celsius (note: no problem with floor division, because temp_f is a float)
    temp_k = temp_c + 273.15    # Convert to Kelvins
    
    # Update message
    message = "Calculating atmospheric absorption loss for an elevation of %0.2f m, air temperature of %0.1f C degrees , and %0.1f" % (elev_m, temp_c, rh) + "% humidity ..."
    print_message(message, point_counter)
    
    # Calculate atmospheric absorption coefficient
    alpha = atmospheric_absorption_loss_core(elev_m, rh, temp_k, freq)
    alpha_ft = alpha / 3.28084
    
    # Create constant raster of atmospheric absorption coefficient
    OutRaster1 = intermediates_dir + "aal/aac" + freq_s
    ccRaster = arcpy.sa.CreateConstantRaster(alpha_ft, "FLOAT", CellSize, RasterExtent)
    ccRaster.save(OutRaster1)
    
    # Calculate atmospheric absorption loss
    OutRaster2 = intermediates_dir + "aal/aal" + freq_s
    aal = arcpy.sa.Times(eucdist_ft, OutRaster1)
    aal.save(OutRaster2)
    
    return(aal)    

# Main calculations for atmospheric absorption. Also used by soundprophlpr.py for screening weather conditions
def atmospheric_absorption_loss_core(elev_m, rh, temp_k, freq):
       # Calculate atmospheric absorption coefficient using ANSI S1.26-1995 standard
    
    # Convert elevation to atmospheric pressure
    p_a = 101.325 * (1 - (2.25577 * (10 ** (-5)) * elev_m)) ** 5.25588
    
    # Convert relative humidity to molar concentration of water vapor
    C = (-6.8346 * ((273.16 / temp_k) ** 1.261)) + 4.6151
    psat_pr = 10 ** C
    h = (rh) * (psat_pr) * ((p_a / 101.325) ** (-1))
    # Calculate derived values for subsequent equations
    pa_pr = p_a / 101.325
    T_Tr = temp_k / 293.15
    e = 2.7182818284
    
    # Calculate frO (equation 3)
    frO = ((pa_pr) * ((24 + (4.04 * 10000)) * h) * (0.02 + h)) / (0.391 + h)
    
    # Calculate frN (equation 4)
    frN = pa_pr * (T_Tr ** (-0.5)) * (9 + (280 * h * (e ** (-4.170 * ((T_Tr ** (-0.33333)) - 1)))))
    
    # Calculate alpha (equation 5)
    term1 = 1.84 * (10 ** (-11)) * (pa_pr ** (-1)) * (T_Tr ** 0.5)
    term2 = (T_Tr ** (-2.5)) * (0.01275 * (e ** (-2239.1 / temp_k)) * (frO / ((frO ** 2) + (freq ** 2))))
    term3 = 0.1068 * (e ** (-3352 / temp_k)) * (frN / ((frN ** 2) + (freq ** 2)))
    alpha = 8.686 * (freq ** 2)*(term1 + term2 + term3)
    
    return(alpha)

# CALCULATE FOLIAGE & GROUND COVER LOSS #
def foliage_groundcover_loss(eucdist_ft, landcover, Sound_Source, freq_s, intermediates_dir, RasterExtent, CellSize):
  
    # Create a constant raster of value = 1
    constant1 = intermediates_dir + "veg/constant1" 
    ccRaster = arcpy.sa.CreateConstantRaster("1", "INTEGER", CellSize, RasterExtent)
    ccRaster.save(constant1)
       
    # Take maximum value of the Euclidean distance grid
    eucdist_max = intermediates_dir + "veg/eucdist_max"
    zStats = arcpy.sa.ZonalStatistics(constant1, "VALUE", eucdist_ft, "MAXIMUM", "DATA")
    zStats.save(eucdist_max)
    
    # Subtract maximum distance from the Euclidean distance grid and take its absolute value
    eucdist_z = intermediates_dir + "veg/eucdist_z"
    #raster1 = arcpy.sa.Raster(intermediates_dir + "eucdist_ft")
    eucdist_ft_raster = arcpy.sa.Raster(eucdist_ft)
    raster2 = arcpy.sa.Raster(intermediates_dir + "veg/eucdist_max")
    calc_result = arcpy.sa.Abs(eucdist_ft_raster - raster2)
    calc_result.save(eucdist_z)
    
    # Reclassify land cover datasets by foliage and ground cover loss rates
    veg_lossrate = intermediates_dir + "veg/veg_lossrate"
    reclass = arcpy.sa.Reclassify(landcover, "SPREADTYPE", "CON 501;HWD 662;SHB 101;HEB 101;BAR 0;WAT 0;URB 0", "DATA")
    reclass.save(veg_lossrate)
    
    # Buffer sound source point (to reduce errors in path distance calculations)
    sound_src_veg = intermediates_dir + "veg/sound_src_veg.shp"
    arcpy.Buffer_analysis(Sound_Source, sound_src_veg, "%s Meters" % CellSize, "FULL", "ROUND", "ALL")
    
    # Calculate the cost raster for vegetation loss by distance
    veg_cost = intermediates_dir + "veg/veg_cost"
    veg_lossrate = arcpy.sa.Raster(intermediates_dir + "veg/veg_lossrate")
    #eucdist_ft = arcpy.sa.Raster(intermediates_dir + "eucdist_ft")
    #eucdist_ft_raster = arcpy.sa.Raster(eucdist_ft) # Defined above (formerly as raster1)
    scale_factor = CellSize * 3.28084 # Conversion from m to ft taken from Google
    #Note: below the / does normal division as at least one input includes a floating point
    calc_result = ((veg_lossrate / eucdist_ft_raster) / scale_factor) + 1 
    calc_result.save(veg_cost)
    
    # Calculate path distance around buffered source point (without vegetation loss)
    pathdist1 = intermediates_dir + "veg/pathdist1"
    #b1 = arcpy.sa.HfBinary("1", "45")
    #b2 = arcpy.sa.VfBinary("1", "-30", "30")
    pathDist = arcpy.sa.PathDistance(sound_src_veg, "", eucdist_z, "", "BINARY 1 45", "", "BINARY 1 -30 30") # This is equivalent to using arcpy.sa.HfBinary(1, 45) and arcpy.sa.VfBinary(1, -30, 30)
    #pathDist = arcpy.sa.PathDistance(sound_src_veg, "", eucdist_z, "", b1, "", b2)
    pathDist.save(pathdist1)
    
    # Calculate path distance around buffered source point (with vegetation loss)
    pathdist2 = intermediates_dir + "veg/pathdist2"
    pathDist = arcpy.sa.PathDistance(sound_src_veg, veg_cost, eucdist_z, "", "BINARY 1 45", "", "BINARY 1 -30 30") # This is equivalent to using arcpy.sa.HfBinary(1, 45) and arcpy.sa.VfBinary(1, -30, 30)
    #pathDist = arcpy.sa.PathDistance(sound_src_veg, veg_cost, eucdist_z, "", b1, "", b2)
    pathDist.save(pathdist2)
    
    # Use ArcGIS 9.3 inputs to test compatibility
    #pathdist1 = "C:/SPreAD-GIS/backup/Arc9_3/intermediate/veg/pathdist1"    
    #pathdist2 = "C:/SPreAD-GIS/backup/Arc9_3/intermediate/veg/pathdist2"
    
    # Subtract path distance 1 from path distance 2 to calculate vegetation loss
    path2minuspath1 = intermediates_dir + "veg/veg_pre" + freq_s # + ".tif"
    veg_pre = arcpy.sa.Minus(pathdist2, pathdist1)
    veg_pre.save(path2minuspath1)
    
    # Pathdist2 gives NoData for cell of origin, need to patch this to be 0 to
    # prevent an inaccurate assessment at the origin cell
    OutRaster1 = intermediates_dir + "veg/veg" + freq_s #+ ".tif"
    veg = arcpy.sa.Con(arcpy.sa.IsNull(veg_pre), 0, veg_pre) # replace NoData with 0
    veg.save(OutRaster1)
    
    return(veg)

# Convert seasonal conditions for windloss calculations
def convert_seasonal_conditions(seas_cond):
    seascond1 = "clear, calm summer day"
    seascond2 = "clear, calm winter day"
    seascond3 = "clear, calm summer night"
    seascond4 = "clear, calm winter night"
    seascond5 = "clear, windy summer day"
    seascond6 = "clear, windy winter day"
    seascond7 = "clear, windy summer night"
    seascond8 = "clear, windy winter night"
    seascond9 = "cloudy, calm"
    seascond10 = "cloudy, windy"
    
    if seas_cond == seascond1:
       phi = 180
    elif seas_cond == seascond2:
       phi = 180
    elif seas_cond == seascond3:
       phi = 0
    elif seas_cond == seascond4:
       phi = 180
    elif seas_cond == seascond5:
       phi = 144
    elif seas_cond == seascond6:
       phi = 144
    elif seas_cond == seascond7:
       phi = 62
    elif seas_cond == seascond8:
       phi = 70
    elif seas_cond == seascond9:
       phi = 90
    elif seas_cond == seascond10:
       phi = 90
    
    return phi

# CALCULATE UPWIND AND DOWNWIND LOSS #
def windloss(wind_dir, wind_sp, seas_cond, eucdist_ft, eucdir, freq_s, intermediates_dir, CellSize, RasterExtent, sound_src, dem_ft):

    # Environment settings
    #arcpy.env.extent = RasterExtent
    #arcpy.env.cellSize = CellSize
    #arcpy.env.snapRaster = arcpy.Raster(dem_ft)
    # Model change: Made no wind loss if no wind! 11/10/2016
    wind_sp_f = float(wind_sp)
    if wind_sp_f < 0:
        raise ValueError("Wind speed cannot be less than 0! Wind speed was %s" % wind_sp_f)

    if wind_sp == 0:
        # Create a constant raster of 0
        wind_file = intermediates_dir + "wind/win" + freq_s
        wind = arcpy.sa.CreateConstantRaster(0, "INTEGER", CellSize, RasterExtent)
        wind.save(wind_file)

    else:
        # Convert frequency to numeric
        freq = float(freq_s) 
        
        # Convert seasonal conditions to phi
        phi = convert_seasonal_conditions(seas_cond)
    
        # Create a constant raster of phi value
        phi_grid_file = intermediates_dir + "wind/phi"
        phi_grid = arcpy.sa.CreateConstantRaster(phi, "FLOAT", CellSize, RasterExtent)
        phi_grid.save(phi_grid_file)
            
        # Compute sound propagation angles away from the source
        plus_result_file = intermediates_dir + "wind/plus1"
        minus_result_file = intermediates_dir + "wind/minus1"
        prop_angle_file = intermediates_dir + "wind/prop_angle"
    
        plus_result = arcpy.sa.Plus(eucdir, 180)
        plus_result.save(plus_result_file)
        minus_result = arcpy.sa.Minus(eucdir, 180)
        minus_result.save(minus_result_file)
        prop_angle = arcpy.sa.Con(eucdir, plus_result, minus_result, "VALUE < 180")
        prop_angle.save(prop_angle_file)
        
        # Subtract prevailing wind direction from 180 and add to sound propagation angles
        wind_dir_f = float(wind_dir)
        diff = 180 - wind_dir_f
        conRaster1_file = intermediates_dir + "wind/conraster1.tif"
        conRaster1 = arcpy.sa.Plus(prop_angle, diff)
        conRaster1.save(conRaster1_file)
        
        # Re-classify propagation values that are less than zero or greater than 360
        trueValue = conRaster1 - 360
        trueValue_2 = conRaster1 + 360
        conresult_file = intermediates_dir + "wind/conresult1.tif"
        con_result = arcpy.sa.Con((conRaster1 >= 360), trueValue, (arcpy.sa.Con((conRaster1 < 0), trueValue_2, conRaster1)))
        con_result.save(conresult_file)
        
        # Re-classify wind angle values to range between 0 and 180
        trueValue = 360 - con_result
        wind_ang_file = intermediates_dir + "wind/wind_ang"
        wind_ang = arcpy.sa.Con((con_result > 180), trueValue, con_result)
        wind_ang.save(wind_ang_file)
        
        # Identify upwind and downwind areas
        updownwind_file = intermediates_dir + "wind/updownwind"
        upwind_file = intermediates_dir + "wind/upwind"
        downwind_file = intermediates_dir + "wind/downwind"
    
        updownwind = arcpy.sa.Minus(phi_grid, wind_ang)    
        updownwind.save(updownwind_file)
        upwind = arcpy.sa.Con(updownwind, updownwind, 0, "VALUE > 0")
        upwind.save(upwind_file)
        downwind = arcpy.sa.Con(updownwind, 1, 0, "VALUE <= 0")
        downwind.save(downwind_file)
        
        # Calculate upwind loss
        falseValue = 5.7642 * arcpy.sa.Ln(upwind) + 2.5664
        upwind_loss_file = intermediates_dir + "wind/upwind_loss"
        upwind_loss = arcpy.sa.Con((upwind <= 0), 0, (arcpy.sa.Con((upwind >= 50), 25, falseValue)))
        upwind_loss.save(upwind_loss_file)
        
        # Calculate shadow zone correction to upwind loss
        x_d_file = intermediates_dir + "wind/x_d"
        upwind_szf_file = intermediates_dir + "wind/upwind_szf"
        upwind_loss_c_file = intermediates_dir + "wind/upwind_loss_c"
    
        table_13 = intermediates_dir + "wind/Table_13.dbf"
        table_freq = intermediates_dir + "wind/table_freq.dbf"
        
        #if wind_sp_f > 0:
        d = 375 * (wind_sp_f ** (-0.85))
        x_d = arcpy.sa.Divide(eucdist_ft, d)
        x_d.save(x_d_file)
        if freq <= 125:
            freq_w = 125
        elif freq >= 2000:
            freq_w = 2000
        else:
            freq_w = freq
        arcpy.TableSelect_analysis(table_13, table_freq, '"FREQ" = ' + str(int(freq_w)))
        upwind_szf = arcpy.sa.ReclassByTable(x_d, table_freq, "FROM_","TO","SZF_100") # This step is coming out wrong! But just in the tool
        upwind_szf.save(upwind_szf_file)
        #upwind_szf = arcpy.sa.ReclassByTable(x_d, table_freq, "FROM_","TO","SZF_100", "NODATA") # This step is coming out wrong! But just in the tool. Wait - it appears to have been working off a corrupted x_d file.
        upwind_loss_c = (upwind_loss * upwind_szf) / 100 #**# This is based on a 100 ft cell size. Does this need to be rescaled to account for the actual cell size? I think the model already accounts for distance above, so a second adjustment for distance is not warranted.
        upwind_loss_c.save(upwind_loss_c_file)
    
        #else:
        #   #arcpy.CopyRaster_management(upwind_loss, upwind_loss_c)
        #   upwind_loss_c = upwind_loss
        
        # Calculate downwind loss
        freq_dist_file = intermediates_dir + "wind/freq_dist"
        freq_dist = arcpy.sa.Times(eucdist_ft, freq)
        freq_dist.save(freq_dist_file)
        falseValue = downwind * (4.2598 * (arcpy.sa.Ln(freq_dist)) - 55.014)
        downwind_loss_file = intermediates_dir + "wind/downwind_loss"
        downwind_loss = arcpy.sa.Con((freq_dist <= 406237), 0, falseValue)
        downwind_loss.save(downwind_loss_file)
        
        # Combine upwind and downwind loss
        wind_loss_file = intermediates_dir + "wind/wind_loss"
        wind_loss = arcpy.sa.Plus(upwind_loss_c, downwind_loss)
        wind_loss.save(wind_loss_file)
        
        # Smooth transition between upwind and downwind areas and calculate wind loss
        # NOTE: The DATA option on FocalStatistics is not working properly here
        # (but this is no longer an issue)
        #BUG -000082753 - need to patch Arc10.3.1
        wind_file = intermediates_dir + "wind/win" + freq_s
        wind = arcpy.sa.FocalStatistics(wind_loss, "Rectangle 9 9 CELL", "MEAN", "DATA")
        wind.save(wind_file)
    
    return(wind)

# CALCULATE TOPOGRAPHIC EFFECTS AND BARRIER LOSS #
def topographic_barrier_effects(barrier, elev, eucdist_ft, eucdir, dem_ft, freq_s, intermediates_dir, CellSize, RasterExtent):
    freq = float(freq_s)
    
    # Formerly used sample tool, but it was giving erroneous results here

    barrier_pts = intermediates_dir + "topo/barrier_pts.shp"
    arcpy.RasterToPoint_conversion(barrier, barrier_pts)
        
    raster_vec = [eucdir, dem_ft, eucdist_ft]    
    field_vec = ["EUCDIR", "DEM_FT","EUCDISTFT"]
    for i in range(len(raster_vec)):
        target_raster = raster_vec[i]
        target_field = field_vec[i]        
        barrier_pts_new = intermediates_dir + "topo/barrier_pts%s.shp" % i
        
        arcpy.sa.ExtractValuesToPoints(barrier_pts,target_raster, barrier_pts_new, "INTERPOLATE", "VALUE_ONLY")
        
        # Extracted value is in a field called RASTERVALU - it needs to be copied out and this field deleted before the next repeat of the tool.
        arcpy.AddField_management(barrier_pts_new, target_field, "LONG")
        arcpy.CalculateField_management(barrier_pts_new, target_field, '!RASTERVALU!',"PYTHON")
        arcpy.DeleteField_management (barrier_pts_new, "RASTERVALU")
        barrier_pts = barrier_pts_new #update barrier_pts for next iteration
    
    # Assign barrier elevation value by Euclidean direction
    elev_bar_dir = intermediates_dir + "topo/elev_bar_dir"
    elev_barrier = intermediates_dir + "topo/elev_barrier"
    reclass_table = arcpy.sa.ReclassByTable(eucdir, barrier_pts, "EUCDIR", "EUCDIR", "DEM_FT", "NODATA")
    reclass_table.save(elev_bar_dir)
    euc_allocation = arcpy.sa.EucAllocation(elev_bar_dir) #Note: This will interpolate values for the missing directions
    euc_allocation.save(elev_barrier)
    
    # Assign barrier distance value by Euclidean direction
    #NOTE: eucdist_ft was previously defined as the path - this reassigns it locally within the function    
    eucdist_ft_raster = arcpy.sa.Raster(eucdist_ft)
    dist_bar_dir = intermediates_dir + "topo/dist_bar_dir"
    dist_barrier = intermediates_dir + "topo/dist_barrier"
    reclass_table = arcpy.sa.ReclassByTable(eucdir, barrier_pts, "EUCDIR", "EUCDIR", "EUCDISTFT", "NODATA")
    reclass_table.save(dist_bar_dir)
    #EucAllocation interpolates for missing values.
    euc_allocation = arcpy.sa.EucAllocation(dist_bar_dir)
    euc_allocation.save(dist_barrier)
        
    # Make constant raster of sound source elevation
    elev_source = intermediates_dir + "topo/elev_source"
    ccRaster = arcpy.sa.CreateConstantRaster(elev, "FLOAT", CellSize, RasterExtent)
    ccRaster.save(elev_source)
    
    # Calculate slope between source and receiver
    slope = intermediates_dir + "topo/slope"
    dem_ft_raster = arcpy.sa.Raster(dem_ft) # Change from working with a path to working with a raster
    elev_source = arcpy.sa.Raster(intermediates_dir + "topo/elev_source")
    calc_result = (dem_ft_raster - elev_source) / eucdist_ft_raster
    calc_result.save(slope)
    
    # Calculate elevation of source-receiver line under barrier
    elev_sr = intermediates_dir + "topo/elev_sr"
    slope = arcpy.sa.Raster(intermediates_dir + "topo/slope")
    dist_barrier = arcpy.sa.Raster(intermediates_dir + "topo/dist_barrier")
    elev_source = arcpy.sa.Raster(intermediates_dir + "topo/elev_source")
    calc_result = (slope * dist_barrier) + elev_source 
    calc_result.save(elev_sr)
    
    # Calculate barrier height
    h_b = intermediates_dir + "topo/h_b"
    minus_result = arcpy.sa.Minus(elev_barrier, elev_sr) 
    minus_result.save(h_b)
    
    # Reclassify negative barrier height values to zero
    h_b_rc = intermediates_dir + "topo/h_b_rc"
    reclass = arcpy.sa.Reclassify(h_b, "VALUE", "-100000 0 0", "DATA")
    reclass.save(h_b_rc)
    
    # Calculate barrier path distance (BPD)
    bar_pathdist = intermediates_dir + "topo/bar_pathdist"
    h_b_rc = arcpy.sa.Raster(intermediates_dir + "topo/h_b_rc")
    dist_barrier = arcpy.sa.Raster(intermediates_dir + "topo/dist_barrier")
    #eucdist_ft = arcpy.sa.Raster(intermediates_dir + "eucdist_ft") # Defined above
    # Break into smaller steps to avoid errors
    # There appear to be some differences from expectation due to rounding
    # in term1
    term1 = arcpy.sa.SquareRoot(arcpy.sa.Square(h_b_rc) + arcpy.sa.Square(dist_barrier))
    term1.save(intermediates_dir + "topo/term1")
    term2 = arcpy.sa.SquareRoot(arcpy.sa.Square(h_b_rc) + arcpy.sa.Square(eucdist_ft_raster - dist_barrier))
    term2.save(intermediates_dir + "topo/term2")
    term3 = eucdist_ft_raster
    term3.save(intermediates_dir + "topo/term3")
    calc_result = term1 + term2 - term3
    calc_result.save(bar_pathdist)
   
    # Calculate barrier factor (N)
    bar_factor = intermediates_dir + "topo/bar_factor"
    L = ((0.0000000000005*freq**4) - (0.000000001*freq**3) - (0.0000004*freq**2) + (0.0028*freq) - (0.3051))
    times_result = arcpy.sa.Times(bar_pathdist, L)
    times_result.save(bar_factor)
                              
    # The coefficients appear to be from a power regression based on Table 14.
    # Note that I got 13.451 and 0.2427 for my coefficients based on medians.
    # Calculate barrier loss
    OutRaster1 = intermediates_dir + "topo/bar" + freq_s
    bar = (13.573) * (arcpy.sa.Power(bar_factor, 0.2299))                    
    bar.save(OutRaster1)
    
    return(bar)    
    
    
# DELINEATE BARRIER (and calculate locations where ground effects dominate)
def delineate_barrier(sound_src, dem_ft, intermediates_dir):
    
    # Delineate basins
    flow_dir = intermediates_dir + "topo/flow_dir"
    all_basins = intermediates_dir + "topo/all_basins"
    flow_direct = arcpy.sa.FlowDirection(dem_ft)
    flow_direct.save(flow_dir)
    basin = arcpy.sa.Basin(flow_dir)
    basin.save(all_basins)
       
    # Identify sound source basin
    sound_src_basin = intermediates_dir + "topo/sound_src_basin.dbf"
    basin = intermediates_dir + "topo/basin"
    #arcpy.sa.Sample(all_basins, sound_src, sound_src_basin, "BILINEAR")
    # Changed to nearest, bilinear interpolation here makes no sense!
    # You are selecting a basin, and need it to remain as an integer value!
    arcpy.sa.Sample(all_basins, sound_src, sound_src_basin, "NEAREST") 

    # Use search cursor to extract the basin corresponding to the source location    
    # Changed to use with syntax to avoid schema locks
    with arcpy.da.SearchCursor(sound_src_basin, "all_basins") as cur:
        for row in cur:
            basin_num = row[0]
    
    extract = arcpy.sa.ExtractByAttributes(all_basins, "VALUE = " + str(basin_num))
    extract.save(basin)
        
    # Delineate barrier (ie, ridgeline) around sound source basin
    basin_exp = intermediates_dir + "topo/basin_exp"
    basin_shr = intermediates_dir + "topo/basin_shr"
    expand = arcpy.sa.Expand(basin, "1", "'" + str(basin_num) + "'")
    expand.save(basin_exp)
    shrink = arcpy.sa.Shrink(basin_exp, "1", "'" + str(basin_num) + "'")
    shrink.save(basin_shr)
    
    basin_shr_rc = intermediates_dir + "topo/basin_shr_rc"
    barrier_raster = intermediates_dir + "topo/barrier"
    basin_plus1 = basin_num + 1
    Reclass = str(basin_num) + " NoData; NoData " + str(basin_plus1)
    reclass_result = arcpy.sa.Reclassify(basin_shr, "VALUE", Reclass)
    reclass_result.save(basin_shr_rc)
    barrier = arcpy.sa.Minus(basin_shr_rc, basin_exp)
    barrier.save(barrier_raster)
    
    # Define areas where ground effects dominate
    ground_raster = intermediates_dir + "topo/ground"
    Reclass = str(basin_num) + " 3; NoData 0"
    ground = arcpy.sa.Reclassify(basin_exp, "VALUE", Reclass)
    ground.save(ground_raster)
    
    return(barrier, ground)    

# IDENTIFY AREAS WHERE GROUND, ATMOSPHERIC, AND BARRIER EFFECTS DOMINATE
def calculate_topozones(Sound_Source, dem_ft, ground, intermediates_dir, sound_src):
    # Environment settings
    #arcpy.env.extent = RasterExtent
    #arcpy.env.cellSize = CellSize
    #arcpy.env.snapRaster = arcpy.Raster(dem_ft)
    
    # This section required patching due to problems with the viewshed tool
    # NOTE: THERE IS A VIEWSHED2 TOOL - BUT I WAS HAVING TROUBLE WITH IT IN A
    # TEST LANDSCAPE and it is only compatible with ArcGIS 10.3+

    # Alternate approach using observer points code
    # potential for efficiency gain as multiple points can be calculated at once.
    viewshed = intermediates_dir + "topo/viewshed"
    agl_raster = intermediates_dir + "topo/agl"    
    dem_ft_copy = arcpy.Raster(dem_ft)
    dem_ft_copy.save(intermediates_dir + "topo/dem_ft_copy")
    view = arcpy.sa.ObserverPoints(dem_ft_copy, Sound_Source, "0.3048", "FLAT_EARTH", "0.13", agl_raster)    
    view.save(viewshed)
      
    '''   
    # Calculate viewshed for source point
    agl_raster = intermediates_dir + "topo/agl"    
    viewshed = intermediates_dir + "topo/viewshed"
    viewshed_result = arcpy.sa.Viewshed(viewshed, Sound_Source, "0.3048","FLAT_EARTH",0.13,agl_raster)
    viewshed_result.save(viewshed)
    '''    
    # ArcGIS 9.3 compatibility version
    #viewshed = "C:/SPreAD-GIS/backup/Arc9_3/intermediate/topo/viewshed"
    
    # NOTE: The ground input is not explicitly used, but see paths below
    # it is implicilty used & must be first created (hence passing it as a
    # function argument)
    # Define areas where ground effects and atmospheric effects dominate
    ground_atmos = intermediates_dir + "topo/ground_atmos"
    highest_position = arcpy.sa.CellStatistics([intermediates_dir + "topo/ground", viewshed], "MAXIMUM")    
    highest_position.save(ground_atmos)
    
    # Define areas where ground, atmospheric, and barrier effects dominate
    topo_zones_raster = intermediates_dir + "topo/topo_zones"
    topo_zones = arcpy.sa.Reclassify(ground_atmos, "VALUE", "0 2; NoData 2")
    topo_zones.save(topo_zones_raster)
    
    return(topo_zones)

# CALCULATE BARRIER EFFECTS IN A WAY THAT AVOIDS IMPLEMENTATION ARTIFACTS
# This implementation will almost certainly be slower than the previous one
# Note: Much of this code was selectively copied and pasted from nmsimhlpr.py
# as the basic concept is the same
def barrier_effects_v2(barrier_path_distance, freq_s, intermediates_dir):
    freq = float(freq_s)

    # Calculate barrier factor (N)
    bar_factor = intermediates_dir + "topo/bar_factor"
    L = ((0.0000000000005*freq**4) - (0.000000001*freq**3) - (0.0000004*freq**2) + (0.0028*freq) - (0.3051))
    times_result = arcpy.sa.Times(barrier_path_distance, L)
    times_result.save(bar_factor)
                              
    # The coefficients appear to be from a power regression based on Table 14.
    # Note that I got 13.451 and 0.2427 for my coefficients based on medians.
    # Calculate barrier loss
    bar_file = intermediates_dir + "topo/bar" + freq_s
    bar = (13.573) * (arcpy.sa.Power(bar_factor, 0.2299))                    
    bar.save(bar_file)

    return bar_file

# Calculate height of largest single barrier for each cell
# This also calculates the maximum vegetation loss (same function because it
# means I only have to loop through the landscape once)
#**# Think about if and how to break landscape into chunks to alleviate memory problems
def calculate_barrier_path_distance_and_vegmax(Sound_Source, dem_ft, landcover, eucdist_ft, source_offset, receiver_offset, idir):
    # For testing purposes
    #dem_ft = r'C:\smt\csu\simple_landscapes\Test4_elevation2.tif'

    m2ft = 3.28084

    # Describe dem_ft to get basic properties
    properties = arcpy.Describe(dem_ft)
    xmin = properties.extent.XMin
    ymin = properties.extent.YMin
    nrow = properties.height
    ncol = properties.width
    cell_size = properties.meancellheight
    # Requires square cells, check that this is true
    if round(properties.meancellheight, 8) != round(properties.meancellwidth,8):
        raise ValueError("Cell height %s must equal cell width %s" % (properties.meancellheight, properties.meancellwidth))

    # Read the elevation raster into memory
    numpy_dem = arcpy.RasterToNumPyArray(dem_ft)

    # Read landcover raster into memory
    numpy_lc = arcpy.RasterToNumPyArray(landcover)

    #**# For efficiency, you might just want to recode the values in numpy_lc to use SPREADTYPE here
    # But only re-code if this appears to be rate limiting, no sense optimizing if it's not a problem.
    # Get attribute table for numpy_lc
    lc_attribute_table = get_attribute_table(landcover)

    # Create a numpy array of the same size as numpy_dem
    dem_dimensions = (nrow, ncol)
    barheight_array = numpy.zeros(dem_dimensions) # starting with 0's
    bardist_array = numpy.zeros(dem_dimensions)
    vegmax_array = numpy.zeros(dem_dimensions)

    #**# Think about adding bilinear interpolation for elevation extraction
    # Find source coordinates on landscape
    src_properties = arcpy.Describe(Sound_Source)
    src_xmin = src_properties.extent.XMin
    src_ymin = src_properties.extent.YMin
    # src_xmin and src_ymin will be greater than xmin and ymin respectively,
    # otherwise the source is not in the landscape!
    # Floor is used because fractions are still within the same cell, regardless
    # of how close they are to the cell border
    source_col = (src_xmin - xmin) * cell_size**-1
    source_col = int(math.floor(source_col))
    source_row = (src_ymin - ymin) * cell_size**-1
    source_row = int(nrow - math.floor(source_row) - 1) # Need to subtract - min gives distance from bottom of raster, but we want distance from the top!

    source_cell_elevation = numpy_dem[source_row, source_col]
    source_elevation = source_cell_elevation + source_offset * m2ft
    xyzsrc_cell = [source_col, source_row, source_cell_elevation]    
    
    # For each cell in the landscape:    
    for m in xrange(ncol):
        for n in xrange(nrow):
            # Calculate height of the relevant barrier
            # Note: SPreAD only looks at barrier height above source height.
            # Here we took a more sophisticated approach and looked at height above the
            # line between the source and the receiver

            # Get cell elevation & put cell information into a single variable
            z = numpy_dem[n, m]
            xyzrec_cell = [m,n,z]

            # Get a terrain cut between the source and the cell
            # BASED ON GET_TRANSECT_V2 function in nmsimhlpr.py
            max_dist, dist_vec, terrain_cut, veg_cut = get_terrain_cut(xyzsrc_cell, xyzrec_cell, cell_size, cell_size, numpy_dem, numpy_lc, lc_attribute_table)

            # Get distance between cell and source
            #distance_cell = math.sqrt((n - source_row)**2 + (m - source_col)**2)
            
            # Convert from m to feet (unfortunately)
            #distance = distance_cell * cell_size * m2ft
            distance = max_dist * m2ft
            
            # At cell of origin, distance will be 0. This will crash the slope calculation
            # Make distance a very small number instead (and note that
            # receiver_elevation - source_elevation for the point of origin is
            # expected to be 0, so it won't matter what distance is specified,
            # because it is multiplied by zero)
            if distance == 0:
                distance = 0.01

            # Get elevation of the cell
            receiver_elevation = z + receiver_offset * m2ft
    
            # Calculate the slope between the source and that cell
            slope = (receiver_elevation - source_elevation) * distance**-1
            
            # Calculate the height above the slope for the terrain cut between the source and that cell
            max_height = 0 # Negative heights do not obstruct sound propagation
            bar_dist = distance
            for i in xrange(len(dist_vec)):
                slope_elevation = slope * dist_vec[i] * m2ft + source_elevation # Slope * distance gives the increase or decrease in height between source and receiver
                height_above_slope = terrain_cut[i] - slope_elevation
                
                # Calculate the maximum height above the slope along the terrain cut
                if height_above_slope > max_height:
                    max_height = height_above_slope
                    bar_dist = dist_vec[i] * m2ft
        
            # Assign that height to the cell
            barheight_array[n,m] = max_height
            bardist_array[n,m] = bar_dist

            ## Calculate maximum vegetation loss for this cell
            max_veg_value = calc_vegmax(veg_cut, cell_size, distance, m2ft)
            vegmax_array[n,m] = max_veg_value

    # Set up for writing to file
    lowerleft = arcpy.Point(xmin, ymin)
    no_data_value = -1

    # Write barrier heights to file
    barrier_height_raster_file = idir + "barrier_heights.tif"
    barrier_height_raster = arcpy.NumPyArrayToRaster(barheight_array, lowerleft, cell_size, cell_size, no_data_value)
    barrier_height_raster.save(barrier_height_raster_file)

    # Write barrier distances to file
    barrier_distance_raster_file = idir + "barrier_distances.tif"
    barrier_distance_raster = arcpy.NumPyArrayToRaster(bardist_array, lowerleft, cell_size, cell_size, no_data_value)
    barrier_distance_raster.save(barrier_distance_raster_file)

    # Define projections for barrier heights & distances
    proj = properties.SpatialReference
    arcpy.DefineProjection_management(barrier_height_raster_file, proj)
    arcpy.DefineProjection_management(barrier_distance_raster_file, proj)

    # Calculate Barrier Path Distance (ft)
    # BPD = sqrt(h**2 + R**2) + sqrt(h**2 + (X - R)**2) - X
    # h = barrier height above source
    # X distance between source and cell
    # R distance between source and barrier
    eucdist_ft_raster = arcpy.Raster(eucdist_ft)
    BPD_pre_file = idir + "bpd_pre.tif"
    term1_file = idir + "term1.tif"
    term2_file = idir + "term2.tif"
    term1 = (barrier_height_raster**2 + barrier_distance_raster**2)**(0.5) # 0.5 = square root
    term1.save(term1_file)
    term2 = (barrier_height_raster**2 + (eucdist_ft_raster - barrier_distance_raster)**2)**(0.5)
    term2.save(term2_file)
    BPD = term1 + term2 - eucdist_ft
    BPD.save(BPD_pre_file)

    #**# NOTE: For barrier heights of 0, BPD should calculate out to be equal to 0
    # and therefore should come out to have no barrier effect.
    # However, the barrier distance and euclidean distance calculations have small errors
    # but these lead to large (4 - 7 dB) barrier effects near the origin.
    # use a pick step to set barrier effects in locations without barriers to 0.
    BPD_file = idir + "bpd.tif"
    BPD2 = arcpy.sa.Con((barrier_height_raster == 0), barrier_height_raster, BPD)
    BPD2.save(BPD_file)


    # Converte veg_max to a raster
    veg_max_file = idir + "vegmax.tif"
    veg_max_raster = arcpy.NumPyArrayToRaster(vegmax_array, lowerleft, cell_size, cell_size, no_data_value)
    veg_max_raster.save(veg_max_file)
    arcpy.DefineProjection_management(veg_max_raster, proj)

    return BPD_file, veg_max_file

# Calculate maximum vegetation loss allowable
def calc_vegmax(veg_cut, cellsize, distance, m2ft):
    # Based on Table 8 of Harrison et al. 1980
    # The grassland frequency-specific results are applied in the check_veg function
    # in the main spreadgishlpr.py workflow

    # According to Harrision et al. 1980, at distances < 75 ft, sound loss = 0
    #total_distance_ft = len(veg_cut) * cellsize * 3.28084
    if distance < 75:
        max_veg_loss = 0
    
    # Otherwise, compute the vegetation loss
    else:
        npts = float(len(veg_cut)) # float is needed to get floating point division
        #print veg_cut
        # Get number of m of conifers
        # filter chosen based on this thread: http://stackoverflow.com/questions/12845288/grep-on-elements-of-a-list
        # [n for n in veg_cut if "HEB" in n] would also work, but makes even less sense to me.
        con = filter(lambda x: "CON" in x, veg_cut)
        
        # Note this is in m here
        if len(con) == 0:
            max_con_loss = 0
        else:
            distance_con = (len(con) / npts) * (distance / m2ft)
            max_con_loss = 5.2504 * math.log(distance_con) - 9.8094 # R2 = 0.99
            if max_con_loss < 0:
                max_con_loss = 0
        
        # Get number of m of hardwoods
        hwd = filter(lambda x: "HWD" in x, veg_cut)
        if len(hwd) == 0:
            max_hwd_loss = 0
        else:            
            dist_hwd = (len(hwd) / npts) * (distance / m2ft)
            max_hwd_loss = 6.6224 * math.log(dist_hwd) - 16.762 # R2 = 0.99
            if max_hwd_loss < 0:
                max_hwd_loss = 0

        # In Table 8, vegeation loss is frequency-specific.
        # The frequency dependence is ignored here, and this may overestimate
        # loss compared to SPreAD model at frequencies >800 Hz
        # (but not by more than 4 dB). The potential for overestimation is
        # greater for higher frequencies. The previous version of SPreAD-GIS
        # also disregarded frequency information in calculating the vegetation
        # loss.

        # Get number of m of grassland or shrub
        heb = filter(lambda x: "HEB" in x, veg_cut) + filter(lambda x: "SHB" in x, veg_cut)
        max_heb_loss = 0        
        if len(heb) > 0:
            max_heb_loss = 4 
        #dist_heb = len(heb) * cellsize
        

        # Add sources of vegetation loss
        max_veg_loss = max_con_loss + max_hwd_loss + max_heb_loss

        # Cap total vegetation loss at 14 dB
        if max_veg_loss > 14:
            max_veg_loss = 14

    return max_veg_loss


# Check if barrier plus wind exceeds 25 dB, if so, truncate to 25 dB per original SPreAD model
def check_barrier_wind(bar, wind_raster, freq_s, intermediates_dir):

        # Add barrier & wind loss rasters together
        barwind_file = intermediates_dir + "topo/barwind" + freq_s + ".tif"
        barwind_prelim_file = intermediates_dir + "topo/barwind_prelim_" + freq_s + ".tif"
        bar_raster = arcpy.Raster(bar)
        #wind_raster = wind # wind is already a raster
        barwind_prelim = bar_raster + wind_raster
        barwind_prelim.save(barwind_prelim_file)

        # Check if barrier plus wind exceeds 25 dB, if so, cap to 25 dB
        barwind = arcpy.sa.Con((barwind_prelim > 25), 25, barwind_prelim)
        barwind.save(barwind_file)

        return barwind_file
    
# CALCULATE SUMMARY NOISE PROPAGATION PATTERNS #
def compute_noise_propagation(freq_s, source, ssl, aal, veg, wind, bar, topo_zones, results_dir, prelim_dir, sound_src, keep_intermediates, intermediates_dir):
    
    # Subtract spherical spreading loss from source sound level
    ssl_loss_file = prelim_dir + "sslloss" + freq_s
    ssl_loss = arcpy.sa.Minus(source, ssl)
    ssl_loss.save(ssl_loss_file)
    
    # Calculate cumulative spherical spreading and atmospheric absorption loss
    sslaal_file = prelim_dir + "sslaal" + freq_s
    sslaal = arcpy.sa.Minus(ssl_loss, aal)
    sslaal.save(sslaal_file)
       
    # Calculate cumulative spherical spreading, atmospheric absorption, and vegetation loss
    salveg_file = prelim_dir + "salveg" + freq_s
    salveg = arcpy.sa.Minus(sslaal, veg)
    salveg.save(salveg_file)
    
    # Calculate cumulative spherical spreading, atmospheric absorption, vegetation, and wind loss
    salvegwin_file = prelim_dir + "salvegwin" + freq_s
    salvegwin = arcpy.sa.Minus(salveg, wind)
    salvegwin.save(salvegwin_file)
                   
    # Calculate cumulative spherical spreading, atmospheric absorption, vegetation, wind, and barrier loss
    salvgwnbr_file = prelim_dir + "salvgwnbr" + freq_s
    salvgwnbr = arcpy.sa.Minus(salvegwin, bar)
    salvgwnbr.save(salvgwnbr_file)
    
    # Pick noise propagation values for areas where ground, atmospheric,
    # or barrier effects dominate
    combined_file = prelim_dir + "combined" + freq_s
    combined = arcpy.sa.Pick(topo_zones, [sslaal, salvgwnbr, salvegwin])
    combined.save(combined_file)
                            
    # Smooth noise propagation patterns
    #smoothed_file = prelim_dir + "smooth" + freq_s
    smoothed_file = results_dir + "smth" + freq_s + ".tif"
    smoothed = arcpy.sa.FocalStatistics(combined, "Rectangle 3 3 CELL", "MEAN", "DATA")
    #smoothed.save(smoothed_file)
    smoothed.save(smoothed_file)
     
    # Patch to prevent smoothing at cell of origin
    pr_file = results_dir + "pr" + freq_s + ".tif"
    pr = arcpy.sa.Con((sound_src > 0), combined, smoothed)
    pr.save(pr_file)
    
    # Convert intermediates to show what was actually used in the final result
    # to improve transparency in displaying the intermediates 
    if keep_intermediates == 1 or keep_intermediates == "YES":

        # include paths to the intermediate outputs
        veg = intermediates_dir + "veg/veg%s" % freq_s
        wind = intermediates_dir + "wind/win%s" % freq_s
        bar = intermediates_dir + "topo/bar%s" % freq_s

        # Define new variables
        veg_final_file = intermediates_dir + "veg/vegfinal_%s.tif" % freq_s
        wind_final_file = intermediates_dir + "wind/windfinal_%s.tif" % freq_s
        bar_final_file = intermediates_dir + "topo/barfinal_%s.tif" % freq_s

        # Areas where vegetation effects are excluded
        veg_final = arcpy.sa.Pick(topo_zones, [0, veg, veg])
        veg_final.save(veg_final_file)
        
        #arcpy.sa.Pick(topo_zones, [sslaal, salvgwnbr, salvegwin]
        # Areas where wind effects are excluded
        wind_final = arcpy.sa.Pick(topo_zones, [0, wind, wind])
        wind_final.save(wind_final_file)

        # Areas where barrier effects are excluded
        bar_final = arcpy.sa.Pick(topo_zones, [0, bar, 0])        
        bar_final.save(bar_final_file)

    #**# Disabled because this lead to non-0 values in the overall summary
    # Remove propagation values that are less than zero
    #pr_file = results_dir + "pr" + freq_s
    #pr = arcpy.sa.Con((smoothed > 0), smoothed, 0)
    #pr.save(pr_file)

# CALCULATE SUMMARY NOISE PROPAGATION PATTERNS #
# Version to correspond to updated barrier calculation
def compute_noise_propagation_v2(freq_s, source, ssl, aal, veg, barwind, results_dir, prelim_dir, sound_src, keep_intermediates, intermediates_dir):

    # NOTE: V2 uses a combined wind/barrier raster, an improved barrier calculation
    # that does not require a distinction between "in basin" and "out of basin"
    # selection

    # The arcpy.sa.Pick step has been completely removed. One potential downside
    # is that intervening vegetation will now influence the result even if
    # there is a direct line of sight between the source and a receiver at a
    # high elevation location (e.g., valley vegetation will be included in the
    # calculation, even if the sound is completely bypassing the valley).
    # To minimize problems associated with this, the caps to vegetation loss
    # have been re-instated in modified form to account for the mixed landscape
    # (the original caps assumed a homogeneous landscape)
    
    # Subtract spherical spreading loss from source sound level
    ssl_loss_file = prelim_dir + "sslloss" + freq_s
    ssl_loss = arcpy.sa.Minus(source, ssl)
    ssl_loss.save(ssl_loss_file)
    
    # Calculate cumulative spherical spreading and atmospheric absorption loss
    sslaal_file = prelim_dir + "sslaal" + freq_s
    sslaal = arcpy.sa.Minus(ssl_loss, aal)
    sslaal.save(sslaal_file)
       
    # Calculate cumulative spherical spreading, atmospheric absorption, and vegetation loss
    salveg_file = prelim_dir + "salveg" + freq_s
    salveg = arcpy.sa.Minus(sslaal, veg)
    salveg.save(salveg_file)
    
    # Calculate cumulative spherical spreading, atmospheric absorption, vegetation, wind, and barrier loss
    salvgwnbr_file = prelim_dir + "salvgwnbr" + freq_s
    salvgwnbr = arcpy.sa.Minus(salveg, barwind)
    salvgwnbr.save(salvgwnbr_file)
                                
    # Smooth noise propagation patterns
    #smoothed_file = prelim_dir + "smooth" + freq_s
    smoothed_file = results_dir + "smth" + freq_s + ".tif"
    smoothed = arcpy.sa.FocalStatistics(salvgwnbr, "Rectangle 3 3 CELL", "MEAN", "DATA")
    #smoothed.save(smoothed_file)
    smoothed.save(smoothed_file)
     
    # Patch to prevent smoothing at cell of origin
    pr_file = results_dir + "pr" + freq_s + ".tif"
    pr = arcpy.sa.Con((sound_src > 0), salvgwnbr, smoothed)
    pr.save(pr_file)
        

# Custom version for SPreAD-GIS
# Create a timer function for timing code
def sg_get_time(do_timing, my_times, my_time_labels, label):
    import time    
    if do_timing == 1:
        this_time = time.time()
        my_times.append(this_time)
        my_time_labels.append(label)

    return(my_times, my_time_labels)

def print_message(message, print_indicator):
    if print_indicator == 1:
        arcpy.AddMessage(message)
        
# Calculate euclidean distance and direction for SPreAD-GIS calculations
# Moved from spreadgishlpr.py on 11/4/2016 to try to prevent schema locks
# from crashing the code, then moved back
def euclidean_dist_dir(Sound_Source, sound_src, eucdist, eucdist_ft, eucdir, CellSize):

    # Check if "ONE" field already exists, if so, delete it
    fields = arcpy.ListFields(Sound_Source)
    for field in fields:
        if field.name == "ONE":
            arcpy.DeleteField_management(Sound_Source, "ONE")

    # Convert sound source to raster dataset
    arcpy.AddField_management(Sound_Source, "ONE", "SHORT")
    arcpy.CalculateField_management(Sound_Source, "ONE", "1", "PYTHON")
    #tempEnvironment0 = arcpy.Extent()
    #arcpy.env.extent = Model_Extent
    arcpy.PointToRaster_conversion(Sound_Source, "ONE", sound_src, "", "", CellSize)
    #arcpy.env.extent = tempEnvironment0
    
    # Calculate Euclidean allocation, distance, and direction around sound source
    eucDirection = arcpy.sa.EucDirection(sound_src, "", CellSize, eucdist)
    eucDirection.save(eucdir)
    
    # Convert Euclidean distance from meters to feet
    times_result = arcpy.sa.Times(eucdist, 3.2808)
    times_result.save(eucdist_ft)

# CALCULATE MEAN ELEVATION
# Function to calculate the mean elevation for the landscape to provide a needed
# input to SPreAD-GIS
def calculate_mean_elevation(dem_ft, sound_src, intermediates_dir):
    # Sample elevation values under sound source
    sound_src_elev_dbf = intermediates_dir + "sound_src_elev.dbf"
    #sound_src_elev_mean_dbf = intermediates_dir + "aal/sound_src_elev_mean.dbf"
    arcpy.sa.Sample(dem_ft, sound_src, sound_src_elev_dbf, "BILINEAR") #**# This is not providing a bilinear interpolation in 10.3.1
    
    # Get mean elevation value of sound source
    # This step does not appear to be useful - you have the same answer as the previous table!
    #arcpy.Statistics_analysis(sound_src_elev_dbf, sound_src_elev_mean_dbf, [["dem_ft", "MEAN"]])
        
    # Get elevation value for sound source from summary table (only has one row)
    with arcpy.da.SearchCursor(sound_src_elev_dbf, "dem_ft") as cur:
        for row in cur:
            elev = row[0]
            
    # OLD VERSION OF CURSOR
    #cur = arcpy.SearchCursor(sound_src_elev_dbf)
    #row = cur.next()
    #
    #try:
    #    while row:
    #        #elev = row.getValue("MEAN_dem_f")
    #        elev = row.getValue("dem_ft")
    #        row = cur.next()
    #except:
    #    del cur, row
    #    raise ValueError(sys.exc_info()[0:2])
    return(elev)


# The following were copied and modified from nmsimhlpr.py on 11/10/16:
#    get_transect_v2 (= get_terain_cut)
#    find_distances
#    find_coords
#    extract_values


# RETURN ELEVATION AND VEGETATION TRANSECTS
def get_terrain_cut(xyzsrc_cell, xyzrec_cell, cell_x, cell_y, numpy_dem, numpy_impedance, lc_attribute_table):

    # calculate max_dist as distance between xyzsrc & xyzrec (without respect to elevation)
    xsource = xyzsrc_cell[0]
    ysource = xyzsrc_cell[1]
    xreceiver = xyzrec_cell[0]
    yreceiver = xyzrec_cell[1]
    xdiff = (xreceiver - xsource) * cell_x # Convert back to m
    ydiff = (yreceiver - ysource) * cell_y # convert back to m
    max_dist = math.sqrt( (xdiff)**2 + (ydiff)**2 )

    # Identify transect points
    npts = 100  # Bruce Ikelheimer said that ONECUT algorithm is ALWAYS given 100 points.
    dist_vec = find_distances(max_dist, npts)

    # Get coordiantes for dist_vec locations
    x_coords = find_coords(xsource, xreceiver, npts)
    y_coords = find_coords(ysource, yreceiver, npts)
    
    
    # use coordinates to extract elevation
    hgt_vec = extract_values(numpy_dem, x_coords, y_coords)

    # use coordiantes to extract landcover
    lc_vec = extract_values(numpy_impedance, x_coords, y_coords, "NEAREST")

    # Reclassify landcover values from raster values to SPREADTYPE values
    lc_vec = reclass_to_spreadtype(lc_vec, lc_attribute_table)

    # Check that everything is the correct length
    if len(dist_vec) != npts or len(hgt_vec) != npts or len(lc_vec) != npts:
        raise ValueError("npts (%s), dist_vec (%s), hgt_vec (%s), and lc_vec (%s) all must be the same length." % (npts, dist_vec, hgt_vec, lc_vec))
    
    return(max_dist, dist_vec, hgt_vec, lc_vec)

# This function is derived from the function of the same name in nmsimhlpr.py
def find_distances(max_dist, npts, n_digits = 8): #offset = 0, 
    increment = max_dist * (npts - 1)**-1 # - 1 accounts for the inclusion of the end points, **-1 does division
    increment = round(increment, n_digits)

    # If the increment is 0, then just create a vector of 0's
    #if round(increment, n_digits) == 0:
    if increment == 0:
        dist_vec = [0] * npts

    else:
        # create dist_vec
        dist_vec = seq(0,max_dist, increment)
        #**# Below is no longer necessary with change in where rounding occurs
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
        
        #**# NOTE: THERE WAS FORMERLY A BUG IN THIS. NOT SURE IF IT IS STILL THERE. THIS OPTION IS NOT CURRENTLY BEING USED
        if extract_type == "NEAREST":
            #raise ValueError("There is an undiagnosed problem where the nearest neighbor calculation fails (e.g., from Source 1,8 to Receiver 0,8 in the designed landscape)")
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

# COPIED FROM NMSIMHLPR.PY
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
    # Added rounding to 8 digits, because I'm getting an error about overflow in short scalars! Which shouldn't be happening, but it is.

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


# Function to extract the corresponding SPREADTYPE for the values in a raster
def get_attribute_table(landcover):
    # initialize a dictionary    
    lc_attribute_table = {}
    
    with arcpy.da.SearchCursor(landcover, ["VALUE","SPREADTYPE"]) as cur:
        for row in cur:
            lc_attribute_table[row[0]] = row[1]
            
    return lc_attribute_table

# Reclassify landcover values from raster values to SPREADTYPE values (note: this is different than adding the SPREADTYPE field
# this is an internal reclassification using the existing SPREADTYPE field, because the text values aren't stored as the raster
# values, so when the raster is converted to a NUMPY array that information is lost)
def reclass_to_spreadtype(lc_vec, lc_attribute_table):
    
    # If this is a rate limiting step, this could probably be made to run faster    
    new_lc_vec = []
    # Go through each vector element
    for lc in lc_vec:

        # Look up value in attribute table        
        new_lc = lc_attribute_table[lc]
                
        # return the match
        new_lc_vec.append(new_lc)
    
    return new_lc_vec
        

# END OF FILE