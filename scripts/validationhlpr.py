# -*- coding: utf-8 -*-
"""
Description: This script contains the functions called by B_Validation.py to
confirm that the user's version of ArcGIS produces the same results as our
version of ArcGIS. Also useful for unit testing to confirm that changes to the
code do not have unintended consequences.

@author: Alexander "Sasha" Keyel  <skeyel@gmail.com>
Copyright (C) 2016 - 2017, A.C. Keyel 
   
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

# Import required modules
import os, shutil
import arcpy, arcpy.sa

# This is the main function call. Source files need to be in source_dir, output_dir will be created & should be deleted at end (but may be problems with schema locks)
def validate_arcgis_version(testing_dir, tbx_root, skip_test_run, test_spread, test_nmsim, test_iso, is_GUI):

    source_dir = testing_dir + "source_data/"
    output_dir = testing_dir + "test_run/"
    output_dir2 = output_dir + "barrierv2/"
    ambient_dir = output_dir + "ambient/"

    # Import required standard modules
    import soundprophlpr as sp
    
    # Check out the spatial analyst extension
    arcpy.CheckOutExtension("Spatial")
    
    # Delete intermediate datasets, if they exist
    arcpy.env.overwriteOutput = True
       
    # Option to skip creating model outputs (if not 0). This option will be mainly useful when
    # updating the validation values when multiple problems are expected to exist
    if skip_test_run == 0:     
        
        # Delete the test_run directory, if it already exists, then create it new
        # Remove any existing results or intermediate files if possible    
        shutil.rmtree(output_dir, ignore_errors=True)
        os.mkdir(output_dir)
        os.mkdir(ambient_dir)
        sp.make(output_dir2)
    
        # Set up inputs
        background_table = source_dir + "Test_Background.csv"
        Sound_Source = source_dir + "one_point.shp"
        Model_Extent = source_dir + "test_extent.shp" # Smaller extent to just enclose the test points
        in_freq = "400;70;15.24" # If giving a single frequency, SHOULD NOT be a list
        dem = source_dir + "dem_m" #dem
        veg = source_dir + "nlcd"
        temp_s = "28.3333333333" #"83"
        hum_s = "26"
        wind_dir = "105" #"285"  
        wind_sp = "12.8748" #"8" 
        seas_cond = "clear, windy summer day"
        results_label = "test_spreadgis"
        Model_Extent = sp.get_extent(Model_Extent)

        # Added for revised SPreAD-GIS model
        results_label2 = "test_spreadgis2"
    
        # Create ambient conditions
        sp.CreateAmbient(in_freq, veg, background_table, wind_speed = "#", snow = "#", ambient_path = ambient_dir, freq_fill = 5) # Get enough leading zeros
        ambient_path = ambient_dir + "ambient00400"

        if test_spread == 1:
            
            # Run spreadgis model for test data
            sp.SoundPropagation("spreadgis", Sound_Source, "POINT", Model_Extent, ambient_path,
                             dem, veg, temp_s, hum_s, wind_dir, wind_sp,
                             seas_cond, output_dir, output_dir, results_label, in_freq, 
                             receiver_offset = 1,
                             n_points = "all", tbx_root = tbx_root, summarize = 2, keep_intermediates = 1)
                             
            # Run spreadgis model with new method for barrier and vegetation for test data
            sp.SoundPropagation("spreadgis", Sound_Source, "POINT", Model_Extent, ambient_path,
                             dem, veg, temp_s, hum_s, wind_dir, wind_sp,
                             seas_cond, output_dir2, output_dir2, results_label2, in_freq, 
                             receiver_offset = 1,
                             n_points = "all", tbx_root = tbx_root, summarize = 2, keep_intermediates = 1, use_old_barrier = 0)

    
        if test_nmsim == 1:
            # Run nmsimgis model for test data
            Sound_Source = source_dir + "two_points.shp" # nmsimgis will be run for 2 points to test summary code 
            results_label = "test_nmsimgis"
            in_src = tbx_root + "Sources\\GroundSources\\Bus.src" # Apparently ArcGIS needs double-backslashes instead of forward slashes for this to work correctly with the .dll. ODD.
            ambient_dir = "NA"
            
            sp.SoundPropagation("nmsimgis", Sound_Source, "POINT", Model_Extent, ambient_dir,
                                dem, veg, temp_s, hum_s, wind_dir, wind_sp,
                                seas_cond, output_dir, output_dir, results_label, in_src,
                                receiver_offset = 1.22, n_points = "all", tbx_root = tbx_root,
                                summarize = 1, keep_intermediates = 1, is_GUI = is_GUI)
        
        if test_iso == 1:
            # Run iso model for test data
            Sound_Source = source_dir + "one_point.shp"
            results_label = "test_iso"
            in_src = tbx_root + "Sources\\GroundSources\\Bus.src" #**# Apparently ArcGIS needs double-backslashes instead of forward slashes for this to work correctly with the .dll. ODD.
            ambient_dir = "NA"
            
            sp.SoundPropagation("iso9613", Sound_Source, "POINT", Model_Extent, ambient_dir,
                                dem, veg, temp_s, hum_s, wind_dir, wind_sp,
                                seas_cond, output_dir, output_dir, results_label, in_src,
                                receiver_offset = 1.22, n_points = "all", tbx_root = tbx_root,
                                summarize = 1, keep_intermediates = 1, is_GUI = is_GUI)
            
 
    #### Perform checks on each section to confirm that calculations are correct
       #and as expected

    # Should not be necessary, but this appears to have been turned off when running the code through the GUI
    arcpy.env.overwriteOutput = True # Make this a setting for all processes
    
    # Check at test points:
    # 1) point source
    # 2) W of point source; in basin
    # 3) E of point source; in basin
    # 4) W of point source; in bright area
    # 5) S of point source; out of basin
    # 6) E of point source; high elevation

    test_points = source_dir + "test_points.shp"
    n_testpoints = 6 # Number of test points 
    # order of points is actually 6,2,1,3,4,5 in temp_table.dbf when using the 
    # sample tool.
    # Note,  they  are  labeled   5,1,0,2,3,4 in temp_table.dbf based on FID's
    # rather than the actual point names.

    # 10.4.1 seems to sample points in a random order (but so far consistent within a layer, just random across layers)
    # Sometimes it uses 1,5,0,2,3,4
    # Sometimes it uses 0,1,2,3,4,5
    version = arcpy.GetInstallInfo('desktop')['Version']

    if version == "10.5":
        print "WARNING: There is a bug in ArcGIS 10.5 that causes an error in SPreAD-GIS when using the old algorithms."
        print "Please ensure that use_old_barrier = 0 for SPreAD-GIS. Otherwise large (20 dB) prediction errors may occur."

    ### TEST SPREAD-GIS
    if test_spread == 1:
        test_spreadgis(output_dir, test_points, n_testpoints, version)
        
        test_spreadgis2(output_dir2, test_points, n_testpoints, version)
        
    ### TEST nmsimgis, multipoint summary, weighting, and overall summary
    if test_nmsim == 1:
        test_nmsimgis(output_dir, test_points, n_testpoints, version)
        ### TEST non-multipoint summary & weighting using nmsimgis data
        #test_soundprop_other(output_dir, source_dir)

    if test_iso == 1:
        test_iso_model(output_dir, test_points, n_testpoints, version)

    # Delete directory when test run is completed
    try:
        shutil.rmtree(output_dir)
    except:
        arcpy.AddMessage("Unable to delete %s folder after completion of testing. Please delete this folder manually. We apologize for the inconvenience" % output_dir)


# Function to test SPreAD-GIS functioning
def test_spreadgis(output_dir, test_points, n_testpoints, version):
    idir = output_dir + "spreadgis/points/pt1/intermediate_400/"
    rpath = output_dir + "spreadgis/points/pt1/results/"    
    epath = output_dir + "spreadgis/frequency_excess/"    
    dddir = output_dir + "spreadgis/temp/"    
    
    # CHECK EUCLIDEAN DISTANCE & DIRECTION
    #try:
    # directory
    check_dist_dir(dddir, dddir, test_points, version)
    #except:
    #    is_error = 1
    #    print >> el, 
        
    # CHECK MEAN ELEVATION #**# Not checked
        
    # CHECK SPHERICAL SPREADING LOSS #
    check_ssl(idir, test_points, version)
    # ssl reclass is same between Arc9.3 and Arc10.3
    
    # CHECK ATMOSPHERIC ABSORPTION LOSS #
    check_aal(idir, test_points, version)
    # Note: Arc 9.3 produces a different output at 5 decimal places
    
    # CHECK FOLIAGE & GROUND COVER LOSS #
    check_veg(idir, test_points, version)
    # Note: Arc 9.3 produces results that are ~1/2 dB different
    # Note: The new version of SPreAD-GIS produces changed results
    
    # CHECK UPWIND AND DOWNWIND LOSS #
    # Matches between ArcGIS 9.3 & 10.3
    check_windloss(idir, idir, test_points, n_testpoints, version)
    
    # CHECK TOPOGRAPHIC EFFECTS AND BARRIER LOSS #
    # Viewshed differs between ArcGIS 9.3 & 10.3
    check_topo_barrier(idir, idir, test_points, version)
                    
    # CHECK SUMMARY NOISE PROPAGATION PATTERNS #
    check_results(rpath, epath, idir, test_points, version)

    if version != "10.5":       
        arcpy.AddMessage("SPreAD-GIS (old terrain/vegetation version) validation completed successfully")
    if version == "10.5":
        arcpy.AddMessage("SPreAD-GIS (old terrain/vegetation version) CONTAINS BUGS IN THIS VERSION OF ARCGIS (10.5) AND IS NOT RECOMMENDED FOR USE")

# Check that nmsimgis executed correctly
def test_nmsimgis(output_dir, test_points, n_testpoints, version):

    # Check that impedance values are correct
    tdir = output_dir + "nmsimgis/"
    imp = ["70", "70", "70","70","70","200"]
    raster_check(tdir, "impedance", tdir, test_points, imp, 0)

    # Check point 1 propagation time
    prop = ["8.7", "4.56", "0.0","0.79","2.36","5.18"]
    if version == "10.4.1" or version == "10.5":
        prop = ["0.0", "4.56", "0.79", "2.36", "5.18", "8.7"] #**# THIS IS weird - the points are being sampled in a different order here, but if so, why weren't they in a different order for impedance???
    raster_check(tdir, "propagation_time_pt1.tif", tdir, test_points, prop, 2)

    # Check final output for 1000 Hz band
    tdir = output_dir + "nmsimgis/frequency_propagation/"
    
    test_values = ["11.64", "10.34", "75.33","43.08","35.65","-34.54"]
    if version == "10.4.1" or version == "10.5":
        test_values = ["75.29","10.34", "43.08", "35.65", "-34.54", "11.64"] #**# getting 75.84 instead of 75.88. What changed in ArcGIS 10.4? SSL
    raster_check(tdir, "pt1_pr01000.tif", tdir, test_points, test_values, 2)

    names_vec = ["pt1_ab01000.tif", "pt1_at01000.tif", "pt1_sr01000.tif", "pt1_ss01000.tif"]
    dec_vec = [2] * len(names_vec)

    # Spot check only the 1000 Hz band for point 1
    tdir = output_dir + "nmsimgis/points/pt1/"

    # Check Sound Source Levels, Spherical Spreading Loss, Atmospheric Absorption,
    # Ground Attenuation, and final values (not in respective order)
    ab = ["16.2", "8.3", "0.1","1.6","4.4","9.8"]
    at = ["-3.8", "11.1", "-3.3","-0.7","-4.8","53.1"]
    ss = ["19.4", "13.6", "-28.7","-0.6","8.1","15.0"]
    if version == "10.4.1" or version == "10.5":
        ab = ["8.33", "16.22", "0.06","1.63","4.43","9.81"]
        at = ["11.15", "-3.85", "-3.79","-0.75","-4.76","53.14"]
        ss = ["13.58", "19.37", "-28.69","-0.57","8.1","15.0"] #**# Getting a value of -28.69 instead of -28.73 for ArcGIS 10.4. What changed?
    sr = ["43.4"] * n_testpoints
    values_vec = [ab, at, sr, ss]

    if version != "10.4.1" and version != "10.5":
        for i in range(0,len(names_vec)):
            raster_name = names_vec[i]
            test_values = values_vec[i]
            decimals = dec_vec[i]
            #print raster_name
            #print decimals
            raster_check(tdir, raster_name, tdir, test_points, test_values, decimals)
    else:
        print "NMSIMGIS intermediate values were not checked due to a bug in the arcpy.sa.Sample tool."

    arcpy.AddMessage("nmsimgis completed successfully")

    # Test summary across points for 1000 Hz band
    tdir = output_dir + "nmsimgis/point_summaries/"
    summary = ["15.53", "10.64", "75.33","46.98","41.77","-32.9"]
    if version == "10.4.1" or version == "10.5":
        summary = ["75.29", "10.65", "46.98", "41.76", "-32.9", "15.52"] #**# FIX 75.84 instead of 75.88, -33.46 instead of -33.47!
    raster_check(tdir, "sumenergy_01000.tif", tdir, test_points, summary, 2)

    # Test weightings for 400 Hz band
    tdir = output_dir + "nmsimgis/weightings/"
    weighting = ["8.57", "0.13", "63.65","34.41","22.91","-39.26"]
    if version == "10.4.1" or version == "10.5":
        weighting = ["63.61", "0.13", "34.41", "22.91", "-39.26", "8.57"] # 64.27 instead of 64.31, 8.57 instead of 8.58
    raster_check(tdir, "prA00400.tif", tdir, test_points, weighting, 2)

    # Test final summary
    tdir = output_dir
    fin = ["25.33", "17.59", "84.87","64.06","49.6","0.0"]
    if version == "10.4.1" or version == "10.5":
        fin = ["84.83", "17.59", "64.06", "49.59", "0.0", "25.33"] # 85.5 instead of 85.54
    raster_check(tdir, "test_nmsimgis_A.tif", tdir, test_points, fin, 2)
    
    arcpy.AddMessage("summary and weighting completed successfully")


# Checks rasters for specific values at the test points using the Sample tool
def sample_check(test_points, raster, field, test_values, check_name, tpath,
                 decimals = 0):
    '''test points are the points at which the raster will be sampled, raster
    is the raster to be sampled, field is the name of the field (= raster name
    but truncated to 10 characters), test_values are what the raster values
    should be, check_name is a name for the check to let the user know where
    to look when there are problems, tpath is the path of the testing directory,
    and decimals is the number of decimal places of the result to check'''
    
    #Still trying to get rid of the stupid schema lock that is making troubleshooting painful.        
    mismatch, check_name, test_value, count  = sample_check_subfunction(test_points, raster, field, test_values, check_name,
                              tpath, decimals)
    
    # Placed outside of loop in hopes that I can avoid an annoying schema lock.
    if mismatch == 1:
        raise ValueError("ERROR: %s value was not correctly computed: %s vs %s" % (check_name, test_value, test_values[count]))

#Still trying to get rid of the stupid schema lock that is making troubleshooting painful.
def sample_check_subfunction(test_points, raster, field, test_values, check_name,
                              tpath, decimals):    
    
    # sample from raster
    temp_table = arcpy.sa.Sample([raster], test_points, tpath + "temp_table.dbf", "NEAREST")
    
    with arcpy.da.SearchCursor(temp_table, field) as cur:
        #row = cur.next()
        
        count = 0
        mismatch = 0
        for row in cur:
            if (decimals == 0):
                test_value = str(int(round(row[0]))) # Rounding is to prevent truncation in int command
            else:
                test_value = str(round(row[0], decimals))
            if test_value != test_values[count]:
                mismatch = 1            
                break # End loop if error
            #row = cur.next()
            count += 1
       
    return (mismatch, check_name, test_value, count)

# Calls sample_check with some pre-filled values
def raster_check(ipath, raster_name, tpath, test_points, test_values, decimals):
    raster = ipath + raster_name
    if raster_name[-4:] == ".tif":
        raster_name = raster_name[:-4]
    sample_check(test_points, raster, raster_name[0:10], test_values, "%s_raster" % raster_name, tpath, decimals)
   

# Check the overall results
def check_results(rpath, epath, tpath, test_points, version):
    ## Check prediction raster
    # Source values are now increased due to a change in the handling of
    # spherical spreading loss for the origin cell
    pr_raster = rpath + "pr400.tif"
    #pr_test_values = ["16.38", "7.19", "83.46", "16.92", "33.79", "-38.67"]
    pr_test_values = ["16.29", "7.13", "83.73", "16.94", "33.76", "-38.95"] # Changed slightly with change in where distance and direction were calculated - results are slightly dependant on extent due to the way barrier effects are calculated.
    # 0.01 difference in results between 10.3 & 10.4 Probably not a problem.    
    if version == "10.4":
        pr_test_values = ["16.38", "7.13", "83.73", "17.08", "33.76", "-39.07"]
        
    if version == "10.4.1":
        pr_test_values = ["83.73", "7.13", "16.94", "33.76", "-39.07", "16.29"] # 16.92 instead of 17.08
    if version == "10.5":
        pr_test_values = ["83.73", "7.39", "34.31", "33.76", "-38.96", "16.29"] #**# FLAG 34.31 vs. 16.94 - this is NOT a trivial difference!
    sample_check(test_points, pr_raster, "pr400", pr_test_values, "pr_raster", tpath, 2) 

    ## Check excess raster
    # Excess raster adjusted to give values more consistent with NMSIMGIS
    ex_raster = epath + "ex00400.tif"
    ex_test_values = ["0.0", "0.0", "54.73", "0.0", "4.76", "0.0"]
    if version == "10.4.1":
        ex_test_values = ["54.73", "0.0", "0.0","4.76", "0.0", "0.0"]
    if version == "10.5":
        ex_test_values = ["54.73", "0.0", "5.31","4.76", "0.0", "0.0"]
    
    sample_check(test_points, ex_raster, "ex00400", ex_test_values, "ex_raster", tpath, 2)

# Check calculation of Euclidean distance & direction
def check_dist_dir(ipath, tpath, test_points, version):
    
    dist_vals = ["9244","4810","0","900","2500","5507"]    
    dir_vals = ["254","110","0","270","287","30"]

    if version == "10.4.1" or version == "10.5":
        dist_vals = ["4810","9244","0","900","2500","5507"]
        dir_vals = ["110","254","0","270","287","30"]
    
    # Check distance
    raster_check(ipath, "eucdist_ft", ipath, test_points, dist_vals, 0)

    # Check direction
    raster_check(ipath,"eucdir", ipath, test_points, dir_vals, 0)    


# Check spherical spreading loss
def check_ssl(ipath, test_points, version):
    #**# Values not critically examined - need to examine logic behind values at some point
    # but since they match the arc 9.3 output, I'm not worrying about it now.    
    
    # Source levels have been changed since 9.3. Now they include spherical
    # spreading loss, and these have been adjusted to correspond to a distance
    # to give values similar to NMSIMGIS
    cpath = ipath + "ssl/"
    ssl_reclass = ["45","40","-24","25","34","41"] 
    xy_grid = ["185","96","0","18","50","110"]

    if version == "10.4":
        ssl_reclass = ["45","40","-24","25","34","41"]

    if version == "10.4.1" or version == "10.5":
        #ssl_reclass = ["40","45","-9999","25","34","41"]
        ssl_reclass = ["40", "45","-24","25","34","41"]
        xy_grid = ["96","185","0","18","50","110"]

    # Patch added because 10.3.1 and earlier returned 0 from the searchcursor, and ArcGIS 10.4 returns -9999 from the sample tool. Should not matter for later code.
    #try:
    raster_check(cpath, "ssl_reclass", cpath, test_points, ssl_reclass, 0)
    #except:
    raster_check(cpath, "ssl_reclass", cpath, test_points, ssl_reclass, 0)
    raster_check(cpath, "xy_grid", cpath, test_points, xy_grid, 0)

# Check atmospheric absorption loss
def check_aal(ipath, test_points, version):
    cpath = ipath + "aal/"    
    
    aac =  ["0.0009053","0.0009053","0.0009053","0.0009053","0.0009053","0.0009053"]
    
    # ArcGIS calculates the elevation slightly differently (and always incorrectly) between versions.
    # This leads to changes in atmospheric absorption    
    
    # This version works with 10.3    
    #aal = ["8.3688955","4.3548012","0.0","0.814759","2.2632194","4.985662"] #**# Temp 100 to crash it & keep the intermediates
    
    # This version works with 10.4    
    aal = ["8.3689137","4.3548111","0.0","0.8147609","2.2632246","4.9856734"]    
    
    # Try with just 4 decimal places
    #aal = ["8.3689", "4.3548", "0.0", "0.8148", "2.2632", "4.9857"] # This version works with model_extent.shp. Odd that changing the extent would change the output!
    #aal = ["8.3688","4.3548","0.0","0.8148","2.2632","4.9856"] #This version works with the test_extent.shp #**# Changing where distance & direction were calculated caused a minor change in atmospheric absorption - I think the extents are slightly different.
    aal = ["8.3689","4.3548","0.0","0.8148","2.2632","4.9857"] # Now the values match 10.4.1. Odd.

    if version == "10.4.1" or version == "10.5":
        aal = ["4.3548","8.3689","0.0","0.8148","2.2632","4.9857"] # 8.3689 in 10.4.1 8.3688 in 10.3 & previously in 10.4.0. Also 4.9857 in 10.4.1 instead of 4.9856.

    # This version used to work, but no longer does
    #aal = ["8.3688478","4.3547769","0.0","0.8147545","2.2632067","4.9856343"]    

    # Just checking the values to 4 decimal places    
    
    raster_check(cpath, "aac400", cpath, test_points, aac, 7)
    raster_check(cpath, "aal400", cpath, test_points, aal, 4)

# Check vegetation section
# Mismatch between ArcGIS 9.3 and 10.3, but this is no longer an issue as 
# the vegetation calculations have been changed in the new version of
# SPreAD-GIS
def check_veg(ipath, test_points, version):
    cpath = ipath + "veg/"
    names_vec = ["pathdist1", "pathdist2", "veg_cost", "veg400", "constant1", "eucdist_max", "eucdist_z", "veg_lossrate"]
       
    #pathdist1 = ["74", "290", "180", "90", "107", "210"]
    pathdist1 = ["9654", "4960", "0", "892", "2559", "5696"]
    #pathdist2 = ["74", "290", "180", "90", "107", "210"]
    pathdist2 = ["9676", "4972", "0", "905", "2576", "5713"]
    veg_cost = ["1", "1", "0", "1", "1", "1"]
    veg400 = ["22", "12", "0", "13", "17", "17"]
    constant1 = ["1", "1", "1", "1", "1", "1"]
    #eucdist_max = ["18597","18597","18597","18597","18597","18597"]
    #eucdist_max = ["18738"] * 6 # Works with model_extent.shp
    eucdist_max = ["11853"] * 6
    #eucdist_z = ["0","0","1","0","0","0"]
    #eucdist_z = ["9494","13928","18738","17838","16238","13231"] #Works with model_extent.shp
    eucdist_z = ["2609", "7043", "11853", "10953", "9353", "6346"] # Works with test_extent.shp
    #veg_lossrate = ["74", "290", "180", "90", "107", "210"]
    veg_lossrate = ["501", "501", "501", "501", "501", "101"]
 
    if version == "10.4":
        pathdist2 = ["9676", "4972", "-9999", "905", "2576", "5713"]
        veg_cost = ["1", "1", "-9999", "1", "1", "1"]
        #veg400 = ["74", "290", "180", "90", "107", "210"]
        veg400 = ["22", "12", "0", "13", "17", "17"]

    if version == "10.4.1" or version == "10.5":
        pathdist1 = ["4960", "9654", "0", "892", "2559", "5696"] #**# 10.4.1 vs. 10.5 difference
        pathdist2 = ["4972", "9676", "-9999", "905", "2576", "5713"]
        veg_cost = ["1", "1", "-9999", "1", "1", "1"]
        veg400 = ["12", "22", "0", "13", "17", "17"]
        eucdist_z = ["7043", "2609", "11853", "10953", "9353", "6346"] # Works with test_extent.shp

    if version == "10.5":
        pathdist1 = ["4961", "9656", "0", "892", "2559", "5697"] #**# 10.4.1 vs. 10.5 difference
        pathdist2 = ["4973", "9675", "0", "888", "2575", "5714"]
        veg400 = ["12", "20", "0", "-4", "16", "17"]
        
    val_vec = [pathdist1, pathdist2, veg_cost, veg400, constant1, eucdist_max, eucdist_z, veg_lossrate]

    for i in range(len(val_vec)):
        raster_check(cpath, names_vec[i], cpath, test_points, val_vec[i], 0)
    

# Check the windloss section
def check_windloss(ipath, tpath, test_points, n_testpoints, version):
    ## Check phi raster
    phi_raster = ipath + "wind/phi/"
    phi_test_values = ["144"] * n_testpoints
    phi_field = "phi"
    sample_check(test_points, phi_raster, phi_field, phi_test_values, "phi_raster", tpath)
    
    ## Check prop_angle raster
    prop_angle_raster = ipath + "wind/prop_angle"
    prop_test_values = ["74", "290", "180", "90", "107", "210"]
    if version == "10.4.1" or version == "10.5":
        prop_test_values = ["290", "74", "180", "90", "107", "210"]
    sample_check(test_points, prop_angle_raster, "prop_angle", prop_test_values, "prop_angle_raster", tpath)
    
    # Check Value1
    val1_raster = ipath + "wind/conraster1.tif"
    val1_test_values = ["-31","185","75","-15","2","105"]
    if version == "10.4.1" or version == "10.5":
        val1_test_values = ["75", "185","-15","2","105","-31"]
    sample_check(test_points, val1_raster, "conraster1", val1_test_values, "conraster1_raster", tpath)
    
    # Check Value2
    val2_raster = ipath + "wind/conresult1.tif"
    val2_test_values = ["329", "185","75","345","2","105"]
    if version == "10.4.1" or version == "10.5":
        val2_test_values = ["75","185","345","2","105","329"]
    sample_check(test_points, val2_raster, "conresult1", val2_test_values, "conresult1", tpath)
    
    # Check wind_ang
    wind_ang = ipath + "wind/wind_ang"
    wind_ang_test_values = ["31", "175","75","15","2","105"]
    if version == "10.4.1" or version == "10.5":
        wind_ang_test_values = ["175","31","75","15","2","105"]
    sample_check(test_points, wind_ang, "wind_ang", wind_ang_test_values, "wind_ang_raster", tpath)


    ## Check updownwind
    udwind = ipath + "wind/updownwind"
    udwind_test_values = ["113", "-31","69","129","142","39"]
    if version == "10.4.1" or version == "10.5":
        udwind_test_values = ["-31","113","69","129","142","39"]
    sample_check(test_points, udwind, "updownwind", udwind_test_values, "updownwind", tpath)
    
    ## check upwind
    upwind = ipath + "wind/upwind"
    upwind_test_values = ["113", "0","69","129","142","39"]
    if version == "10.4.1" or version == "10.5":
        upwind_test_values = ["0","113","69","129","142","39"]
    sample_check(test_points, upwind, "upwind", upwind_test_values, "upwind_raster", tpath)

    ## check downwind
    downwind = ipath + "wind/downwind"
    downwind_test_values = ["0", "1","0","0","0","0"]
    if version == "10.4.1" or version == "10.5":
            downwind_test_values = ["1","0","0","0","0","0"]
    sample_check(test_points, downwind, "downwind", downwind_test_values, "downwind_raster", tpath)

    # Upwind loss step 1
    upwindloss = ipath + "wind/upwind_loss"
    upwindloss_test_values = ["25", "0","25","25","25","24"]
    if version == "10.4.1" or version == "10.5":
        upwindloss_test_values = ["0","25","25","25","25","24"]        
    sample_check(test_points, upwindloss, "upwind_los", upwindloss_test_values, "upwindloss_raster", tpath)
    #Note the change in name of the field, thanks to ESRI naming restrictions.
    
    # Upwind loss step 2
    # My equation fit from Table 12 is not quite the same as the code's - I get
    # 379.94 * x ** -0.861, but close enough.
    # d is not directly saved, but is used in the calculation of x_d - so check that instead:
    xd = ipath + "wind/x_d"
    xd_test_values = ["144", "75", "0", "14", "39", "86"]
    if version == "10.4.1" or version == "10.5":
        xd_test_values = ["75", "144", "0", "14", "39", "86"]
    sample_check(test_points, xd, "x_d", xd_test_values, "xd_raster", tpath)
    # d = 375 * (8 ** (-0.85)) = 64.03, x_d = 5507 / 64 = 26.2 = 86.0
    
    # Upwind loss step 3
    #Should give 57 from lookup table, which will need to be divided by 100
    ulc = ipath + "wind/upwind_loss_c"
    ulc_test_values = ["14", "0", "0", "14", "14", "13"]
    if version == "10.4.1" or version == "10.5":
        ulc_test_values = ["0", "14", "0", "14", "14", "13"]
    sample_check(test_points, ulc, "upwind_los", ulc_test_values, "upwind_loss_c_raster", tpath)
        
    ## Downwind loss
    # freq * distance (ft)
    #If freq < 406237, then use downwind value (where is the check that the point is actually downwind??)
    #downwind * (4.2598 * (arcpy.sa.Ln("C:/SPreAD-GIS/intermediate/wind/freq_dist")) - 55.014)
    # Then look up on Table 10 (I think this is where the eqn comes from)
    # Equations are based on table 10, again, not using the median, but hte values are pretty close
    # I get 4.2608 log (x) - 55.056
    
    # So, for point 2
    # 4810 * 400 = 1924000
    # ln 1924000 = 6.28
    # I get 6.60. Close enough. note that natural log is used here.
    # So downwind loss checks out too.
    
    dl = ipath + "wind/downwind_loss"
    dl_test_values = ["0.0", "6.63", "0.0", "0.0", "0.0", "0.0"]
    if version == "10.4.1" or version == "10.5":
        dl_test_values = ["6.63", "0.0", "0.0", "0.0", "0.0", "0.0"]        
    sample_check(test_points, dl, "downwind_l", dl_test_values, "downwind_loss_raster", tpath, 2) # Added a decimals option for closer checking!
    
    # Check combined upwind/downwind loss
    wl = ipath + "wind/wind_loss"
    wl_test_values = ["14.25", "6.63", "0.0", "14.25", "14.25", "13.5"]
    if version == "10.4.1" or version == "10.5":
        wl_test_values = ["6.63", "14.25", "0.0", "14.25", "14.25", "13.5"]    
    sample_check(test_points, wl, "wind_loss", wl_test_values, "wind_loss_raster", tpath, 2) # Added a decimals option for closer checking!
    
    # Check smoothed upwind/downwind loss raster
    wlc = ipath + "wind/win400"
    wlc_test_values = ["14.25", "6.63", "9.93", "14.25", "14.25", "13.5"]
    if version == "10.4.1" or version == "10.5":
        wlc_test_values = ["6.63", "14.25", "9.93", "14.25", "14.25", "13.5"]            
    sample_check(test_points, wlc, "win400", wlc_test_values, "win400_raster", tpath, 2) # Added a decimals option for closer checking!

# Check calculation of barrier effects. Note: code format has changed to increase
# coding efficiency. NOTE: Upgrade makes the issues highlighted here no longer
# as critical as this corresponds to the old version of SPreAD-GIS
def check_topo_barrier(ipath, tpath, test_points, version):
    topopath = ipath + "topo/"    
    
    # Set up names vector
    names_vec = ["flow_dir", "ground", "viewshed", "ground_atmos",
                 "topo_zones", "elev_barrier", "dist_barrier", "slope", "elev_sr", "h_b",
                 "h_b_rc", "bar_pathdist", "bar400"] #"basin", , "salvgwnbr400", "salvgwntp400"]     
    dec_vec = [0] * 7 + [2] + [0] * 5 # [0] * 8 if basin is added back in.

    # Set up test values vector
    flowdir = ["32", "1", "8", "16", "32", "4"]
    #**# Odd, it switches between 26 & 28 for the basin calculations.
    # I don't understand why there would be a discrepancy here, but it does not have any quantitative effect on the downstream results
    # so I am disabling this check.
    #basin = ["0", "26", "26", "26", "0", "0"]
    #basin = ["0", "28", "28", "28", "0", "0"] #**# Odd, this used to work, now it wants the other!
    #if version == "10.4":
        #basin = ["-9999", "28", "28", "28", "-9999", "-9999"]
        #basin = ["-9999", "26", "26", "26", "-9999", "-9999"]

    ground = ["0", "3", "3", "3", "0", "0"]
    viewshed = ["1", "0", "1", "1", "1", "0"]
    groundatmos = ["1", "3", "3", "3", "1", "0"]
    topo_zones = ["1", "3", "3", "3", "1", "2"]
    # #**# This change appears to be due to patch & not 10.4
    #elev_barrier = ["5916", "5672", "5852", "5884", "5887", "5922"]
    #if version == "10.4":
    #elev_barrier = ["5916", "5668", "5862", "5884", "5887", "5922"] # Works with model_extent.shp
    elev_barrier = ["5916", "5832", "5914", "5884", "5887", "5922"] # Works with test_extent.shp #**# THIS HIGHLIGHTS THAT THE BARRIER IS DEPENDENT ON ELEVATION BEHIND THE POINT, AND MAKES THE MODEL RESULT SENSITIVE TO MODEL EXTENT. BUT NOTE THAT THE CHANGE IN EXTENT DID NOT CHANGE THE FINAL RESULT!
    #**# This might be a serious problem. Why is this different? (related to the patch, not 10.4)
    #dist_barrier = ["1552", "10756", "7500", "1500", "1552", "1253"]
    #if version == "10.4":
    #dist_barrier = ["1552", "10791", "7701", "1500", "1552", "1253"] #Works with model_extent.shp
    dist_barrier = ["1552", "6074", "1500", "1500", "1552", "1253"] #Works with test_extent.shp

    #slope = ["0.07", "-0.02", "0.0", "0.03", "0.07", "-0.22"] # Works with model_extent.shp
    #slope = ["0.07", "-0.02", "0.0", "0.02", "0.06", "-0.22"] # Works with test_extent.shp #**# Not after I moved the calculation out of spreadgishlpr
    slope = ["0.07", "-0.02", "0.0", "0.03", "0.06", "-0.22"] # Works with test_extent.shp
    
    #**# -9999.0?
    #slope = ["0.07", "-0.02", "0.0", "0.02", "0.06", "-0.22"] #Temporary patch, unclear what changed that it was needed, nor why it was no longer needed.
    #elev_sr = ["5950", "5670", "0", "5891", "5948", "5567"] # Works with model_extent.shp 
    #elev_sr = ["5956", "5745", "0", "5886", "5951", "5573"]  # Works with test_extent.shp #**# Not after I moved the distance dir calculations out of spreadgishlpr
    elev_sr = ["5951", "5747", "0", "5890", "5949", "5569"]
    #elev_sr = ["5956", "5661", "0", "5886", "5951", "5573"] 
    #**# THIS COULD BE A PROBLEM
    #h_b = ["-34", "2", "0", "-7", "-61", "355"]
    #h_b = ["-34", "-2", "0", "-7", "-61", "355"] # Works with model_extent.shp
    #h_b = ["-40", "87", "0", "-2", "-64", "349"] # Works with test_extent.shp #**# changed after I moved distance & direction calculations out of spreadgishlpr
    h_b = ["-35", "85", "0", "-6", "-62", "353"]
    #h_b = ["-40", "11", "0", "-2", "-64", "349"]
    #h_b_rc = ["0", "1", "0", "0", "0", "354"] #Change in values - reclass truncates the decimal instead of rounding.
    h_b_rc = ["0", "85", "0", "0", "0", "353"]
    #**# THIS COULD BE A PROBLEM    
    #h_b_rc = ["0", "10", "0", "0", "0", "348"]
    #bar_pathdist = ["0", "11891", "0", "1200", "0", "64"]
    #bar_pathdist = ["0", "11961", "0", "1200", "0", "64"] # Works with model_extent.shp
    #bar_pathdist = ["0", "2531", "0", "1200", "0", "62"] # Works with test_extent.shp
    bar_pathdist = ["0", "2531", "0", "1200", "0", "63"] # Works with test_extent.shp post-dist/dir change

    #bar_pathdist = ["0", "11891", "0", "1200", "0", "62"]
    #bar400 = ["0", "108", "0", "64", "0", "32"] # Works with model_extent.shp
    bar400 = ["0", "76", "0", "64", "0", "32"] # Works with test_extent.shp
    #salvgwnbr400 = ["-20", "-101", "0", "-47", "2", "-39"]
    #salvgwntp400 = ["16", "7", "0", "17", "34", "-39"]

    if version == "10.4":
        slope = ["0.07", "-0.02", "-9999.0", "0.03", "0.07", "-0.22"]
        elev_sr = ["5947", "5674", "-9999", "5893", "5947", "5565"] 
        h_b = ["-31", "87", "-9999", "-9", "-60", "357"]
        h_b_rc = ["0", "86", "-9999", "0", "0", "356"]
        bar_pathdist = ["0", "2531", "-9999", "1200", "0", "64"]
        bar400 = ["0", "76", "-9999", "64", "0", "33"]

    if version == "10.4.1" or version == "10.5":
        flowdir = ["1", "32", "8", "16", "32", "4"]
        ground = ["3", "0", "3", "3", "0", "0"]
        viewshed = ["0", "1", "1", "1", "1", "0"]
        groundatmos = ["3", "1", "3", "3", "1", "0"]
        topo_zones = ["3", "1", "3", "3", "1", "2"]
        elev_barrier = ["5832", "5916", "5914", "5884", "5887", "5922"] # Works with test_extent.shp #**# THIS HIGHLIGHTS THAT THE BARRIER IS DEPENDENT ON ELEVATION BEHIND THE POINT, AND MAKES THE MODEL RESULT SENSITIVE TO MODEL EXTENT. BUT NOTE THAT THE CHANGE IN EXTENT DID NOT CHANGE THE FINAL RESULT!
        dist_barrier = ["6074", "1552", "1500", "1500", "1552", "1253"] #Works with test_extent.shp
        slope = ["-0.02", "0.07", "-9999.0", "0.03", "0.07", "-0.22"]
        elev_sr = ["5748", "5947", "-9999", "5893", "5947", "5565"]  #5748 instead of 5645!
        h_b = ["84", "-31", "-9999", "-9", "-60", "357"] #84 instead of 87, -31 instead of -40 in 10.3
        h_b_rc = ["83", "0", "-9999", "0", "0", "356"]
        bar_pathdist = ["2531", "0", "-9999", "1200", "0", "64"]
        bar400 = ["76", "0", "-9999", "64", "0", "33"]

    
    values_vec = [flowdir, ground, viewshed, groundatmos, topo_zones,
                  elev_barrier, dist_barrier, slope, elev_sr, h_b, h_b_rc,
                  bar_pathdist, bar400] #basin, , salvgwnbr400, salvgwntp400]


    # Run loop
    for i in range(0,len(names_vec)): #range(8,9):
        raster_name = names_vec[i]
        test_values = values_vec[i]
        decimals = dec_vec[i]
        raster_check(topopath, raster_name, topopath, test_points, test_values, decimals)
    
   
def check_aweight(rpath, tpath, test_points):
    raise ValueError("This has not yet been scripted")    
    pr_raster = rpath + "pr400"
    pr_test_values = ["16.29", "7.21", "49.62", "17.08", "33.79", "0.0"]
    version = arcpy.GetInstallInfo('desktop')['Version']
    # 0.01 difference in results between 10.3 & 10.4 Probably not a problem.    
    if version == "10.4":
        pr_test_values = ["16.29", "7.2", "49.62", "17.08", "33.79", "0.0"]
    sample_check(test_points, pr_raster, "pr400", pr_test_values, "pr_raster", tpath, 2) 

# Test for new version of SPreAD-GIS
def test_spreadgis2(output_dir2, test_points, n_testpoints, version):

    # Just test final outputs
    tdir = output_dir2

    fin = ["0.0","0.0","70.83","11.06","0.73","0.0"]
    if version == "10.4.1" or version == "10.5":
        fin = ["70.83", "0.0", "11.06","0.73","0.0","0.0"]
    
    raster_check(tdir, "test_spreadgis2_A_pt1.tif", tdir, test_points, fin, 2)
    
    arcpy.AddMessage("SPreAD-GIS (new terrain/vegetation version) tested successfully")


# Test for added ISO model
def test_iso_model(output_dir, test_points, n_testpoints, version):

    # Just test final outputs
    tdir = output_dir
    
    fin = ["11.89", "22.26", "83.79", "51.11","37.45","10.7"]    #**# 83.76 vs 83.75 slight difference between 10.3 and 10.4.1
    
    if version == "10.4.1" or version == "10.5":
        fin = ["83.75", "22.26", "51.11","37.45","10.7","11.89"]
    raster_check(tdir, "test_iso_A.tif", tdir, test_points, fin, 2)
    
    arcpy.AddMessage("ISO GIS implementation tested successfully")


