# -*- coding: utf-8 -*-
"""
Description: To Calculate The Auditory Dectectability Index d' as a Function of
1/3 Octave Band Target and Background Spectra. This has been adapted to
calculate dPrime separately for 1/3 octave bands, so that raster tools can be
applied to the calculations

Authors: A.C. Keyel derived this script from previous work done by
R. Horonjeff and KJP (assumed to be Keith Plotkin).

Authorship details:
This script was adapted from dprime.py, which was translated into Python
from Dprime.for from 09/2016 - 02/2017 by A.C. Keyel. The results were tested
for equivalence using dprime_test.py and a .dll of the dPrime algorithm created
by Bruce Ikelheimer. The two versions matched to at least 3 decimal places.

Dprime.for was provided to A.C. Keyel under the GPL 2.0 license by
Bruce Ikelheimer, consequently this version is also licensed under GPL 2.0.

The author of Dprime.for is R. Horonjeff, 29 September 2000, with revisions on
4 Dec 2004 and 17 Dec 2004 by KJP (assumed to be Keith Plotkin).

Dprime.for was derived from 1983 Acoustic Range Prediction Program (ARPP)
software developed under USAF contract F33615-83-C-3216 by Bolt Beranek and
Newman, Inc (BBN) for:                        *

FLIGHT DYNAMICS LABORATORY
AIR FORCE WRIGHT AERONAUTICAL LABORATORIES
AIR FORCE SYSTEMS COMMAND
WRIGHT-PATTERSON AFB,  OHIO  45433

Documentation for the original may by found in USAF Technical Report
No. AFWAL-TR-83-3115 which contains the Program Description and Users' Guide.

KJP revised the EASN method, the hearing threshold (EARSPC from 2002 ISO
standard). EASN is now based on d' = 1.75, instead of 1.5.

METHOD:   Uses basic auditory signal detection equation:
              D' = ETA * SQRT(BW) * S/N                       
                 where: ETA = observer efficiency             
                        BW  = masking bandwidth (usually      
                              1/3 octave bandwidth)           
                        S   = target signal intensity,        
                              exp10(SPL/10)                   
                        N   = background intensity,           
                              exp10(SPL/10)                   
                                                             
CALLING ARGUMENTS:
point_source_file: A shape file containing sound source points.

source_id_field: A unique ID field for the point sources

n_points: The number of points for which to run the model. "all" will run the
        model for all points

point_fill: The number of leading zeros before the point in the propagation
        results file names

bands: Audibility will be calculated based on the 1/3 octave bands specified
        in this field. NMSIMGIS will run for the 1/3 octave bands from
        50 - 10,000 Hz, SPreAD-GIS will run for 125 - 2000 Hz, and HARRISON
        will run for 400 - 2000 Hz. Frequency bands < 50 or >10,000 Hz are not\
        supported

band_fill: The number of leading zeros before the point in the propagation
        results file names

onethird_dir: The directory containing the sound propagation results
        corresponding to the point specified above. Typically, this will be
        in the frequency_propagation subfolder in the model results folder.
        (e.g., C:/smt/spreadgis/frequency_propagation)

ambient_dir: The directory containing ambient values. Ambient values must exist
        for each band in the bands parameter.

dprime_path: The output directory.

RETURNED VARIABLES:                                         
Spatially-explicit rasters depicting the cumulative
         dprime value across all specified octave
         bands are produced in the dprime_path
         directory.                                


"""


import shutil
import arcpy, arcpy.sa
import soundprophlpr as sp
arcpy.CheckOutExtension("Spatial")

# Calculate audibility (dprime) for one or more sound sources
def calculate_dprime(point_source_file, source_id_field, n_points, point_fill, bands, band_fill, onethird_dir, ambient_dir, dprime_path):
    '''
    point_source_file:  A list of points at which to evaluate audibility
    source_id_field:    A field containing a unique ID for each source
    n_points:           The number of points for which to run the model. "all"
                        will run for all points in the shapefile.
    point_fill:         The number of leading zeros before the point label.
    bands:              The frequency bands to use in the dPrime analysis.
                        Frequency bands < 50 or >10,000 Hz are not supported
                        In soundprophlpr.py these are called frequencies.
    band_fill:          The number of leading zeros before the frequency band
    onethird_dir:       The directory with the 1/3 octave band sound levels
    ambient_dir:        The directory with the 1/3 octave band ambient levels
                        Currently, the ambient files MUST be in ESRI GRD format
    dprime_path:        The directory to contain the audibility results
    '''
        
    # Create a temporary path as the sum is computed
    # Needed because over-writing rasters here has led to crashes
    temp_path = dprime_path + "temp/"
    sp.make(temp_path)
    
    # Create a path for band-specific audibilities
    dp_path = dprime_path + "dp/"
    sp.make(dp_path)    
    
    # Extract point list from point source file
    point_lst = sp.get_points(point_source_file, source_id_field, n_points)
    
    # Calculate dPrime for each point    
    for point in point_lst:
    
        # Get correct number of leading zeros for point
        point_label = str(point).zfill(point_fill)
    
        # Calculate dPrime for each 1/3 octave band
        for bb in xrange(len(bands)):
            band = bands[bb]
            band_label = str(band).zfill(band_fill)
            propagation_values = onethird_dir + "pt%s_pr%s.tif" % (point_label, band_label)
            ambient_values = ambient_dir + "ambient%s" % band_label #**# May want to allow formats other than grd
            
            out_raster = calculate_band_audibility(band, propagation_values, ambient_values, dp_path)

            #**# Not sure why it is squared, but that is what the original code does
            power_raster = arcpy.sa.Power(out_raster, 2) 
            power_raster.save(temp_path + "pow_%s_%s.tif" % (point, band))

            new_sum_raster = temp_path + "temp_pt%s_band%s.tif" % (point, band_label)
            if bb == 0:
                sum_raster = power_raster                
                sum_raster.save(new_sum_raster)
            else:
                sum_raster = sum_raster + power_raster
                sum_raster.save(new_sum_raster)
            
        # Save the final sum with a more intuitive name
        dpss = dprime_path + "DPSS_pt%s.tif" % point_label
        dpss_raster = arcpy.sa.Power(sum_raster, 0.5) # 0.5 is same as square root
        dpss_raster.save(dpss)
        
        # Fidell, Sanford, K. Pearsons, and M. Sneddon. Evaluation of the
        # effectiveness of SFAR 50-2 in restoring natural quiet to Grand Canyon
        # National Park. Report NPOA-93-1, BBN-7197 (NTIS Number: PB95-195202),
        # 1994.
        # Fidell et al. propose a threshold of 10 log10(d') of 7 (See section
        # 4.8 on p.55 of the report (p. 67 of the pdf)
        # NOTE: They round 7.3 down to 7
        # Here we return 10 log10(d') as well, as this will be the more useful
        # quantity.
        dpss_10log10_file = dprime_path + "DPSS_10LOG10_pt%s.tif" % point_label
        dpss_10log10 = 10 * arcpy.sa.Log10(dpss_raster)
        dpss_10log10.save(dpss_10log10_file)

# Calculate audibility (dprime) for a summary file of multiple sound sources
def calculate_dprime_summary(bands, band_fill, onethird_dir, ambient_dir, dprime_path):
    '''
    bands:        The frequency bands to use in the dPrime analysis.
                  Frequency bands < 50 or >10,000 are not supported
                  In soundprophlpr.py these are called frequencies.
    band_fill:    The number of leading zeros before the frequency band
    onethird_dir: The directory with the 1/3 octave band sound levels
    ambient_dir:  The directory with the 1/3 octave band ambient levels
    dprime_path:  The directory to contain the audibility results
    '''
        
    # Create a temporary path as the sum is computed (over-writing rasters here
    # has led to crashes)    
    temp_path = dprime_path + "temp/"
    sp.make(temp_path)
    
    # Create a bath for band-specific audibilities
    dp_path = dprime_path + "dp/"
    sp.make(dp_path)    
        
    # Calculate dPrime for each 1/3 octave band
    for bb in xrange(len(bands)):
        band = bands[bb]
        band_label = str(band).zfill(band_fill)
        propagation_values = onethird_dir + "sumenergy_%s.tif" % (band_label)
        ambient_values = ambient_dir + "ambient%s" % band_label #**# May want to allow formats other than grd
        
        out_raster = calculate_band_audibility(band, propagation_values, ambient_values, dp_path)

        #**# Not sure why it is squared, but that is what the original code does
        power_raster = arcpy.sa.Power(out_raster, 2) 
        power_raster.save(temp_path + "pow_%s.tif" % (band))

        new_sum_raster = temp_path + "temp_band%s.tif" % (band_label)
        if bb == 0:
            sum_raster = power_raster                
            sum_raster.save(new_sum_raster)
        else:
            sum_raster = sum_raster + power_raster
            sum_raster.save(new_sum_raster)
        
    # Save the final sum with a more intuitive name
    dpss = dprime_path + "DPSS_summary.tif"
    dpss_raster = arcpy.sa.Power(sum_raster, 0.5) # 0.5 is same as square root
    dpss_raster.save(dpss)

    # Fidell, Sanford, K. Pearsons, and M. Sneddon. Evaluation of the
    # effectiveness of SFAR 50-2 in restoring natural quiet to Grand Canyon
    # National Park. Report NPOA-93-1, BBN-7197 (NTIS Number: PB95-195202),
    # 1994.
    # Fidell et al. propose a threshold of 10 log10(d') of 7 (See section
    # 4.8 on p.55 of the report (p. 67 of the pdf)
    # NOTE: They round 7.3 down to 7
    # Here we return 10 log10(d') as well, as this will be the more useful
    # quantity.
    dpss_10log10_file = dprime_path + "DPSS_10LOG10_summary.tif"
    dpss_10log10 = 10 * arcpy.sa.Log10(dpss_raster)
    dpss_10log10.save(dpss_10log10_file)


    arcpy.AddMessage("Dprime summary raster created successfully")

# Calculate dprime for each band
def calculate_band_audibility(band, propagation_values, ambient_values, dp_path):
    '''
    band: target 1/3 octave band
    propagation_values: Raster with sound levels
    ambient_values: Raster with background levels
    dp_path: path for band-specific audibility
    '''

    # Create a temp path to hold intermediates
    temp_path = dp_path + "temp/"
    sp.make(temp_path)

    # Check that this band is supported, and if so, get an index value for it

    ab_part1 = [50, 63, 80, 100, 125, 160, 200, 250]
    ab_part2 = [315, 400, 500, 630, 800, 1000, 1250, 1600]
    ab_part3 = [2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000]
    allowed_bands = ab_part1 + ab_part2 + ab_part3
    
    # Watch to ensure that both band and a_band are in the correct format!
    index = "NA"
    for i in xrange(len(allowed_bands)):
        a_band= allowed_bands[i]
        if int(a_band) == int(band):
            index = i
            break # Stops loop, such that i corresponds to index

    # check if index was assigned, otherwise give a value error
    if index == "NA":
        m1 = "Only 1/3 octave bands between 50 Hz and 10 kHz are supported. "
        m2 = "Your input was %s. If this corresponds to a supported 1/3" % band
        m3 = " octave band, please check your input format."
        raise ValueError(m1 + m2 + m3)

    # NOTE: The use of lists below is an artifact of when the 1/3 octave values
    # were computed within a loop. It was faster to keep the structure, and just
    # extract the results corresponding to the appropriate index, i

    # initialize variables
    easn = ["NA"]*24

    # Threshold of hearing array from ISO 389-7
    part1 = [43.5, 37.5, 31.5, 26.5, 22.0, 18.0, 14.5, 11.0, 8.5, 5.5, 3.5, 1.5]
    part2 = [1.0, 0.5, 0.0, -1.0, -1.5, -3.0, -4.5, -5.0, -3.5, 0.5, 5.3, 9.5]
    earspc = part1 + part2
        
    # Compute EASN, the effective noise associated with EARSPC
    dpnoise = 1.75

    # use the earspc function to return a sound pressure level, and convert it to energy
    s = 10.**(earspc[i] * 0.1)

    # Assign the value to the corresponding effective noise for the frequency band
    I = i + 1 # Convert to a 1 - 24 range from 0 - 23
    easn[i] = etabw(I,0.0) * (s / dpnoise)
   
    # Calculate value of d' for each cell in the raster
    asl_raster = arcpy.Raster(propagation_values)
    bkgnd_raster = arcpy.Raster(ambient_values)

    # New method - noise is combination of background noise and ear noise

    #sig = 10.**(0.1 * asl)
    asl_tenth = asl_raster * 0.1
    asl_tenth.save(temp_path + "asl_tenth.tif")    
    sig_raster = arcpy.sa.Power(10, asl_tenth)
    sig_raster.save(temp_path + "sig.tif")

    #bnoise = 10**(0.1 * bkgnd)    
    bkgnd_tenth = bkgnd_raster * 0.1
    bkgnd_tenth.save(temp_path + "bkgnd_tenth.tif")    
    bnoise_raster = arcpy.sa.Power(10, bkgnd_tenth)
    bnoise_raster.save(temp_path + "bnoise.tif")
    
    #dp = etabw(I, bkgnd) * (sig / (bnoise + easn[i]))
    term1 = etabw_raster(I, bkgnd_raster, temp_path)
    term2 = sig_raster / (bnoise_raster + easn[i])
    dp = term1 * term2
    dp.save(dp_path + "dp_band%s.tif" % band)

    # Remove temporary intermediates
    #shutil.rmtree(temp_path)

    return (dp)



# Product of Detector Efficiency function
# This is used for calculating easn variable.
# A derivative is used in calculation of dp
'''
C      ****************************************************************
C      *                                                              *
C      *  PURPOSE:  To Generate the Product of Detector Efficiency    *
C      *            (ETA) and the Square Root of the Bandwidth As a   *
C      *            Function of Band Number and Background Sound      *
C      *            Level                                             *
C      *                                                              *
C      *  METHOD:   The Auditory Masking Bandwidth Can Be Greater     *
C      *            Than a 1/3 Octave Band, Especially at Low         *
C      *            Frequencies and High Background Sound Levels      *
C      *            Most of the Code Has To Do With Determining When  *
C      *            This Occurs and Making Appropriate Adjustments    *
C      *                                                              *
C      *  CALL SEQUENCE:                                              *
C      *            VAR = ETABW (I, BKGLVL)                           *
C      *                                                              *
C      *  CALLING ARGUMENTS:                                          *
C      *            I      = Band Number (1 = 50 Hz, 24 = 10 KHz)     *
C      *            BKGLVL = 1/3 Octave Band Background SPL           *
C      *                                                              *
C      ****************************************************************
'''
def etabw(I, BKGLVL):
    
    # Fortran requires specification of the dimensions before the variable
    # is given values. Python does this on the fly. This part of the code was
    # not translated
    # DIMENSION   NP(24), TLETA(24), TLBW(2,6,24)
    # DIMENSION   TLBW1(12), TLBW2(12), TLBW3(12), TLBW4(12)
    # DIMENSION   TLBW5(12), TLBW6(12), TLBW7(12), TLBW8(12)
    # DIMENSION   TLBW9(12), TLBW10(12), TLBW11(12), TLBW12(12)
    # DIMENSION   TLBW13(12), TLBW14(12), TLBW15(12), TLBW16(12)
    # DIMENSION   TLBW17(12), TLBW18(12), TLBW19(12), TLBW20(12)
    # DIMENSION   TLBW21(12), TLBW22(12), TLBW23(12), TLBW24(12)

    # EQUIVALENCE (TLBW(1,1,1),TLBW1(1))
    # EQUIVALENCE (TLBW(1,1,2),TLBW2(1))
    # EQUIVALENCE (TLBW(1,1,3),TLBW3(1))
    # EQUIVALENCE (TLBW(1,1,4),TLBW4(1))
    # EQUIVALENCE (TLBW(1,1,5),TLBW5(1))
    # EQUIVALENCE (TLBW(1,1,6),TLBW6(1))
    # EQUIVALENCE (TLBW(1,1,7),TLBW7(1))
    # EQUIVALENCE (TLBW(1,1,8),TLBW8(1))
    # EQUIVALENCE (TLBW(1,1,9),TLBW9(1))
    # EQUIVALENCE (TLBW(1,1,10),TLBW10(1))
    # EQUIVALENCE (TLBW(1,1,11),TLBW11(1))
    # EQUIVALENCE (TLBW(1,1,12),TLBW12(1))
    # EQUIVALENCE (TLBW(1,1,13),TLBW13(1))
    # EQUIVALENCE (TLBW(1,1,14),TLBW14(1))
    # EQUIVALENCE (TLBW(1,1,15),TLBW15(1))
    # EQUIVALENCE (TLBW(1,1,16),TLBW16(1))
    # EQUIVALENCE (TLBW(1,1,17),TLBW17(1))
    # EQUIVALENCE (TLBW(1,1,18),TLBW18(1))
    # EQUIVALENCE (TLBW(1,1,19),TLBW19(1))
    # EQUIVALENCE (TLBW(1,1,20),TLBW20(1))
    # EQUIVALENCE (TLBW(1,1,21),TLBW21(1))
    # EQUIVALENCE (TLBW(1,1,22),TLBW22(1))
    # EQUIVALENCE (TLBW(1,1,23),TLBW23(1))
    # EQUIVALENCE (TLBW(1,1,24),TLBW24(1))
    
    # 1/3 octave band values of 10*log10(eta)
    part1 = [-6.96, -6.26, -5.56, -5.06, -4.66, -4.36, -4.16, -3.96, -3.76, -3.56, -3.56, -3.56]
    part2 = [-3.56, -3.56, -3.76, -3.96, -4.16, -4.36, -4.56, -4.96, -5.36, -5.76, -6.26, -6.86]
    TLETA = part1 + part2
    
    # Masking bandwidths as a function of frequency band and background
    # spectrum level (dB/Hz)
    NP = [2] * 3 + [4] * 9 + [5] * 3 + [3] * 2 + [4] * 5 + [3] * 2
    
    # Here I modify the FORTRAN to be easier to follow by adding a function with
    # a lot of cumbersome if statements
    TLBW_x, TLBW_y = TLBW_LOOKUP(I)
    
    # Calculate 10*LOG(BANDWIDTH) for 1/3 and 1/1 Octave Bands
    # Center Frequency Multipliers are:
    # 1/3OB = 10^(+0.5/10) - 10^(-0.5/10) = 0.230768
    # 1/1OB = 10^(+1.5/10) - 10^(-1.5/10) = 0.704592
    
    TL13BW = float(I + 16) - 6.368253 # 
    TL11BW = float(I + 16) - 1.520624
    
    #
    # Calculate Average Spectrum Level in Band By Subtracting
    # 10*LOG(BANDWIDTH) From 1/3 Octave SPL
    #
    
    BN0 = BKGLVL - TL13BW

    #
    # Look Up 10*LOG(MASKING BANDWIDTH) As A Function of
    # Spectrum Level (N0)
    #

    # Curve function is defined below
    TLMBW = CURVE_function(BN0, TLBW_x, TLBW_y, NP[I-1]) # I - 1 to get back to i, for use in indexing
    
    # Bound the Masking Bandwidth by the 1/3 Octave Bandwidth
    # (lower bound) and the 1/1 Octave Bandwidth (upper bound)
    # to form usable bandwidth.

    # AMIN1 and AMAX1 are just min & max
    TLUBW = min(max(TLMBW, TL13BW), TL11BW)
    
    # Calculate ETA * SQRT(BANDWIDTH)
    ETABW = 10. ** ((TLETA[I - 1] + 0.5 * TLUBW) / 10.) # I - 1 converts to i

    return ETABW

# Raster version of etabw
def etabw_raster(I, BKGLVL_raster, temp_path):
        
    # 1/3 octave band values of 10*log10(eta)
    part1 = [-6.96, -6.26, -5.56, -5.06, -4.66, -4.36, -4.16, -3.96, -3.76, -3.56, -3.56, -3.56]
    part2 = [-3.56, -3.56, -3.76, -3.96, -4.16, -4.36, -4.56, -4.96, -5.36, -5.76, -6.26, -6.86]
    TLETA = part1 + part2
    
    # Masking bandwidths as a function of frequency band and background
    # spectrum level (dB/Hz)
    NP = [2] * 3 + [4] * 9 + [5] * 3 + [3] * 2 + [4] * 5 + [3] * 2
    
    # Here I modify the FORTRAN to be easier to follow by adding a function with
    # a lot of cumbersome if statements
    TLBW_x, TLBW_y = TLBW_LOOKUP(I)
    
    # Calculate 10*LOG(BANDWIDTH) for 1/3 and 1/1 Octave Bands
    # Center Frequency Multipliers are:
    # 1/3OB = 10^(+0.5/10) - 10^(-0.5/10) = 0.230768
    # 1/1OB = 10^(+1.5/10) - 10^(-1.5/10) = 0.704592
    
    TL13BW = float(I + 16) - 6.368253 # 
    TL11BW = float(I + 16) - 1.520624
    
    #
    # Calculate Average Spectrum Level in Band By Subtracting
    # 10*LOG(BANDWIDTH) From 1/3 Octave SPL
    #
    
    BN0_raster = BKGLVL_raster - TL13BW
    BN0_raster.save(temp_path + "BN0_%s.tif" % I)

    #
    # Look Up 10*LOG(MASKING BANDWIDTH) As A Function of
    # Spectrum Level (N0)
    #

    # Curve function is defined below
    #TLMBW = CURVE_function(BN0, TLBW_x, TLBW_y, NP[I-1]) # I - 1 to get back to i, for use in indexing
    TLMBW_raster_file = CURVE_raster_function(BN0_raster, TLBW_x, TLBW_y, NP[I-1], temp_path) # I - 1 to get back to i, for use in indexing
    
    # Bound the Masking Bandwidth by the 1/3 Octave Bandwidth
    # (lower bound) and the 1/1 Octave Bandwidth (upper bound)
    # to form usable bandwidth.

    #Re-create with a con statement to compare two rasters
    # AMIN1 and AMAX1 are just min & max
    #TLUBW = min(max(TLMBW, TL13BW), TL11BW)

    TLMBW_raster = arcpy.Raster(TLMBW_raster_file)
    TLUBW_term1 = arcpy.sa.Con(TLMBW_raster > TL13BW, TLMBW_raster, TL13BW)
    TLUBW_term1.save(temp_path + "TLUBW_term1.tif")
    TLUBW = arcpy.sa.Con(TLUBW_term1 < TL11BW, TLUBW_term1, TL11BW)
    TLUBW.save(temp_path + "TLUBW.tif")

    # Calculate ETA * SQRT(BANDWIDTH)
    #ETABW = 10. ** ((TLETA[I - 1] + 0.5 * TLUBW) / 10.) # I - 1 converts to i
    ETABW_term1 = ((TLETA[I - 1] + 0.5 * TLUBW) / 10.) # I - 1 converts to i
    ETABW_raster = arcpy.sa.Power(10, ETABW_term1)
    ETABW_raster.save(temp_path + "etabw_raster.tif")

    return ETABW_raster


# Add a function to look up the threshold level for each background level
def TLBW_LOOKUP(I):
    # Check if I is above or below the allowed values
    if I > 24 or I < 1:
        raise ValueError("I = %s. This value is not allowed" % I)

    # NOTE: There were trailing 0's in the Fortran version. These crash the
    # Python implementation, and were removed.

    # each I corresponds to a frequency band    
    if I == 1:    
        TLBW_x = [40, 100] # + [0.0] * 4
        TLBW_y = [12.1, 21.1] # + [0.0] * 4
    if I == 2:
        TLBW_x = [40, 100] # + [0.0] * 4
        TLBW_y = [12.3, 21.3] # + [0.0] * 4
    if I == 3:
        TLBW_x = [40, 100] # + [0.0] * 4
        TLBW_y = [12.6, 21.6] # + [0.0] * 4
    if I == 4:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [13.1, 14.4, 15.8, 24.8] # + [0.0] * 2
    if I == 5:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [13.5, 14.6, 16.1, 25.1] # + [0.0] * 2
    if I == 6:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [14.1, 15.1, 16.4, 25.4] # + [0.0] * 2
    if I == 7:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [14.7, 15.5, 16.8, 25.8] # + [0.0] * 2
    if I == 8:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [15.4, 16.2, 17.2, 26.2] # + [0.0] * 2
    if I == 9:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [16.1, 16.8, 17.8, 26.8] # + [0.0] * 2
    if I == 10:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [16.9, 17.6, 18.4, 27.4] # + [0.0] * 2 
    if I == 11:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [17.6, 18.3, 19.1, 28.1] # + [0.0] * 2
    if I == 12:
        TLBW_x = [40, 50, 60, 120] # + [0.0] * 2
        TLBW_y = [18.4, 19.2, 20.1, 29.1] # + [0.0] * 2
    if I == 13:
        TLBW_x = [30, 40, 50, 60, 120] # + [0.0]
        TLBW_y = [19.1, 19.3, 20.2, 21.2, 31.3] # + [0.0]
    if I == 14:
        TLBW_x = [30, 40, 50, 60, 120] # + [0.0]        
        TLBW_y = [19.8, 20.2, 21.2, 22.3, 31.3] # + [0.0]
    if I == 15:
        TLBW_x = [30, 40, 50, 60, 120] # + [0.0]        
        TLBW_y = [20.5, 21.3, 22.3, 23.6, 32.6] # + [0.0]
    if I == 16:
        TLBW_x = [30, 40, 120] # + [0.0] * 3
        TLBW_y = [21.4, 22.1, 34.1] # + [0.0] * 3
    if I == 17:
        TLBW_x = [30, 40, 120] # + [0.0] * 3
        TLBW_y = [22.1, 23.5, 35.5] # + [0.0] * 3
    if I == 18:
        TLBW_x = [20, 30, 40, 120] # + [0.0] * 2
        TLBW_y = [22.9, 23.4, 24.8, 36.8] # + [0.0] * 2
    if I == 19:
        TLBW_x = [20, 30, 40, 120] # + [0.0] * 2
        TLBW_y = [24.0, 24.7, 26.2, 38.2] # + [0.0] * 2
    if I == 20:
        TLBW_x = [20, 30, 40, 120] # + [0.0] * 2
        TLBW_y = [25.1, 26.2, 27.6, 39.6] # + [0.0] * 2
    if I == 21:
        TLBW_x = [20, 30, 40, 120] # + [0.0] * 2
        TLBW_y = [26.2, 27.5, 29.0, 41.0] # + [0.0] * 2
    if I == 22:
        TLBW_x = [20, 30, 40, 120] # + [0.0] * 2
        TLBW_y = [27.5, 28.9, 30.4, 42.4] # + [0.0] * 2
    if I == 23:
        TLBW_x = [20, 40, 120] # + [0.0] * 3
        TLBW_y = [28.8, 31.8, 43.8] # + [0.0] * 3
    if I == 24:
        TLBW_x = [20, 40, 120] # + [0.0] * 3
        TLBW_y = [30.1, 33.1, 45.1] # + [0.0] * 3

    return (TLBW_x, TLBW_y)


# This version is used by the original etabw function. Below is a derivative 
# for working with rasters
#      ****************************************************************
#      *                                                              *
#      *  PURPOSE:  To Look Up a Y value as a Function of X along a   *
#      *            Curve Defined By a Series of X,Y Points Connected *
#      *            By Straight Line Segments                         *
#      *                                                              *
#      *  METHOD:   Simple Linear Interpolation Between Coordinate    *
#      *            Pairs Whose X values Bound the Calling Argument   *
#      *            If X is Less Than the First X-coordinate in the   *
#      *            Curve Array, Then the First Y-coordinate is       *
#      *            Returned.  If X is Greater Than the Last          *
#      *            X-coordinate, Then the Last Y-coordinate is       *
#      *            Returned.                                         *
#      *                                                              *
#      *  CALL SEQUENCE:                                              *
#      *            VAR = CURVE (XL, DUMY, NPTS)                      *
#      *                                                              *
#      *  CALLING ARGUMENTS:                                          *
#      *            XL    = Lookup Value of X                         *
#      *            DUMY  = Two-Dimensional Array Defining the        *
#      *                     Curve                                    *
#      *            NPTS  = Number of Points Defining the Curve       *
#      *                                                              *
#      *  NOTE:                                                       *
#      *            Please pardon the antiquated IF statements.  The  *
#      *            routine is over 20 years old, but it works!       *
#      *                                                              *
#      ****************************************************************
#
# NOTE: instead of a 2 dimensional array, I just split DUMY into x and y vectors
# Probably not the cleverest way to do it, but it seemed intuitive
def CURVE_function(XL, DUMY_x, DUMY_y, NPTS):

    # The fortran, below, seems unnecessarily complicated to me (but maybe it's faster that way?)
    # I tried to get teh same thing using Python
    # initialize curve at 0
    CURVE = 0

    # Check if it falls off the left
    if XL < DUMY_x[0]:
        CURVE = DUMY_y[0]

    # Check if it falls off the right #**# NEED TO GET RID OF THE TRAILING 0's for this version
    elif XL > DUMY_x[-1]:
        CURVE = DUMY_y[-1]
    
    # Otherwise, interpolate the value
    else:
        CURVE = curve_interpolate(XL, DUMY_x, DUMY_y)

    return CURVE

# function to do the interpolation. #**# This could be optimized
def curve_interpolate(XL, DUMY_x, DUMY_y):

    # Loop through values of x, until one exceeds XL's value
    # Then calculate the interpolation with that value and the previous value
    for j in xrange(len(DUMY_x)):
        value = DUMY_x[j]
        if XL < value:
            i = j - 1
            # Get the difference between the y values, multiply by the  upper and lower bounds, multiply it by 
            CURVE = (DUMY_y[j] - DUMY_y[i]) * (XL - DUMY_x[i]) / (DUMY_x[j] - DUMY_x[i]) + DUMY_y[i]
            # End for loop when the if criterion is met            
            break
        
    return CURVE

# Raster version of the curve function
# First pass of writing this function aimed for getting code that works, not
# code optimized for speed. There are some serious inefficiencies in the code
# below, but I could not think of an easy way to overcome them.
# HOWEVER: the curve_raster_file vs. curve_raster distinction is not an
# inefficiency, it is a work-around for a bug in ArcGIS. One would assume that
# you could just pass a single curve_raster variable, and have it update with
# each iteration of the loop. This did not consistently work, hence the need
# to pass the file names and re-rasterize them on each iteration.
def CURVE_raster_function(XL_raster, DUMY_x, DUMY_y, NPTS, temp_path):

    # Check for raster values below the minimum x value. #the 0 if false is a place-holder
    #if XL < DUMY_x[0]:
    #    CURVE = DUMY_y[0]
    curve_raster_low = arcpy.sa.Con((XL_raster < DUMY_x[0]), DUMY_y[0], 0)
    curve_raster_low.save(temp_path + "curve_step_low.tif")

    # Check if raster values are larger than the maximum x value. The false value keeps the original raster value
    #elif XL > DUMY_x[-1]:
    #    CURVE = DUMY_y[-1]
    curve_raster_high_file = temp_path + "curve_step_high.tif"
    curve_raster_high = arcpy.sa.Con((XL_raster > DUMY_x[-1]), DUMY_y[-1], curve_raster_low)
    curve_raster_high.save(curve_raster_high_file)
    
    # Otherwise, interpolate the value
    #else:
    #    CURVE = curve_interpolate(XL, DUMY_x, DUMY_y)
    curve_raster_new_file = curve_raster_high_file
    for j in xrange(1, len(DUMY_x)): # Skip first entry - we checked that one above (0)
        curve_raster_old_file = curve_raster_new_file    # Reset the curve raster file each time to try to bypass an ArcGIS glitch    
        k = j - 1
        new_value_raster = ((DUMY_y[j] - DUMY_y[k]) * (XL_raster - DUMY_x[k]) / (DUMY_x[j] - DUMY_x[k]) + DUMY_y[k])   
        new_value_raster.save(temp_path + "nvr_%s.tif" % j)
        
        # This SQL expression looks valid to me, but does not execute properly - it gives true as long as the first criterion is met and disregards the second criterion!
        # This is true when I enter it into ArcGIS, or if I remove the '' around VALUE, or if I change the order
        #expression = "'VALUE' < %s AND 'VALUE' >= %s" % (DUMY_x[j], DUMY_x[k])
        #curve_raster = arcpy.sa.Con(XL_raster, new_value_raster, curve_raster, expression)
        
        # This also did not work:
        #curve_raster = arcpy.sa.Con((XL_raster >= DUMY_x[k] and XL_raster < DUMY_x[j]), new_value_raster, curve_raster)

        # Test each condition separately, add them to find where they are both true, and then use that as the condition
        con1 = arcpy.sa.Con(XL_raster < DUMY_x[j], 1, 0)
        con1.save(temp_path + "con1_%s.tif" % (j + 2))
        con2 = arcpy.sa.Con(XL_raster >= DUMY_x[k], 1, 0)
        con2.save(temp_path + "con2_%s.tif" % (j + 2))
        con3 = con1 + con2
        con3.save(temp_path + "con3_%s.tif" % (j + 2))
        curve_raster = arcpy.Raster(curve_raster_old_file)
        curve_raster_new_file = temp_path + "curve_step%s.tif" % (j + 2)
        curve_raster_new = arcpy.sa.Con(con3 == 2, new_value_raster, curve_raster)
        curve_raster_new.save(curve_raster_new_file)

    return curve_raster_new_file


