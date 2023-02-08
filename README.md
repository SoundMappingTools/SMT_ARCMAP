# SMT_ARCMAP
Sound Mapping Tools Code for ArcMap

Toolbox: Sound Mapping Tools
Title: Geographic Information System for Prediction of Acoustic Detectability
Version: 4.4.2
Date: 2017-06-22
Authors: Alexander "Sasha" Keyel and Sarah Reed
         NMSim code contributed by Bruce Ikelheimer, see source code for
         original authors. Jessica Sushinski updated the toolbox to
         ArcGIS 10.X
	 Some code based on code originally written by Damon Joyce, National Park Service.
Maintainer: Sasha Keyel <skeyel@gmail.com>
License: GPL 2
Depends:
    ArcGIS >=10.3
    ArcGIS Spatial Analyst Extension
    Python 2.7

DESCRIPTION: 
Sound Mapping Tools are tools "for modeling spatial patterns of anthropogenic
noise propagation in natural ecosystems. SPreAD-GIS incorporates commonly
available datasets on land cover, topography, and weather conditions to
calculate noise propagation patterns and excess noise above ambient conditions
for one-third octave frequency bands around one or multiple sound sources.
User-specified noise source characteristics, ambient sound conditions, and
frequency-weighting make SPreAD-GIS flexible to incorporate field measurements
and model noise propagation for any type of source, environment, or species.
SPreAD-GIS is a free, open-source application written in Python and implemented
as a toolbox in ArcGIS software."

Sound Mapping Tools also includes NMSIMGIS and a GIS implementation of
ISO 9613-2. All three models produce different predictions, consequently
careful model choice is warrented.
- quote from Reed et al. 2012, see citation below.

CITATION:
Reed, S.E., Boggs, J.L., Mann, J.P. 2012. A GIS tool for modeling anthropogenic 
    noise propagation in natural ecosystems. Environmental Modelling & Software
    37: 1-5

Keyel, A.C. and Reed, S.E. 2017. Sound Mapping Tools: an ArcGIS toolbox for
    modeling the propagation of sounds in a wildland setting. Version 4.4.
    Colorado State University, Fort Collins, CO.

TEST ENVIRONMENTS (Windows 7, Python 2.7 for <=10.3, 2.8 for 10.4)
ArcGIS 10.1: THIS VERSION NOT TESTED
ArcGIS 10.2: THIS VERSION NOT TESTED
ArcGIS 10.3.0: PASSED 06/22/2017)
ArcGIS 10.4.1: PASSED 06/22/2017)
ArcGIS 10.5: NOT TESTED (LAST TEST V4.4.0 PASSED 03/29/2017 EXCEPT SEE NOTE 1)
ArcGIS Pro: NOT COMPATIBLE WITH ARCGIS PRO

NOTE 1: Using the old barrier settings with ArcGIS 10.5 led to a 20 dB
        difference in prediction for one test point due to a bug in the
        pathdistance calculation in ArcGIS 10.5. It is not recommended to use
        ArcGIS 10.5 with this setting (use_old_barrier = 1, default for
        SPreAD-GIS to maintain backwards compatibility)
NOTE 2: Slightly different model results between 10.3, 10.4, and 10.5

ISSUES AND BUGS ARE DESCRIBED IN Details.docx

Major NEWS: (more details given in Details.docx under Development History)
* Four stars (****) indicate changes that will lead to different model results
* Two stars (**) indicate changes that may change model results.

VERSION 4.4.2
CHANGES
  * Added 10 log10(d') as an output of the d' tool for audibility calculations
  * Added an option to obtain the maximum value across points instead of the sum
  * Default audibility threshold in AssessLengthAudible changed to 7.3 instead
    of 7 (based on Fidel, S., K. Pearsons, and M. Sneddon. Evaluation of the
    Effectiveness of SFAR 50–2 in Restoring Natural Quiet to Grand Canyon National
    Park. BBN Report 7197. No. 93-1. NPOA Report, 1994.)
  * Patched validation script to run correctly from the GUI for the ISO model test
  * Patched a bug where the 0 summarize option would crash the code
  * Minor patches to improve error handling of code
  * Patched this bug: "NLCD wetland codes 95, 96, 97, 98, 99 are not supported due to
    a bug. Please reclassify as NLCD code 90 (the wetland codes all map to the same
    impedance of 100,000)."
  * Added a Allow Defaults box to the NMSIMGIS and ISO 9613-2 dialogue boxes to
    improve awareness and control over the fields that control the source characteristics

VERSION 4.4.1 (Released 19 April 2017)
MAJOR CHANGES
  **** The NMSIMGIS .dll and the rounding approach used in bilinear interpolation
       were patched to fix a bug where sloped flat terrain was not giving
       consistent results with horizontal flat terrain

VERSION 4.4 (Released 12 April 2017)
  MAJOR ADDITIONS
  * Added a GIS implementation of ISO 9613-2 (spherical spreading, atmospheric
    absorption, ground effects, and barrier effects only)
  * Added a tool for calculating (human) audibility of sound given 1/3 octave
    levels and 1/3 octave background levels
  * Added a tool for calculating how loud a signal must be to be audible
    given background noise levels and species-specific critical ratio/band

  MAJOR CHANGES
  **** Changed the calculation of area for the AssessImpacts tool
  information
  **** Changed spherical spreading loss equation in SPreAD-GIS to correspond to
       ISO 9613-2
  **** Changed spherical spreading loss calculation for cell of origin for
    SPreAD-GIS to correspond to a 1 m distance (this may give different results
    than NMSIMGIS, which uses the difference between the origin cell coordinates
    and the source point coordinates)
  ** Changed wind speed of 0 to give 0 wind loss. Previously, 0 wind speed gave
    erroneously high levels of sound loss
  **** Added a new option for calculation of barrier loss which corrected many
    issues associated with the old calculation
  **** New option changed calculation of vegetation loss in the new option to be
    more in line with Harrison et al. 1980. Specifically, the 14 dB maximum loss
    is now included, and the vegetation loss is based on equations that better
    fit Table 8 from Harrison et al. 1980.        
  **** New option removes step where wind and vegetation are excluded when there
    is line of sight
  **** Moved calculate_mean_elevation and euclidean_dist_dir functions
    outside the frequency loop to try to avoid a schema lock and increase model
    efficiency. This appears to have slightly changed the model extent, and has
    led to very slight changes in model results under the old option (a change of
    0.21 dB for one test point).

VERSION 4.3 (patch of public release, 10/31/2016)
  MAJOR ADDITIONS
  * Added very limited support for landfire landcover to Add_SPREADTYPE and
    Add_NMSIMTYPE functions (contains syntax bugs that need to be fixed before it
    is functional, only a subset of the included landcover types are included).

  MAJOR CHANGES
  **** Patched Spherical Spreading Loss and Vegetation modules to eliminate
    NoData point at source. This may fix the the issue identified above where the
    maximum sound value does not appear at the sound source.
  **** Patched final smoothing step to exclude the cell of origin, as the
    smoothing changes greatly reduces the values relative to their true sound
    levels 

VERSION 4.2 ("soft" public release, 10/10/2016)
  MAJOR ADDITIONS
  * Re-organized the Toolbox and added an Extra Tools toolset to introduce
    additional useful tools without distracting users from the main model tools.

  MAJOR CHANGES
  **** Changed atmospheric absorption in nmsimgis to use three-dimensional
    distance instead of two-dimensional distance
  * Changed toolbox directory input to select the "toolbox" folder instead of the
    folder that contains it.

VERSION 4.1 (limited release 9/15/2016 at Technical Workshop)
  MAJOR ADDITIONS
  * nmsimgis results found to match main NMSim model (without Nord2000) and
    cleared for general use!

VERSION 4.0 (9/8/2016)
  MAJOR ADDITIONS
  * nmsimgis model now functional (but currently giving incorrect results)
  * Added background ambient comparison to nmsimgis
  * Added code to truncate overall results to 0 dB (toolbox, other options in
    Python)

  MAJOR CHANGES
  * Corrected a bug in basin extraction (changed to "NEAREST" from "BILINEAR")
  * Added ambient comparison to nmsimgis, SMT now supports comparison with
    overall background layers such as the Geospatial Sound Model outputs
  * Changed nmsimgis inputs (head, roll, pitch, velocity, engine power, source
    offset) to come from the input shapefile.
  * Changed SPreAD-GIS to no longer truncate individual 1/3 octave bands at 0 dB

VERSION 3.7 (8/13/2016, given to one person, not a stable release)
  PATCHES, BUT NO MAJOR CHANGES OR ADDITIONS

VERSION 3.6 (invitation only beta version, 6/3/2016)
  MAJOR ADDITIONS
  * Added a tool (AssessImpact) to assess sound impacts on a focal area
    (python only)

VERSION 3.5 (invitation only beta version, 5/13/2016)
  PATCHES, BUT NO MAJOR CHANGES OR ADDITIONS
 
VERSION 3.4 (Internal release 4/15/2016)
  MAJOR CHANGES
  * Toolbox name changed from SPreAD-GIS to Sound Mapping Tools to reflect
    addition of NMSimGIS
  * Tool interface changed to have separate tools for validation, spread-gis,
    and nmsim-gis due to a bug when trying to run a combined tool (only present
    when running script via toolbox directly in ArcGIS).
  * **** Patched bug (patch in all versions, bug only in 10.4). In 10.4,
    arcpy.sa.Sample produced erroneous -9999 values. The patch appears to have
    let to some slight changes in calculated intermediates but does not appear
    to have affected the overall results.
  **** Changed validation to only check aal raster to 4 decimal places to mask
    minor differences between 10.3 & 10.4
  * Source code documentation improved & flags (#**#) checked for inclusion in
    known issues

VERSION 3.3 (Never released)
  MAJOR CHANGES
  * Wind direction is now the direction the wind is blowing FROM not the
    direction the wind is blowing TO.
  * Simplified tool interface & updated ArcGIS tool documentation (in progress)
  * Changed Frequency, Sound Level of Source & Measurement distance to be input
    as "String" instead of "Double" to accommodate possibility of inputting
    tables.
  * Changed tool to allow selection of output directory
  * Model inputs are now metric instead of Old English (but some calculations
    still convert to Old English)
  * Changed to allow flexible choice of cell size (conversion for nmsimgis in
    progress)

VERSION 3.2
  MAJOR CHANGES  
  **** Changed viewshed tool to ObserverPoints tool to allow user to specify
    observer offset & starting elevation. This now includes a default 1 unit
    elevation offset & will produce different results than the 9.3 SPreAD-GIS
    model. This also patches a critical bug in Arc9.3 where unrealistic results
    were occuring in a simulated landscape.
  * Multifrequency code was modified to allow calculation of remaining
    frequencies when a frequency cannot be processed successfully
  * Added option to get timing of code sections for benchmarking purposes
  * Major changes to G_NMSim.py and nmsimhlpr.py scripts
  * A_AmbientSoundConditions.py can now be called as a function
  * (python only) Multipoint function changed to allow outputting results to a
    custom directory.
  ** Changed model extent to only use a minimum bounding rectangle and not
    custom shapes.

VERSION 3.1
  MAJOR CHANGES
  * Corrects bugs in Wind Loss, Topography and Barrier Loss, and overall
    summation sections
  * Corrected a bug in Wind Loss section where negative downwind values were
    being generated
  * Changed code to accept user-specified input paths
  * Removed old ArcGIS 9.3 code
  * Restricted allowed frequencies to those given in the drop-down menu in
    the ArcGIS toolbox
  * Added F_CodeValidation.py and validationhlpr.py to provide a worked example
    with the code and to provide test code to ensure consistent outputs.
  * Converted each Loss Section to a function and moved them to the helper
    script spreadhlpr
  * B_NoisePropagationforSinglePoint is now a function, with a small piece of
    code at the bottom of the script that only executes if it is the main
    function
  * NoisePropagationforSinglePoint.py (no prefix) was removed from the toolbox.
  * Scripts now start with letters instead of numbers to improve Python syntax
  * This README was added.

VERSION 3.0 (Never released)
  **** Atmospheric absorption differs between ArcGIS 9.3 and 10.x
     (but at 5 decimal places - very minor!)
  * Code updated to be compatible with ArcGIS 10.x.
     (Runs, but does not produce correct outputs)
