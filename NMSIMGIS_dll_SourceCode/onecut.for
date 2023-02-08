      subroutine onecut(attbnd,npts,prof_dist,hgt,imp,hsrc,hrec,
     1                  nfreq,freq)

      !DEC$ ATTRIBUTES DLLEXPORT::onecut

c This routine performs propagation loss calculations along a single
c terrain cut.  It is derived from the ONECUT routine developed by
c J. Page and K. Plotkin for Plovsing simplified algorithms, and further
c developed by J. Czech to include the full Rasmussen algorithms.
c
c This version varies from the Fall 1995 Czech version (as used in CHKPROP
c and the corresponding developmental version of NMAP) in that the profile
c has already been extracted.
c
c The output is the attenuation and the type of modeled terrain,
c obtained by:
c   a) specifying the source and receiver coordinates
c   b) specifying the source model P (Plovsing algorithms) or R (Ras-
c      mussen algorithms).
c
c If the user chooses P, a single a-weighted ATTEN value is returned.
c If the user chooses R, the full third-octave attenuation spectrum
c ATTBND is returned.
c
c 28 January 2016
c Modifications to OneCut to allow access by Python, and removed Plovsing routines
c
c ******************* VARIABLE DICTIONARY *********************************
c
c  Arguments
c              i/o
c attbnd(24)    o  band-by-band terrain attenuation

c prof(2,npts)  i  array of x and z data for the cross-section of terrain
c                  between the source & receiver
c                  1,* = x (radii from source) in ft
c                  2,* = z (elevation) in ft, MSL
c npts          i  number of points in the elevationprofile PROF
c flowpr(2,npts) i  Array of radii and flow resistivity from the source
c                  to the receiver
c hrec          i  height of receiver AGL, feet
c hsrc          i  height of source AGL, feet
c nfreq         i  number of frequency bands to compute

c
c  Local variables
c
c alp1         angle of plane of source
c alp2         angle between plane of source and plane of slope relative to
c              plane of source
c alpha        angle between plane of source and plane of slope
c ampf         Conversion multiplier from feet to meters
c atthrd       attenuation due to totally hard terrain
c attsft       attenuation due to totally soft terrain
c ax           dummy for source x-coordinate (ft)
c ay           dummy for source y-coordinate (ft)
c az           dummy for source elevation
c d0           distance parameter (m) [see subroutine calls]
c d1           distance parameter (m) [see subroutine calls]
c d2           distance parameter (m) [see subroutine calls]
c d3           distance parameter (m) [see subroutine calls]
c dist         dummy for distance between two points
c dsl          distance between two points (2D or 3D)
c dumf         dummy for frequency (Hz)
c duml         dummy for level (dB) relative to the free field with
c              turbulence [EGAL subroutine]
c hgtmm(2)     maximum and minimum heights in the profile as examined by
c              MAXMIN
c hgts(3)      heights of the three important points in the profile
c hillxz(2,5)  x,z of all five critical points.  The ends are in (i,1)
c              and (i,5), while the middle three are the three important
c              points at locs.done; i,3 is the top of the hill
c hrdum        dummy for height of receiver (ft or m, AGL)
c hsdum        dummy for height of source (ft or m, AGL)
c htr          screen height rel. to the ground @ receiver (m)
c hts          screen height rel. to the ground @ source (m)
c ignd(100)    dummy grid of terrain type; value of 0=soft; 1=hard
c igrnd        flag for terrain type; 0=soft; 1=hard
c isoft        integer counter of intervening ground type of soft
c ihard         o  integer counter of intervening ground type of hard
c k1           dummy for klocs(1)
c k2           dummy for klocs(2)
c k3           dummy for klocs(3)
c klocs(5)
c locmm(2)     indices of maximum and minimum heights of the profile as
c              examined by MAXMIN
c locs(3)      indices of the three important points in the profile;
c              (2) points to the top of the hill
c n            dummy for the number of points in the profile
c nlocs
c nmm          number of important points in the profile; 1=level;
c              2 or more=hill or valley
c ntemp
c ox           dummy for receiver x-coordinate (ft)
c oy           dummy for receiver y-coordinate (ft)
c oz           dummy for receiver elevation
c pi           Constant, value of pi=3.14....
c radius(100)  dummy array of flowpr(1,*) data
c radum        dummy for radius(i)
c SI(2)        array containing values for flow resistivity in kNs/m^4
c              (1) = soft ground
c              (2) = hard ground
c attbnd(24)   propagation loss by 1/3 octave band
c v            wind velocity at 4m AGL (m/s)
c vl           effective wind velocity returned by function EFFV
c z0           height of wind measurement near ground (m, AGL)
c zcrit        critical altitude at important point in profile
c zz1          altitude at 1st important point in profile
c zz2          altitude at 2nd important point in profile
c zz3          altitude at 3rd important point in profile
c *******************************************************************
c  Calls to:
c     VARYSURF subroutine -- mixed terrain (hard/soft) prop. model
c     EGAL subroutine -- level topography propagation model for
c                         spectral sound levels
c     BAKKE subroutine -- hill topography propagation model for
c                         spectral sound levels
c     DAL subroutine -- valley topography propagation model for
c                         spectral sound levels
c     MAXMIN function -- examines profile PROF and returns the heights
c                        at several critical locations, indices LOCS
c *******************************************************************
c
      implicit none

c     Output
      real,intent(inout) :: attbnd(nfreq) !This is the spectra attenuation

c     Input
      integer,intent(in) :: npts        !Number of points in the profile
      real,intent(in) :: prof_dist(npts)     !Array of distances
      real,intent(in) :: hgt(npts)      !Array of heights
      real,intent(in) :: imp(npts)      !Array of surface impedances
      real,intent(in) :: hsrc, hrec     !Source and receiver heights above the terrain
      integer,intent(in) :: nfreq       !The number of frequencys to analyize
      real,intent(in) ::  freq(nfreq)   !Array of one-third octave band center frequencies

c     Internal variables     

c Arrays for profile max/min data and inputs to hill model
      real prof(2,npts)   !The terrain profile in meters
      real flowpr(2,npts) !The impedance profile
      real hgts(3)
      integer locs(3)
      real hgtmm(2),locmm(2),klocs(5),hillxz(2,5),hillm(2,5)
c
c Source/receiver actual/mirror locations.  First index is x,z, second
c is actual/mirror
c
      real srcloc(2,2),recloc(2,2)
      real radius(npts),radum(npts),SI(2)
      real ampf,pi,vl,floval,flosoft,flohard
      real dx,dy,ax,ay,az,ox,oy,oz
      real hsdum,hrdum,dist,atthrd,atten,attsft,d1,d2
      real dumf,duml,k1,k2,k3
      real dsl,zz1,zz2,zz3,zcrit
      real costh,sinth,delx,thet
      real hsperp,hrperp,disdum,dsldum,alp1,alp2,alpha,d0

      integer ignd(npts)
      integer imdltyp,ihard,isoft,iflosoft,i,nlocs,ntemp,nmm
      integer igrnd,n,j

      data ampf /.30480/
      data pi /3.141592654/
c
c Hard-wire effective wind to be zero.  Corresponding code (which calls
c EFFV) is commented out.
      data vl/0./
c
c Flow resistivity threshold for soft/hard ground
      data floval/1000./
c
c
c ======== Section to Process Profile ===========
c
c Transfer data from flow resistivity profile to RADIUS and IGND.
c Flow resitivity can have any number of values, but the algorithms
c can handle only two: hard or soft.  A threshold of 1000 rayls divides
c hard from soft.  All flow values below that are considered to be soft,
c and all above are hard.  The program will use the average of all soft
c values in the current cut for soft, and the average of all hard values
c in the current cut for hard.
c
c Note that this bi-modal form is local along each cut, so the program
c accomodates more than just two values of flow resistivity.
c
c In addition to copying the data and determining the average hard/soft
c values for this cut, a count of hard/soft is made.  Calculations are done
c by doing all-hard and all-soft analysis, and proprtioning via VARYSURF.
c
c     Step one is to copy the data into the prof and flowpr arrays
      do i=1,npts
        prof(1,i)   = prof_dist(i)
        flowpr(1,i) = prof_dist(i)
        prof(2,i)   = hgt(i)
        flowpr(2,i) = imp(i)
      end do
c
      ihard = 0
      isoft = 0
      floSoft = 0.0
      floHard = 0.0

      do 18, i = 1,npts
         radius(i) = flowpr(1,i)
         ignd(i) = 1
         if(flowpr(2,i).le.floval) ignd(i) = 0
         if(ignd(i).eq.0) then
             isoft = isoft+1
             floSoft = floSoft + flowpr(2,i)
           else
             ihard = ihard+1
             floHard = floHard + flowpr(2,i)
          endif
 18   continue
      floSoft = floSoft/(max(isoft,1))
      floHard = floHard/(max(ihard,1))
c
c If either comes out to be zero (i.e., no values in one of the ranges)
c assign default value.  That's just a filler to avoid calling the ground
c effects routine with zero flow resistivity.
c
      if(floSoft.eq.0.0) floSoft=200.
      if(floHard.eq.0.0) floHard=1000000.
c
c Insert these values into array SI, which is passed to the ground effects
c routines.
c
      si(1) = floSoft
      si(2) = floHard
c
c=======================================================================
c
c The following routine returns the max and min heights and their indices,
c as the first and second values in hgtmm and locmm.  For this pass, we're
c interested in the highest point.
c
      call maxmin(prof,npts,hsrc,hrec,hgtmm,locmm)
c
c Profiles are characterized by three important points: the highest, and the
c lowest to either side, plus the end points.  This set of five can represent
c from two to five distinct points.
c
c Get the three important points.
c
c Store the max as the second important point, then get the other two
c points.
c
      hgts(2) = hgtmm(1)
      locs(2) = locmm(1)
      nlocs = 3
c     The first important point is the minimum in the near-source half
      call maxmin(prof,locs(2),hsrc,0.,hgtmm,locmm)
      hgts(1) = hgtmm(2)
      locs(1) = locmm(2)
c     The third important point is the minimum in the near-receiver half
      ntemp = npts - locs(2) + 1
      call maxmin(prof(1,locs(2)),ntemp,0.,hrec,hgtmm,locmm)
      hgts(3) = hgtmm(2)
      locs(3) = locs(2) + locmm(2) -1
c
c Set up an array with 2 to 5 points.  The set is source end, the three
c locs, the receiver end.  Redundant points are deleted.  Also set up an
c array with all 5, for input to the hill model.
c
      hillxz(1,1) = prof(1,1)
      hillxz(2,1) = prof(2,1)
      klocs(1) = 1
      nmm = 1
      do 15, i = 1,3
         hillxz(1,i+1) = prof(1,locs(i))
         hillxz(2,i+1) = prof(2,locs(i))
         if(locs(i).eq.klocs(nmm)) go to 15
         nmm = nmm+1
         klocs(nmm) = locs(i)
 15   continue
c
      hillxz(1,5) = prof(1,npts)
      hillxz(2,5) = prof(2,npts)
      if(klocs(nmm).lt.npts) then
         nmm = nmm+1
         klocs(nmm) = npts
      endif
c
c hillxz(j,3) is the top of the hill.  LOCS(2) points to it.  If either
c LOCS(1) or LOCS(3) point to LOCS(2), make the corresponding hillxz point
c one of the ends.  (Need to make sure we're not trampling a flat or uphill
c case here.)
c
      if(locs(1).eq.locs(2)) then
         hillxz(1,2) = hillxz(1,1)
         hillxz(2,2) = hillxz(2,1)
      endif
c
      if(locs(3).eq.locs(2)) then
         hillxz(1,4) = hillxz(1,5)
         hillxz(2,4) = hillxz(2,5)
      endif

c Correct the hillxz to eliminate a case without a flat at either end.
c This is a correction made 8/1/03.

      if(hillxz(1,1).eq.hillxz(1,2)) then
        dx = hillxz(1,3)-hillxz(1,1)
        dy = hillxz(2,3)-hillxz(2,1)
        hillxz(1,2) = hillxz(1,1) + 0.1*dx
        hillxz(2,2) = hillxz(2,1) + 0.1*dy
      end if

      if(hillxz(1,4).eq.hillxz(1,5)) then
        dx = hillxz(1,5)-hillxz(1,3)
        dy = hillxz(2,5)-hillxz(2,3)
        hillxz(1,4) = hillxz(1,3) + 0.1*dx
        hillxz(2,4) = hillxz(2,3) + 0.1*dy
      end if

c
c =========== Profile and Supporting Stuff Established ========
c =========== Go to Various Propagation Models =============
c Two points in the model means it's level.  More and we skip it.
c
      if(nmm .gt. 2) go to 100
c
c -----------------Level Model-----------------
c
c Input parameters for LEVEL
      ax = prof(1,1)
      ay = 0.0
      az = prof(2,1)
      ox = prof(1,npts)
      oy = 0.0
      oz = prof(2,npts)
c
c Copy geometry into metric working storage
c
      hsdum = hsrc  !*ampf  3/29/02  Everything should already be metric
      hrdum = hrec  !*ampf
      dist = sqrt((ax-ox)**2 + (ay-oy)**2)  !*ampf
c
c Rasmussen model for level ground called EGAL which requires:
c      1. distance from source over which the flow resistivity
c         near the source applies (m)
c      2. distance from receiver over which the flow resistivity
c         near the receiver applies (m)
c      3. source height (m, AGL)
c      4. receiver height (m, AGL)
c      5. flow resistivity near the source (kNs/m^4)
c      6. flow resistivity near the receiver (kNs/m^4)
c      7. effective wind velocity (m/s) at 4m AGL [VL]
c      8. Altitude at which the air flow transitions to constant
c         sound speed (m, AGL) [APH].  Must be >=10m for non-zero wind
c         velocity
c      9. height of wind measurement near ground (m, AGL) [Z0].
c         Specify as negative to request linear wind profile.
c     10. turbulence paramter (m/s) [MN].  Must be between 1 and 2.
c     11. Frequency of interest (Hz) [DUMF]
c
c EGAL does its work in double precision but returns single precision.
c EGAL returns the level relative to the free field with turbulence
c and without turbulence; the former being ignored with the use of
c the variable duml.
c
c Although EGAL allows for a wind velocity profile and two values
c of flow resistivity, we employ it with no wind and and a single
c flow resistivity.  Winds are not addressed, and VARYSURF will be
c used for mixed terrain.
c
c The distance from the source to the receiver is modelded as two-
c dimensional in x and y, just like inputs to LEVEL.
c
         D1=dist/2.
         D2=D1
         do 50 j=1,nfreq
            dumf=freq(j)
c
            if(ihard.gt.0) then
               igrnd = 1
               call EGAL(D1,D2,hsdum,hrdum,SI(igrnd+1),
     1              SI(igrnd+1),VL,0.,0.,dumf,duml,atthrd)
               atten = atthrd
            endif
c
            if(isoft.gt.0) then
               igrnd = 0
               call EGAL(D1,D2,hsdum,hrdum,SI(igrnd+1),
     1              SI(igrnd+1),VL,0.,0.,dumf,duml,attsft)
               atten = attsft
            endif
c
            if(ihard.gt.0 .and. isoft.gt.0) then
               n = npts
               do 35, i = 1,n
                  radum(i) = radius(i)  !*ampf
  35           continue
               call varysurf(radum,ignd,n,hsdum,hrdum,
     1                       attsft,atthrd,atten)
            endif
            attbnd(j)=atten
  50     continue

c
c  Label 100 is the target of the not-level skipover.
 100  continue
c
c
c Hill only if there are at least three points in the profile and the
c middle one is higher.
c
c --------------  Hill Model ---------------
c
c Selection of this is via the following IF.  ZCRIT is the altitude, at the
c middle point, of a line connecting the end points.  Let the "equals" case
c be included as a hill.  A 4 or 5 point case will always be a hill.  A 3
c point case that's exactly on-line would never occur: the lines would have
c to be parallel, in which case MAXMIN would select an end, not the middle.
c
      if(nmm.eq.3) then
         k1 = klocs(1)
         k2 = klocs(2)
         k3 = klocs(3)
c      Total distance dist and distance to second point dsl
         dist = prof(1,k3) - prof(1,k1)
         dsl = prof(1,k2) - prof(1,k1)
c      Altitudes at the three points
         zz1 = prof(2,k1)
         zz2 = prof(2,k2)
         zz3 = prof(2,k3)
         zcrit = zz1 + (zz3-zz1)*dsl/dist
      endif
c
      if(nmm.gt.3 .or. (nmm.eq.3 .and. zz2.ge.zcrit)) then
c
c Copy geometry into metric working storage
c
            hsdum = hsrc  !*ampf
            hrdum = hrec  !*ampf
            do 140, i = 1,2
               do 130, j = 1,5
                  hillm(i,j) = hillxz(i,j)  !*ampf
 130           continue
 140        continue
c
c Rasmussen model for hills called BAKKER which requires:
c      1. distance from source to beginning of wedge (m)
c      2. distance from beginning of wedge to apex (m)
c      3. distance from apex to end of wedge (m)
c      4. distance from end of wedge to receiver (m)
c      5. source height (m, AGL)
c      6. receiver height (m, AGL)
c      7. screen height rel. to the ground @ source (m)
c      8. screen height rel. to the ground @ receiver (m)
c      9. flow resistivity between source and wedge (kNs/m^4)
c     10. flow resistivity for wedge (kNs/m^4)
c     11. flow resistivity between wedge and receiver (kNs/m^4)
c     12. frequency of interest (Hz)
c VARYSURF will account for mixed terrain surfaces.
c BAKKE returns the level relative to the free field.
c
c Set up the source and receiver locations, normal to the corresponding
c plateaus.
            costh = 1.
            sinth = 0.
            delx = hillm(1,2) - hillm(1,1)
            if(delx.gt.0.) then
                thet = atan2(hillm(2,2)-hillm(2,1),delx)
                costh = cos(thet)
                sinth = sin(thet)
             endif
c
ckjp July 04.  When this was first written, slopes were considered to be
c modest and, and the original source and receiver locations were moved
c such that they were on lines perpendicular to the end flats, those lines
c passing through the end points of hillm.  That caused problems when the
c ground slopes were not small.  The correct geometry is to keep the
c source and receiver positions as is, and to make the images on reflections
c from those.
c
c           srcloc(1,1) = hillm(1,1) - sinth*hsdum
c           srcloc(2,1) = hillm(2,1) + costh*hsdum
c           srcloc(1,2) = hillm(1,1) + sinth*hsdum
c           srcloc(2,2) = hillm(2,1) - costh*hsdum
c
c Source location is hs above the start of the terrain cut
            srcloc(1,1) = hillm(1,1)
            srcloc(2,1) = hillm(2,1) + hsdum
c New "reflect the original" source image
            hsperp = hsdum*costh
            srcloc(1,2) = srcloc(1,1) + 2*sinth*hsperp
            srcloc(2,2) = srcloc(2,1) - 2*costh*hsperp
c
            costh = 1.
            sinth = 0.
            delx = hillm(1,5) - hillm(1,4)
            if(delx.gt.0.) then
                thet = atan2(hillm(2,5)-hillm(2,4),delx)
                costh = cos(thet)
                sinth = sin(thet)
             endif
c
c Old "move them both" receiver and receiver image
c           recloc(1,1) = hillm(1,5) - sinth*hrdum
c           recloc(2,1) = hillm(2,5) + costh*hrdum
c           recloc(1,2) = hillm(1,5) + sinth*hrdum
c           recloc(2,2) = hillm(2,5) - costh*hrdum
c
c New "right over the end of the cut" receiver
            recloc(1,1) = hillm(1,5)
            recloc(2,1) = hillm(2,5) + hrdum
c New "reflect the original" receiver image
            hrperp = hrdum*costh
            recloc(1,2) = recloc(1,1) + 2*sinth*hrperp
            recloc(2,2) = recloc(2,1) - 2*costh*hrperp
c
            do 150 j=1,nfreq
               dumf=freq(j)
c
               if(ihard.gt.0) then
                  igrnd = 1
                  call BAKKER(hsdum,hrdum,hillm,srcloc,recloc,
     1                 SI(igrnd+1),SI(igrnd+1),SI(igrnd+1),dumf,atthrd)
                  atten = atthrd
               endif
c
               if(isoft.gt.0) then
                  igrnd = 0
                  call BAKKER(hsdum,hrdum,hillm,srcloc,recloc,
     1                 SI(igrnd+1),SI(igrnd+1),SI(igrnd+1),dumf,attsft)
                  atten = attsft
               endif
c
               if(ihard.gt.0 .and. isoft.gt.0) then
                  n = npts
                  do 135, i = 1,n
                    radum(i) = radius(i)  !*ampf
 135              continue
                  call varysurf(radum,ignd,n,hsdum,hrdum,
     1                          attsft,atthrd,atten)
               endif
            attbnd(j)=atten
 150        continue
         endif
c
c--------------------- Valley Model -------------------------
c
c Valley model only if there are exactly three points in the profile and
c the middle one is lower than the ends.  This is a condition excluded
c by the test for the hill model, so only one or the other will apply.
c
      if(nmm.eq.3 .and. zz2.lt.zcrit) then
c
c Copy geometry into metric working storage and compute angles
c
         hsdum = hsrc  !*ampf
         hrdum = hrec  !*ampf
         disdum = dist !*ampf
         dsldum = dsl  !*ampf
         alp1 = atan2((prof(2,k2) - prof(2,k1)),
     1                (prof(1,k2) - prof(1,k1)))
         alp2 = atan2((prof(2,k3) - prof(2,k2)),
     1                (prof(1,k3) - prof(1,k2)))
         alpha = (alp2-alp1)

c Rasmussen model for valley which requires:
c      1. distance from source to beginning of slope (m)
c      2. distance from beginning of slope to the receiver (m)
c      3. source height (m, AGL)
c      4. receiver height (m, AGL)
c      5. up-slope angle (radians) computed CCW from plane of source
c      6. flow resistivity between source and slope (kNs/m^4)
c      7. flow resistivity between begin. of slope and receiver
c         (kNs/m^4)
c      8. frequency of interest (Hz).
c VARYSURF will account for mixed terrain surfaces.
c DAL returns the level relative to the free field.
            alpha=alpha+pi
            D0=sqrt((prof(1,k2)-prof(1,k1))**2+
     1              (prof(2,k2)-prof(2,k1))**2)  !*ampf
            D1=sqrt((prof(1,k3)-prof(1,k2))**2+
     1              (prof(2,k3)-prof(2,k2))**2)  !*ampf
            do 250 j=1,nfreq
               dumf=freq(j)
c
               if(ihard.gt.0) then
                  igrnd = 1
                  call DAL(D0,D1,hsdum,hrdum,alpha,SI(igrnd+1),
     1                     SI(igrnd+1),dumf,atthrd)
                  atten = atthrd
               endif
c
               if(isoft.gt.0) then
                  igrnd = 0
                  call DAL(D0,D1,hsdum,hrdum,alpha,SI(igrnd+1),
     1                     SI(igrnd+1),dumf,attsft)
                  atten = attsft
               endif
c
               if(isoft.gt.0 .and. ihard.gt.0) then
                  n = npts
                  do 235, i = 1,n
                     radum(i) = radius(i)  !*ampf
 235               continue
                  call varysurf(radum,ignd,n,hsdum,hrdum,
     1                          attsft,atthrd,atten)
               endif
               attbnd(j)=atten
 250        continue
         endif

      return
      end
