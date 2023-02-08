      subroutine NMSimSource(dbband,refdist,distance,xyzsrc,xyzrec,zsrc,
     &                       zrec,heading,roll,vel,pitch,engpow,srcfile,
     &                       len_str,nfreq)

      !DEC$ ATTRIBUTES DLLEXPORT::NMSimSource
      !!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, DECORATE, ALIAS:"NMSIMSOURCE" :: NMSimSource
     
!     This will provide the output for a NMSim style source definition.
!     The input parameters will be generic to allow for future more complex source definitions
!
!  Written by Bruce J. Ikelheimer
!  Blue Ridge Research and Consulting
!  6/24/2009
!
!  29 January 2016
!  Moodifications to allow direct acces from Python
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Input Parameters
!
!  xyzsrc    - an array that contains the x,y,z position of the source in meters, with source height in AGL
!  xyzrec    - an array that contains the x,y,z position of the receiver in meters, with source height in AGL
!  zsrc      - the elevation of the ground under the source in meters MSL
!  zrec      - the elevation of the ground under the receiver in meters MSL
!  heading   - the true-north heading for the vehilce - not necessarily its direction of travel, in degrees
!  roll      - the roll angle in degrees
!  pitch     - the climb angle in degrees
!  vel       - the velocity of the vehicle in m/s
!  engpow    - the engine power setting in % of max.  Can be greater than %100
!  srcfile   - source file name with full path
!  nfreq     - number of frequencies
!  reset     - flag asking the program to reset for a new source.
!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Output Parameters
!
!  dbband    - output of the nfreq 1/3 octave band for this source
!  refdist   - reference distance for this source
!  distance  - Distance between the source and the receiver
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
c     
      Implicit none
      
c     Input
      real,intent(in) :: xyzsrc(3),xyzrec(3),zsrc,zrec
      real,intent(in) :: heading,roll,vel,pitch,engpow
      character*(len_str),intent(in) :: srcfile
      integer,intent(in) :: nfreq ,len_str

c     Output
      real,intent(inout) :: dbband(nfreq)
      real,intent(inout) :: refdist,distance
      
c     Internal
      integer nfreq_tmp,npower,n,m,j,i,npow
      integer ithet,ibnd
      real thrad,phi,elvang,thdeg,db1,db2
      real fracth,fmet,gamma,fracp
      real power(20)
	real bands(nfreq,0:36,20)
      character*255 srcpath,lastpath,src
	character*15 fmtout
      character*25 infile
	character*5 srcinfo

      refdist=1000.0/3.2808  !NMsim sources have a reference distance of 1000 feet.
      dbband=-99.9

      write(fmtout,'(''(20x,'',i2,''f5.1)'')') nfreq
      srcpath=srcfile(1:index(srcfile,'\',back=.true.)-1)
      
      call chdir(srcpath)

      src=srcfile(index(srcfile,'\',back=.true.)+1:
     &            len_trim(srcfile))
      
      open(unit=17,file=src,shared)
	read(17,'(a5)') srcinfo
	read(17,*)
	read(17,*)
	read(17,*)
      read(17,*)nfreq_tmp
      if(nfreq.ne.nfreq_tmp) return !An error with the number of frequencies
	read(17,*)npower
      do n = 1,npower
        read(17,*)infile,power(n)
        open(unit=18,file=infile,shared)
	  read(18,*)
	  read(18,*)
        read(18,fmtout)
     &       ((bands(m,j,n),m=1,nfreq),j=0,36)                   
        close(18)
      enddo
      close(17)
      call chdir(lastpath)

c Convert input theta to degrees, find the index of the beginning of the
c interval in which it falls.
c
      call incidence(xyzsrc,xyzrec,zsrc,zrec,heading,pitch,roll,
     &               distance,thrad,phi,elvang)
      
      thdeg = thrad*57.2957795                 
      ithet = thdeg/5.
      ithet = min(ithet,35)
      fracth = (thdeg/5.) - float(ithet)
c
c Find the power interval in which the input engpow falls
c
      if(srcinfo.ne.'Autom' .and. srcinfo.ne.'Train') then
	  fmet=engpow
	else if (srcinfo.eq.'Autom') then
	  fmet=vel*2.369
	else 
	  fmet=gamma*180.0/3.14159
	end if
c
      if(npower.gt.1) then
        do n = 1,npower-1
          if(fmet.lt.power(n+1)) exit
        enddo
        npow = min(n,npower-1)
        fracp = (fmet-power(npow))/(power(npow+1)-power(npow))
c Limit things to the range of input powers
c        fracp = max(fracp,0.)
c        fracp = min(fracp,1.)
	else
	  fracp=power(1)
	  npow=1
	end if
c
c Two-way linear interpolation on Theta and engpow.  DB1 and DB2 are values
c at Theta on the bracketing power valuse.  These are then interpolsted
c on engpow.
c
      do 50, ibnd = 1,nfreq
      db1 = bands(ibnd,ithet,npow)*(1.-fracth) +
     +      bands(ibnd,ithet+1,npow)*fracth
      if(npower.gt.1) then
        db2 = bands(ibnd,ithet,npow+1)*(1.-fracth) +
     +      bands(ibnd,ithet+1,npow+1)*fracth
	else
	  db2=db1
	end if

c
      dbband(ibnd) = db1*(1.-fracp) + db2*fracp

c     Limit the extrapolation to 5 dB

      if(dbband(ibnd).gt.min(bands(ibnd,ithet,1),
     &                       bands(ibnd,ithet+1,1))-5.0) 
     &   dbband(ibnd)=max(dbband(ibnd),min(bands(ibnd,ithet,1),
     &                    bands(ibnd,ithet+1,1))-5.0)

      if(dbband(ibnd).gt.max(bands(ibnd,ithet,npower),
     &                       bands(ibnd,ithet+1,npower))+5.0) 
     &   dbband(ibnd)=min(dbband(ibnd),max(bands(ibnd,ithet,npower),
     &                    bands(ibnd,ithet+1,npower))+5.0)
      
 50   continue
            
      end     