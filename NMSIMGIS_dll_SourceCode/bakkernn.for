      subroutine BAKKER(HS,HR,hill,srcloc,recloc,SI0,SI1,SI3,F,LEVEL)
c
c BAKKE as written assumes level plateaus on either side of wedge, and the
c D's and H's in the call adequately describe that.  Variation BAKKER (named
c after Jim & Tammy) handles the general case where nothing is on the level.
c The actual geometry, in HILL, is passed instead of the D's and H's.
c (HILL is a shorter, easier-to-type, name for HILLXZ.)
c
c Geometric calculations are essentially exact, using actual/reflected
c source/receiver locations set up in the matchng version of ONECUT.
c The exception is that mirror locations are taken to be at -HS and -HR
c relative to the ground points in an early calculation to get generic
c reflection coefficients for the flats.
c
c(* BEREGNING AF LYDUBREDELSE HEN OVER FORHINDRING I TERRíNET EFTER *)
c(* PRINCIPPER BESKREVET I LI RAPPORT 111 FRA 1984. DEN GEOMETRISKE *)
c(* SITUATION ER BESKREVET I FIGUR B1 OG B2 SIDE 53 *)
c
c(* KILE SVARER TIL JONW I RAP.111 P.59 MEN MED TO FLADE SEKTIONER*)
c(* I STEDET FOR AT HAVE Z SOM INDDATA HAR 'KILE' DOG SI
c (STRùMNINGSMODSTAND)*)
c(* SOM INDDATA*)
c(* VARIABELNAVNE FRA RAP. 111 ER STORT SET BIBEHOLDT*)
c(* REALDEL ER MARKERET MED ET EFTERHíNGT R OG IMAGINíRDEL*)
c(* MED ET EFTERHíNGT I *)
c
c*******************************************************************************
c
c     Input Variables
c       HS          Height of the source in meters
c       HR          Height of the receiver in meters
c       hill        x,z of all five critical points.  The ends are in (i,1)
c                     and (i,5), while the middle three are the three important
c       srcloc      Location of the source
c       recloc      Location of the receiver
c       SI0         Groud Flow Resistivity
c       SI1         Groud Flow Resistivity
c       SI3         Groud Flow Resistivity
c       F           Frequency of interest
c     Output:
c       Level       Resulting attenuation level for that frequency (dB)


      real hill(2,5),HS,HR,SI0,SI1,SI3,F
      dimension srcloc(2,2),recloc(2,2)
      integer I
      real A,HX,K,RR,R1,RD,F0,F1
      real AL1,T,T1,T2,TH,RH0,RH1,ANY,level
      real PI
      complex PO,PT,QS,QR,Q,Q1,Q2,QN,PL,AC,B,C,VD
      logical spew
c     spew = f.eq.400.
      spew = .false.
c
      PI = 4.*ATAN(1.)
      radeg = 180./pi
c
c        (*DELANY-BAZLEY BELOW SOURCE*)
      CALL DELBAZ(F,SI0,AC)
c
c        (*DELANY-BAZLEY BELOW BERM*)
      CALL DELBAZ(F,SI1,B)
c
c        (*DELANY-BAZLEY BELOW REC.*)
      CALL DELBAZ(F,SI3,C)
c
c          (* K IS WAVENO.*)
      K=2.0*PI*F/340.0
      PL=CMPLX(0.,0.)
      PT=CMPLX(0.,0.)
c
c          {MIRROR SOURCE}
ckjp Changed July 2004.  Original was based on the flat being near horizontal,
c    so vertical height from terrain point was good enough.  Current version
c    works relative to the angle of the flat.
c
c Distance from image source to top of hill
      delx = hill(1,3) - srcloc(1,2)
      dely = hill(2,3) - srcloc(2,2)
      rr = sqrt(delx**2 + dely**2)
c Angle from image to top of hill
      theti = atan2(dely,delx)
c Angle of the hillside
      dxh = hill(1,3) - hill(1,2)
      dyh = hill(2,3) - hill(2,2)
      theth = atan2(dyh,dxh)
c
c If the angle from the image is steeper than the slope of the hill, then
c this path is blocked.  If it's not steeper, then we compute spherical
c reflection based on angles and heights relative to the flat.
c
      if(theti.ge.theth) then
           qs = cmplx(0.,0.)
         else
c          Angle of the flat
           dxf = hill(1,2) - hill(1,1)
           dyf = hill(2,2) - hill(2,1)
           thetf = atan2(dyf,dxf)
c          Angle of the image path relative to the flat
           thet = theti-thetf
c          Net image source to receiver height, as needed by QQ
           delz = rr*sin(thet)
c          QQ needs horizontal distance along reflecting plane
c kjp 2 Sep 04 Take absolute distance, in case hill overhangs
c     the flat when in the flat's coordinates
           rr   = abs(rr*cos(thet))
           call QQ(RR,delz,K,AC,QS)
        endif
      if(spew) write(2,*)'delz,qs=',delz,qs
c
c          {MIRROR RECEIVER}
ckjp July 2004.  Idedntical logic to source.  Just swap point 5 for 1,
c    point 4 for 2, rec for src.  Also have X positive in the reverse
c    direction
c Distance from image receiver to top of hill
      delx = recloc(1,2) - hill(1,3)
      dely = hill(2,3) - recloc(2,2)
      rr = sqrt(delx**2 + dely**2)
c Angle from image to top of hill
      theti = atan2(dely,delx)
c Angle of the hillside
      dxh = hill(1,4) - hill(1,3)
      dyh = hill(2,3) - hill(2,4)
      theth = atan2(dyh,dxh)
c
      if(theti.ge.theth) then
           qr = cmplx(0.,0.)
         else
c          Angle of the flat
           dxf = hill(1,5) - hill(1,4)
           dyf = hill(2,4) - hill(2,5)
           thetf = atan2(dyf,dxf)
c          Angle of the image path relative to the flat
           thet = theti-thetf
c          Net image source to receiver height, as needed by QQ
           delz = rr*sin(thet)
c kjp 2 Sep 04 - overhang protection, as noted for source image
           rr   = abs(rr*cos(thet))
c          write(*,*)'theti,thetf,thet',theti,thetf,thet
c          write(*,*)'rr,delz,k',rr,delz,k
           call QQ(RR,delz,K,AC,QR)
        endif
      if(spew) write(2,*)'delz,qr=',delz,qr
c
c          {WEDGE ANGLE}
      delx = hill(1,3) - hill(1,2)
      delz = hill(2,3) - hill(2,2)
      T1=ATAN2(delx,delz)   !Angle of up slope
      delx = hill(1,4) - hill(1,3)
      delz = hill(2,3) - hill(2,4)
      T2=ATAN2(delx,delz)
      T=T1+T2
      if(spew) write(2,*)'T1,T2,T3=',t1*radeg,t2*radeg,t*radeg
c
c The 100 loop adds up diffracted contributions via four ray paths:
c     1. Direct source-vertex, direct vertex-receiver
c     2. Mirrored source-vertex, direct vertex-receiver
c     3. Direct source-vertex, mirrored vertex-receiver
c     4. Mirrored on both sides of vertex
c
      do 100 I=1,4
c         {mirror source and mirror receiver VIA WEDGE VERTEX}
c        HS=-ABS(HS)
c        IF((I.EQ.1).OR.(I.EQ.3)) HS=ABS(HS)
         isrc = 2
         IF((I.EQ.1).OR.(I.EQ.3)) isrc = 1
c
c        HR=-ABS(HR)
c        IF((I.EQ.1).OR.(I.EQ.2)) HR=ABS(HR)
         irec = 2
         IF((I.EQ.1).OR.(I.EQ.2)) irec = 1
c
         QN=CMPLX(1.,0.)
c
         IF(I.EQ.2) QN=QS
         IF(I.EQ.3) QN=QR
         IF(I.EQ.4) QN=QS*QR
c
c Get length and angle of path from source to top of wedge
         delz = hill(2,3) - srcloc(2,isrc)
         delx = hill(1,3) - srcloc(1,isrc)
c Distance
         rh0 = sqrt(delx**2 + delz**2)
c Angle, clockwise from straight down
         TH=ATAN2(delx,delz)   !Angle from source to vertex
c         {ANGLE FROM FIRST WEDGE-LEG TO SOURCE-RAY}
         F0=TH-T1
c
c Get length and angle of path from receiver to top of wedge.  Note the
c reversal of the X direction, so that the angle will be positive counter-
c clockwise.
         delz = hill(2,3) - recloc(2,irec)
         delx = recloc(1,irec) - hill(1,3)
         rh1 = sqrt(delx**2 + delz**2)
c Angle, counterclockwise from straight down
         TH=ATAN2(delx,delz)
c         {ANGLE FROM FIRST WEDGE-LEG TO RECEIVER-RAY}
         F1=2.*PI-T1-TH

c
c Save the F0 and F1 angles for later use.  I = 1 is both direct, I = 4
c is both reflected
c
         if(i.eq.1) then
              f0dir = f0
              f1dir = f1
           endif
         if(i.eq.4) then
              f0rfl = f0
              f1rfl = f1
           endif
c
c Total propagation path
         R1=RH0+RH1

c
           if(spew)write(2,*)'a.I,F0,F1=',i,f0*radeg,f1*radeg
         IF((F0.GT.0.).AND.((F1+T).LT.2.0*PI)) THEN
c              {HEIGHT OVER WEDGE-LEG}
            HX=SIN(F0)*R1

c              {FINITE WEDGE IMPEDANCE}
            call QQ(R1,HX,K,B,Q1)
c              {HEIGHT OVER WEDGE-LEG}
            HX=SIN(2.*PI-F1-T)*R1
c              {FINITE WEDGE IMPEDANCE}
            call QQ(R1,HX,K,B,Q2)
            if(spew) write(2,*)'Q1,Q2 =',q1,q2

            A=RH0*RH1/R1
            ANY=2.0-T/PI
            call DIF(R1,A,F1-F0,-1.0,ANY,K,VD)

            if(spew) write(2,*)'VD1 =',vd
            PL=VD
            call DIF(R1,A,F1+F0,-1.0,ANY,K,VD)
            if(spew) write(2,*)'VD2 =',vd
            VD=VD*Q1

c
c             {DIFFRACTION}
            PL=PL+VD
            call DIF(R1,A,F1+F0,1.0,ANY,K,VD)
            if(spew) write(2,*)'VD3 =',vd

            VD=VD*Q2
            PL=PL+VD
c
            call DIF(R1,A,F1-F0,1.0,ANY,K,VD)
            if(spew) write(2,*)'VD4 =',vd
            VD=VD*Q1*Q2
            PL=PL+VD
c
c             {MULTIPLY BY RELEVANT REFLECTION FACTOR}
            PL=PL*QN
            PT=PT+PL
            if(spew)write(2,*)'PT, |PT| =',pt,cabs(pt)
c             {END OF CONTRIBUTION LOOP}
         ENDIF
c
 100  continue
c                  {END OF DIFFRACTION CONTRIBUTIONS}
c
c The rest of this routine accounts for various direct and reflected
c ray paths over the vertex, should a line of sight exist between the
c source and receiver.  Angle AL1 measures how much of an upward bend
c a source-vertex-receiver path takes at the vertex.  Negative AL1 means
c a downward bend, i.e., the vertex sticks up into the path.
c
c First of three possible paths is direct from source to receiver.
c
      isrc = 1
      irec = 1
c
      f0 = f0dir
      f1 = f1dir
c
      AL1=PI+F0-F1
c
      if(spew)write(2,*)'b.I,F0,F1,AL1=',i,f0*radeg,f1*radeg,al1*radeg
c
      IF(AL1.GT.0.) THEN
c        {ADD DIRECT G.O. FIELD}
         delx = recloc(1,1) - srcloc(1,1)
         delz = recloc(2,1) - srcloc(2,1)
         RD=SQRT(delx**2.+delz**2.)
         PO=CEXP(CMPLX(0.,K*RD))/CMPLX(RD,0.)
         if(spew)write(2,*)'POb =',po
         PT=PT+PO
         if(spew)write(2,*)'PT, |PT| =',pt,cabs(pt)
      ENDIF
c
c Second possible path is from a mirrored source to the reciver.  Exactly
c the same formulae for AL1, except that -HS (the source mirrored about
c its ground plane) is used instead of HS.
c
      f0 = f0rfl
      f1 = f1dir
      AL1=PI+F0-F1
c
      if(spew)write(2,*)'c.I,F0,F1,AL1=',i,f0*radeg,f1*radeg,al1*radeg
      IF(AL1.GT.0.) THEN
c        {ADD REFLECTED G.O. FIELD OVER WEDGE}
         delx = recloc(1,1) - srcloc(1,2)
         delz = recloc(2,1) - srcloc(2,2)
         RD=SQRT(delx**2.+delz**2.)
         PO=CEXP(CMPLX(0.,K*RD))
         call QQ(RD,delz,K,AC,Q)
         PO=PO*Q/CMPLX(RD,0.)
         if(spew)write(2,*)'POc =',po
         PT=PT+PO
         if(spew)write(2,*)'PT, |PT| =',pt,cabs(pt)
      ENDIF
c
c Third possible path is from the source to a mirrored receiver.  Again,
c same formulae except for negative receiver height.
c
c     TH=ATAN2(D0+D1,HTS-HS)
c        {ANGLE FROM FIRST WEDGE-LEG TO SOURCE-RAY}
c     F0=TH-T1
c     TH=ATAN2(D2+D3,HTR+HR)
c        {ANGLE FROM FIRST WEDGE-LEG TO RECEIVER-RAY}
c     F1=2.*PI-T1-TH
      f0 = f0dir
      f1 = f1rfl
      AL1=PI+F0-F1
c
      if(spew)write(2,*)'d.I,F0,F1,AL1=',i,f0*radeg,f1*radeg,al1*radeg
      IF(AL1.GT.0.) THEN
c        {ADD REFLECTED G.O. FIELD OVER WEDGE}
         delx = recloc(1,2) - srcloc(1,1)
         delz = recloc(2,2) - srcloc(2,1)
         RD=SQRT(delx**2.+delz**2.)
         PO=CEXP(CMPLX(0.,K*RD))
         call QQ(RD,delz,K,C,Q)
c        if(spew)write(2,*)'HS,HR,HTR,HTS,sum=',
c    1              hs,hr,htr,hts,hs+hr+htr-hts
         if(spew)write(2,*)'RD,K,C,Q',rd,k,c,q
         PO=PO*Q/CMPLX(RD,0.)
         if(spew)write(2,*)'POd =',po
         PT=PT+PO
         if(spew)write(2,*)'PT, |PT| =',pt,cabs(pt)
      ENDIF
c
      delx = recloc(1,1) - srcloc(1,1)
      delz = recloc(2,1) - srcloc(2,1)
      RR=SQRT(delx**2.+delz**2.)
      level=4.34*LOG((RR*CABS(PT))**2.)
      if(spew)write(2,*)'RR, |PT| =', rr,cabs(pt)
c
      return
      end
