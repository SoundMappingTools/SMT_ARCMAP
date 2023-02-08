      subroutine EGAL(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,
     1                arg10,arg11,arg12)
c        i  D1,D2,HS,HR,SI1,SI2,VA,APH,MN,F
c        o  LEVTURB,LNOT
c{CALCULATION FOR POINT SOURCE OVER IMPEDANCE PLANE ACCORDING
cTO CHIEN AND SOROKA I AND DELANY-BAZLEY AND KBR. Empirical Wind
c model taken into acc.}
c{12/12/1991 version with limited gradient height}
c{kbr-turbulence model included}
c{freq-dependant gradient included}
c
c     This is the main routine for calculating the ground effect for 
c     flat ground.  It takes multi-path interacting into account.

c     Inputs:
c       arg1              Half of the distance from the source to rhe receiver
c       arg2              Half of the distance from the source to rhe receiver
c       arg3              Height of source in meters
c       arg4              Height of receiver in meters
c       arg5              Ground flow resistivity in kN*s/m^4
c       arg6              Ground flow resistivity in kN*s/m^4
c       arg7              Effective wind velocity (usually 0.0)
c       arg8              Air turbulence parameter (usually 0.0)
c       arg9              Air turbulence parameter (usually 0.0)
c       arg10             Frequncy of interest (Hz)
c       arg11             dummy for level (dB) relative to the free field with
c                           turbulence [EGAL subroutine]

c
c     Outputs:
c       arg12             Ground attenuation at the given frequency in dB
c
      real arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11
      real arg12
      real*8 D1,D2,VA,APH,HS,HR,SI1,SI2,F,MN
      real*8 LEVTURB,LNOT
      real*8 RR,RD,K,AR,PI,LA,HA,RR1,RR2,RR3,RR4,CR,CRH
      real*8 HL,L,KV(3),V(3),NFF
      integer M,J
      complex*16 RX1,RX2,DB(3),JARAY(3),JH,CARAY(3),CH,Q1,Q2,C

      PI = 4.D0*DATAN(1.D0)
      D1=dble(arg1)
      D2=dble(arg2)
      HS=dble(arg3)
      HR=dble(arg4)
      SI1=dble(arg5)
      SI2=dble(arg6)
c  KBR initializes V to -4
c  Non-zero V requires a value for APH
      VA=dble(arg7)
      APH=dble(arg8)
      MN=dble(arg9)
      F=dble(arg10)


c          Delany-Bazley for Source side
      call DELBAZ2(F,SI1,RX1)

c          Delany-Bazley for Receiver side
      call DELBAZ2(F,SI2,RX2)

      L=D1+D2
      RD=DSQRT(L**2.0D0+(HS-HR)**2.0D0)
      RR=DSQRT(L**2.0D0+(HS+HR)**2.0D0)

      LA=340.0D0/F
      K=2.0D0*PI*F/340.0D0
      IF(VA.NE.0.0D0) THEN
c           {APERTURE HEIGHT: APH} {NOTE: 6 POINTS PER WVLENGTH}
c           M is made an integer.
ckjp 10/28/96 changed dnint to idnint
         M=iDnint(APH/(LA/6.0D0))
c           {FIXED SPEED}
         K=2.0D0*PI*F/(340.0D0+((DBLE(M-1)*(LA/6.0D0)+LA/10.0D0)
     1    /2.0D0)*VA/10.0D0)
c           {HL IS effective TURBULENCE SCALE}
c           {TURB SCALE INV. PROP WITH FREQ}
         HL=1000.0D0/F
c           Initialize Variables
         HL=0.0D0
         AR=0.0D0
         DO 150 I=1,3
 150     DB(I)=DCMPLX(0.D0,0.D0)
      ENDIF
      IF(VA.EQ.0.0D0) GOTO 300

c ************************ PRIMARY LOOP *****************************
      DO 200 J=1,M

         HA=DBLE(J-1)*(LA/6.0D0)+LA/10.0D0

c           {wind speed is varying}
         V(1)=VA+MN/10.0D0*DCOS(HA*2.0D0*PI/HL)
c           {-COS INSTEAD OF COS}
         V(2)=2.0D0*VA-V(1)
c           {THIRD PASS: NO TURB}
         V(3)=VA

         DO 175 I=1,3

            KV(I)=2.0D0*PI*F/(340.0D0+(HA/2.0D0)*V(I)/10.0D0)

            RR1=DSQRT(D1**2.0D0+(HS-HA)**2.0D0)
            RR2=DSQRT(D1**2.0D0+(HS+HA)**2.0D0)
            RR3=DSQRT(D2**2.0D0+(HA-HR)**2.0D0)
            RR4=DSQRT(D2**2.0D0+(HA+HR)**2.0D0)

            JARAY(I)=CDEXP(DCMPLX(0.0D0,-KV(I)*(RR1+RR3)))
            JH=CDEXP(DCMPLX(0.D0,-K*(RR1+RR3)))
            JARAY(I)=JARAY(I)-JH

            CR=RR3*DSQRT(RR3*RR1*(RR1+RR3))
            JARAY(I)=JARAY(I)/DCMPLX(CR,0.0D0)
            DB(I)=DB(I)+JARAY(I)

c              finite wedge impedance
c              {Q IS AS USUAL CALC FOR FIXED K}
            call QQ2(D2,HA,HR,F,RX2,Q2)
            CARAY(I)=CDEXP(DCMPLX(0.0D0,-KV(I)*(RR1+RR4)))
            CH=CDEXP(DCMPLX(0.0D0,-K*(RR1+RR4)))
            CARAY(I)=CARAY(I)-CH
            JARAY(I)=CARAY(I)*Q2

            CR=RR4*DSQRT(RR4*RR1*(RR1+RR4))
            JARAY(I)=JARAY(I)/DCMPLX(CR,0.0D0)
            DB(I)=DB(I)+JARAY(I)

c              finite wedge impedance
c              { Q IS AS USUAL CALC FOR FIXED K}
            call QQ2(D1,HS,HA,F,RX1,Q1)
            CARAY(I)=CDEXP(DCMPLX(0.0D0,-KV(I)*(RR2+RR3)))
            CH=CDEXP(DCMPLX(0.0D0,-K*(RR2+RR3)))
            CARAY(I)=CARAY(I)-CH
            JARAY(I)=CARAY(I)*Q1

            CR=RR3*SQRT(RR3*RR2*(RR2+RR3))
            JARAY(I)=JARAY(I)/DCMPLX(CR,0.0D0)
            DB(I)=DB(I)+JARAY(I)

            CARAY(I)=CDEXP(DCMPLX(0.0D0,-KV(I)*(RR2+RR4)))
            CH=CDEXP(DCMPLX(0.0D0,-K*(RR2+RR4)))
            CARAY(I)=CARAY(I)-CH
            JARAY(I)=CARAY(I)*Q1*Q2

            CR=RR4*DSQRT(RR4*RR2*(RR2+RR4))
            JARAY(I)=JARAY(I)/DCMPLX(CR,0.0D0)
            DB(I)=DB(I)+JARAY(I)

 175     continue

 200  continue
c ********************** END OF PRIMARY LOOP ************************

c        {+jwt convention}
      C=CDEXP(DCMPLX(0.0D0,PI/4.0D0))
c        { from JSV/KBR 1982 let.}
      CRH=(LA/6.0D0)*D2*DSQRT(8.0D0*PI*K)/(16.0D0*PI**2.0D0)
      C=C*DCMPLX(CRH,0.0D0)
c        {END OF APERTURE CALC}
      DB(1)=DB(1)*C
      DB(2)=DB(2)*C
      DB(3)=DB(3)*C

c        finite wedge impedance
 300  call QQ2(L,HS,HR,F,RX2,Q2)
      C= CDEXP(DCMPLX(0.0D0,-K*RR))/DCMPLX(RR,0.0D0)*Q2
      CH=CDEXP(DCMPLX(0.0D0,-K*RD))/DCMPLX(RD,0.0D0)
      C=C+CH

c        {homogeneous still atmosphere contrib.}
      C=C/DCMPLX(4.0D0*PI,0.0D0)

      NFF=16.0D0*PI**2.0D0*RD**2.0D0
      IF(VA.NE.0.0D0) THEN
c        {APERTURE CONTRIB ADDED}
         DB(1)=DB(1)+C
         DB(2)=DB(2)+C
         DB(3)=DB(3)+C
         AR=CDABS(DB(1))**2.0D0+CDABS(DB(2))**2.0D0+CDABS(DB(3))**2.0D0
c           {result with turbulence taken into acc.}
c           {16*PI*PI*RD*RD Represents norm with free field}
         LEVTURB=4.3429D0*DLOG((AR*NFF)/3.0D0)
c           {Level without turb}
         LNOT=4.3429D0*DLOG(CDABS(DB(3))**2.0D0*NFF)
      ELSE
c           {STILL ATMOS}
         LEVTURB=0.0D0
         LNOT=4.3429D0*DLOG(CDABS(C)**2.0D0*NFF)
      ENDIF

      arg11=sngl(LEVTURB)
      arg12=sngl(LNOT)

      return
      END
