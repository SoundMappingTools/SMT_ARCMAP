      subroutine DAL(D0,D1,HS,HR,TR,SI0,SI1,F,level)
c(* BEREGNING AF LYDUBREDELSE HEN OVER FORHINDRING I TERRíNET EFTER *)
c(* PRINCIPPER BESKREVET I LI RAPPORT 111 FRA 1984. DEN GEOMETRISKE *)
c(* SITUATION ER BESKREVET I FIGUR B1 OG B2 SIDE 53 *)

c(* KNAEK SVARER TIL JONW I RAP.111 P.57, JVF FIG.B1 P. 53 *)
c(* I STEDET FOR AT HAVE Z SOM INDDATA HAR 'KILE' DOG SI (STRùMNINGSMODSTAND)*)
c(* SOM INDDATA*)
c(* VARIABELNAVNE FRA RAP. 111 ER STORT SET BIBEHOLDT*)
c(* REALDEL ER MARKERET MED ET EFTERHíNGT R OG IMAGINíRDEL*)
c(* MED ET EFTERHíNGT I *)
c
c     Inputs:
c       D0          Distance from the first point in the profile to the second
c       D1          Distance from the second point in the profile to the third
c       HS          Height of the source in meters
c       HR          Height of the receiver in meters
c       TR          Angle between plane of source and plane of slope
c       SI0         Ground Flow Resistivity in kN*s/m^4
c       SI1         Ground Flow Resistivity in kN*s/m^4
c       F           Frequency of interest (Hz)
c     Outputs:
c       Level       Ground attenuation in dB

      real D0,D1,HS,HR,SI0,SI1,F,level
      real R1,RS,RB
      real A,HX,K,RR
      real RD,F0,F1,TR
      real RH0,RH1,ANY,THETA
      real PI
      complex AC,B,PL,P0,Q1,Q2,PS,PR,PB,Q,VD
      
      PI = 4.*ATAN(1.)
      
c      (*DELANY-BAZLEY FOR SOURCE LEG*)
      call DELBAZ(F,SI0,AC)

c      (*DELANY-BAZLEY FOR RECEIVER LEG*)
      call DELBAZ(F,SI1,B)
      
c          (* K IS WAVENO.*)
      K=2.0*PI*F/340.0
      PL=CMPLX(0.,0.)
      RH0=SQRT(HS**2.+D0**2.)
      RH1=SQRT(HR**2.+D1**2.)
      R1=RH0+RH1
      F0=ATAN2(HS,D0)
      F1=ATAN2(HR,D1)
      F1=2.*PI-TR-F1

      RD=SQRT(RH0**2.+RH1**2.-2.*RH0*RH1*COS(F1-F0))
c      {DIRECT FIELD}
      P0=CEXP(CMPLX(0.,K*RD))/CMPLX(RD,0.)

c        {HEIGHT OVER WEDGE-LEG}
      HX=SIN(F0)*R1

c        {FINITE WEDGE IMPEDANCE}
      call QQ(R1,HX,K,AC,Q1)
c        {HEIGHT OVER WEDGE-LEG}
      HX=SIN(2.*PI-F1-TR)*R1

c        {FINITE WEDGE IMPEDANCE}
      call QQ(R1,HX,K,B,Q2)
      A=RH0*RH1/R1
      ANY=2.0-TR/PI

      call DIF(R1,A,F1-F0,-1.0,ANY,K,VD)
      PL=VD
      call DIF(R1,A,F1+F0,-1.0,ANY,K,VD)
      VD=VD*Q1

c                {DIFFRACTION}
      PL=PL+VD
      call DIF(R1,A,F1+F0,1.0,ANY,K,VD)
      VD=VD*Q2
      PL=PL+VD
      
      call DIF(R1,A,F1-F0,1.0,ANY,K,VD)
      VD=VD*Q1*Q2
      PL=PL+VD

C           {ADDING DIRECT FIELD IF TEST IS TRUE}
      IF((F1-F0).LT.PI) PL=PL+P0

      IF((F1+F0).LT.PI) THEN
         RS=SQRT(RH0**2.+RH1**2.-2.0*RH0*RH1*COS(F1+F0))
c           {REFLECTION IN WEDGE LEG BELOW SOURCE}
         PS=CEXP(CMPLX(0.,K*RS))/CMPLX(RS,0.)
         call QQ(RS,HS+RH1*SIN(PI-F1),K,AC,Q)
         PS=PS*Q
         PL=PL+PS
      ENDIF

      THETA=ATAN2(HR,D1)

      IF((F1-F0+2.*THETA).LT.PI) THEN
         RR=SQRT(RH0**2.+RH1**2.-2.*RH0*RH1*COS(F1-F0+2.*THETA))
c            {REFLECTION IN WEDGE LEG BELOW RECEIV}
         PR=CEXP(CMPLX(0.,K*RR))/CMPLX(RR,0.)
         call QQ(RR,HR+RH0*SIN(F0+TR-PI),K,B,Q)
         PR=PR*Q
         PL=PL+PR
      ENDIF

      IF((F1+F0+2.*THETA).LT.PI) THEN
         RB=SQRT(RH0**2.+RH1**2.-2.*RH0*RH1*COS(F1+F0+2.*THETA))
c           {REFLECTION IN BOTH WEDGE LEGS}
         PB=CEXP(CMPLX(0.,K*RB))/CMPLX(RB,0.)
         call QQ(RB,HS+RH1*SIN(2.*TR-3.*PI+F1),K,AC,Q)
         PB=PB*Q
         call QQ(RB,HR+RH0*SIN(-F0+TR-PI),K,B,Q)
         PB=PB*Q
         PL=PL+PB
      ENDIF

      level=4.34*LOG((RD*CABS(PL))**2.)
      
      return
      END
      
      