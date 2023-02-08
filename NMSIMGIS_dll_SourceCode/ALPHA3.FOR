      FUNCTION ALPHA3(F,P,T,CON,A)
C
****************************************************************************
*                                                                          *
*  Function ALPHA3(F,P,T,CON,A) calculates absorption of sound by the      *
*  atmosphere using the formulas found in ANSI s1.26-1978 (ASA 23-1978)    *
*                                                                          *
****************************************************************************
*                                                                          *
*                          VARIABLE DICTIONARY                             *
*                                                                          *
*  ARGUMENTS     i/o             DESCRIPTION OF VARIABLE                   *
*                                                                          *
*   F            i     acoustic frequency, Hz                              *
*   P            i     atmospheric pressure, Pa                            *
*   T            i     atmospheric temperature, K                          *
*   CON          i     molar concentration of water vapor in percent, or,  *
*                      for a given temperature, 100 times the ratio of     *
*                      the number of moles of water present in an          *
*                      atmospheric sample to the total number of moles in  *
*                      that sample, %                                      *
*   ALPHA        o     absorption coefficient, Np m^-1 (nepers per meter)  *
*   A(4)         o     classical, oxygen relaxation, and nitrogen relax-   *
*                      ation, respectively, components of ALPHA, and ALPHA *
*                      itself, all normalized by frequency squared.        *
*                                                                          *
*                                                                          *
*  LOCAL VARIABLES               DESCRIPTION OF VARIABLE                   *
*                                                                          *
*   Po                 reference atmospheric pressure 101.325 kPa (1atm)   *
*   To                 reference ambient atmospheric temperature 293.15 K  *
*   Fro                vibration relaxation frequency of O2                *
*   Frn                vibration relaxation frequency of N2                *
*                                                                          *
*                                                                          *
****************************************************************************
C
      save  !This way we only compute them once.
      
      dimension a(4)
      data po,t0/101.325,293.15/
C  Time saving switch: use previous FRs if negative humidity
      if(con.lt.0.) go to 10
C  Calculate the relaxation frequencies given in ANSI S1.26-1978,
C  section 4.0, equations 8 and 9.
c      fro = (p/po)*(24. + 4.41e04 * con * ((0.05 + con)/(0.391 + con)))
c      frn = (p/po)*sqrt(t0/t)*(9. + 350. * con *
c     1        exp(-6.142*(((t/t0)**(-1./3.)) - 1.)))

C  Calculate the relaxation frequencies given in ANSI S1.26-1995 (R2004),
C  section 4.0, equations 8 and 9.
      fro = (p/po)*(24. + (4.04e04 * con *((0.02 + con)/(0.391 + con))))
      frn = (p/po)*sqrt(t0/t)*(9. + 280. * con *
     1        exp(-4.170*(((t/t0)**(-1./3.)) - 1.)))

C  Calculate the absoption coefficient according to ANSI S1.26-1995 (R2004),
C  section 4.0, equation 10.
 10   ff = f*f
C
C Equation 10 as published:
C     alpha = ff * (1.84e-11 * (po/p) * sqrt(t/to)
C    1        + ((t/to)**(-2.5))*(.01278 * (exp(-2239.1/t))/
C    2        (fro + (ff/fro)) + .1068 * (exp(-3352./t))/
C    3        (frn + (ff/frn))))
C
C The three components, normalized by ff:
      a(1) = 1.84e-11 * (po/p) * sqrt(t/t0)
      a(2)=((t/t0)**(-2.5))*(.01275*(exp(-2239.1/t))/(fro + (ff/fro)))
      a(3)=((t/t0)**(-2.5))*(.1068*(exp(-3352./t))/(frn + (ff/frn)))
C The sum
      a(4) = a(1)+a(2)+a(3)

C And alpha itself
      alpha3 = a(4)*ff
      return
      end