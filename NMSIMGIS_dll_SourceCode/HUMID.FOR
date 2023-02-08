      function humid(relhum,p,t)
C
****************************************************************************
*                                                                          *
*  Function HUMID(RELHUM,P,T) calculates molar concentration of water      *
*  vapor using the formulas found in ANSI s1.26-1978 (ASA 23-1978).        *
*                                                                          *
****************************************************************************
*                                                                          *
*                          VARIABLE DICTIONARY                             *
*                                                                          *
*  ARGUMENTS     i/o             DESCRIPTION OF VARIABLE                   *
*                                                                          *
*   RELHUM        i    relative humidity, either percent or fraction       *
*   P             i    atmospheric pressure, Pa                            *
*   T             i    atmospheric temperature, K                          *
*   HUMID         o    molar concentration of water vapor in percent or    *
*                      fraction, corresponding to the units of RELHUM      *
*                                                                          *
*  LOCAL VARIABLES               DESCRIPTION OF VARIABLE                   *
*                                                                          *
*   Po                 reference atmospheric pressure 101.325 kPa (1atm)   *
*   To                 reference ambient atmospheric temperature 293.15 K  *
*                                                                          *
****************************************************************************
*                                                                          *
C Triple point temperature (deg K) and reference pressure (Pa)
      data t01,ps0/273.16,101.325/
C
C Equation D11 for the common log of saturation pressure ratio
C
      
       arg=10.79586*(1.-t01/t) -5.02808*alog10(t/t01)
     1 + .000150474*(1.-10.**(-8.29692*((t/t01)-1.)))
     2 + .00042873*(-1+10.**(4.76955*(1.-(t/t01))))
     3 -2.2195983    
      
C
C Antilog for the ratio
      satrat=10.**arg
C
C The concentration
      humid=relhum*satrat*p/ps0
      return
      end
