      subroutine QQ2(D,HS,HR,F,Z,Q)
c         i   D,HS,HR,F,Z
c         o   Q
      real*8 D,HS,HR,F
      real*8 R2,RR,RI,FR,FI,C,N,NYR,NYI,K
      real*8 PI,DUMR,DUMI, time0,time1,time2
      complex*16 Z,Q,T,W,ZERFE

      PI = 4.0D0*DATAN(1.0D0)
      K=2.0D0*PI*F/340.0D0
c      R1=DSQRT(D**2.0D0+(HS-HR)**2.0D0)
      R2=DSQRT(D**2.0D0+(HS+HR)**2.0D0)
      C=(HS+HR)/R2

      N=(DREAL(Z)*C+1.0D0)**2.0D0+(DIMAG(Z)*C)**2.0D0
c            {plane wave reflection coeff.}
      RR=((C*CDABS(Z))**2.0D0-1.0D0)/N
      RI=2.0D0*DIMAG(Z)*C/N
      
      NYR=DREAL(Z)/(CDABS(Z)**2.0D0)
      NYI=-DIMAG(Z)/(CDABS(Z)**2.0D0)
      
      DUMR=DSQRT(K*R2/4.0D0)*(NYR+C+NYI)
      DUMI=DSQRT(K*R2/4.0D0)*(-NYR-C+NYI)
      T=DCMPLX(DUMR,DUMI)
      
c          {EXP(-T2)*ERFC(JT)}
      call WW2(-T,W)

      FR=1.0D0+DSQRT(PI)*DIMAG(T*W)
      FI=     -DSQRT(PI)*DREAL(T*W)

c           {spherical wave refl.coef.accord. to c&s I}
      DUMR=RR+(1.0D0-RR)*FR+FI*RI
      DUMI=RI+FI*(1.0D0-RR)-RI*FR
      Q=DCMPLX(DUMR,DUMI)
      
      return
      END
