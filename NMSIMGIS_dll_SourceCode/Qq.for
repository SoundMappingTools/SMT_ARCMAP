      subroutine QQ(R2,Hin,K,Z,Q)
c         i    R2,H,K,Z
c         o    Q
c(* QQ UDREGNER SFíRISK BùLGE REFLEKSIONSKOEFF. Q *)
      real R2,H,K,C,N,DUMR,DUMI,GAR,GAI,NYR,NYI,FR,FI
      real PI
      complex Z,Q,T,W,CERFE,value
	real*8 time0,time1,time2

      h=abs(hin)

      PI = 4.*ATAN(1.)
      C=H/R2
      N=(REAL(Z)*C+1.0)**2.+(AIMAG(Z)*C)**2.

      GAR=((C*CABS(Z))**2.-1.0)/N
      GAI=2.0*AIMAG(Z)*C/N

      NYR=  REAL(Z)/(CABS(Z)**2.)
      NYI=-AIMAG(Z)/(CABS(Z)**2.)

      DUMR=SQRT(K*R2/4.0)*(NYR+C-NYI)
      DUMI=SQRT(K*R2/4.0)*(NYR+C+NYI)
      T=CMPLX(DUMR,DUMI)

      call WW(T,W)

      FR=1.0+SQRT(PI)*(-AIMAG(T*W))
      FI=    SQRT(PI)*REAL(T*W)

      DUMR=GAR+(1.-GAR)*FR+FI*GAI
      DUMI=GAI+FI*(1.-GAR)-GAI*FR
      Q=CMPLX(DUMR,DUMI)

      IF(REAL(Z).GE.1000.0) Q=CMPLX(1.,0.)

      return
      end
