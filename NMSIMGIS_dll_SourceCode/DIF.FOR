      subroutine DIF(R1,A,AL2,PM,ANY,K,VD)
c         i   R1,A,AL2,PM,ANY,K
c         o   VD
c(* DF SVARER TIL DF I RAP. 111*)
c(*REGNER Pè KILEDIFFRAKTION*)
      real R1,A,AL2,PM,ANY,K
      real AA,XV,PI
      complex VD,DF,Y,FR
      integer NP,NM,N
c
      PI = 4.*ATAN(1.)
      DF=-CEXP(CMPLX(0.,K*R1+PI/4.0))/CMPLX(R1,0.)
      tangnt = TAN((PI+PM*AL2)/(2.0*ANY))
c     write(9,*)'pm,al2,any:',pm,al2,any
c     write(9,*)'tan:',tangnt
c Anti-blowup thing.  TAN = 0 is just grazing over the top of the barrier.
      if(tangnt.ne.0.) then
           AA=1.0/tangnt/(2.0*ANY)/SQRT(2.0*PI*K*A)
         else
           aa = aalast
        endif
      aalast = aa
      DF=DF*CMPLX(AA,0.)
c
      NP=0
      IF(AL2.GT.ANY*PI-PI) NP=1
      NM=-1
      IF(AL2.GT.PI-ANY*PI) NM=0
      IF(AL2.GT.PI+ANY*PI) NM=1
      N=NP
      IF(PM.LT.-0.9) N=NM
      XV=2.0*K*A*(COS((2.0*float(N)*ANY*PI-AL2)/2.0))**2.
      call FRES(SQRT(XV),FR)
c
      AA=-2.0*SQRT(XV)
      Y=CEXP(CMPLX(0.,-XV))
      Y=Y*FR
      VD=DF*Y
      VD=VD*CMPLX(0.,AA)
c
      return
      end
