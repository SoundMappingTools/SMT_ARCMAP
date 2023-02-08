      subroutine WW2(T,W)
c         i   T
c         o   W

c(* W ER UNDERPROGRAM TIL QQ*)
      real*8 A1,B1,C1,D1,P2,Q2,BB,AA,AN,AH,AK,A,B,C,D
      real*8 WR,WI,PI,dum
      integer N
      complex*16 T,W
      
      PI = 4.0D0*DATAN(1.0D0)
      A=DABS(DREAL(T))
      B=DABS(DIMAG(T))
c      A=(DREAL(T))
c      B=(DIMAG(T))
      
      
      IF((A.LE.3.9D0).AND.(B.LE.3.0D0)) THEN
         A1=DCOS(2.0D0*A*B)
         B1=DSIN(2.0D0*A*B)
         C1=DEXP(-2.0D0*B*PI/0.8D0)
         C1=C1-DCOS(2.0D0*A*PI/0.8D0)
         D1=DSIN(2.0D0*A*PI/0.8D0)
         P2=2.0D0*DEXP(-(A**2.0D0+2.0D0*B*PI/0.8D0-B**2.0D0))
         P2=P2*((A1*C1-B1*D1)/(C1**2.0D0+D1**2.0D0))
         Q2=2.0D0*DEXP(-(A**2.0D0+2.0D0*B*PI/0.8D0-B**2.0D0))
         Q2=Q2*((A1*D1+B1*C1)/(C1**2.0D0+D1**2.0D0))
         BB=B**2.0D0
         AA=A**2.0D0
         AH=0.8D0*B/(PI*(AA+BB))
         AK=0.8D0*A/(PI*(AA+BB))
         
         do 100 N=1,5
            dum=DBLE(N)
            AN=(BB-AA+dum**2.0D0* 0.64D0)**2.0D0 + 4.0D0*AA*BB
            AH=AH+
     1         2.0D0*B*0.8D0/PI*DEXP(-dum**2.0D0*0.64D0)*
     2         (BB+AA+0.64D0*dum**2.0D0)/AN
            AK=AK+
     1         2.0D0*A*0.8D0/PI*DEXP(-dum**2.0D0*0.64D0)*
     2         (BB+AA-0.64D0*dum**2.0D0)/AN
 100     continue
 
         IF(B.LT.PI/0.8D0) THEN
            AH=AH+P2
            AK=AK-Q2
         ENDIF
         
         IF(B.EQ.PI/0.8D0) THEN
            AH=AH+P2/2.0D0
            AK=AK-Q2/2.0D0
         ENDIF
         
         W=DCMPLX(AH,AK)
         
      ELSE
      
         W=DCMPLX(0.0D0,0.0D0)
         C=0.4613135D0
         D=0.1901635D0
         call SUBW2(A,B,C,D,W)
         C=0.09999216D0
         D=1.7844927D0
         call SUBW2(A,B,C,D,W)
         C=0.002883894D0
         D=5.5253437D0
         call SUBW2(A,B,C,D,W)
         
      ENDIF
      
      
      IF(DREAL(T).LT.0.0D0) THEN
         W=DCONJG(W)
         A=-A
      ENDIF
      
      IF(DIMAG(T).LT.0.0D0) THEN
         WR=2.0D0*DEXP(B**2.0D0-A**2.0D0)*DCOS(2.0D0*A*B)-DREAL(W)
         WI=2.0D0*DEXP(B**2.0D0-A**2.0D0)*DSIN(2.0D0*A*B)+DIMAG(W)
         W=DCMPLX(WR,WI)
      ENDIF
      
      return
      end
