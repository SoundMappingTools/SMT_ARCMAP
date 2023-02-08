      subroutine WW(T,W)
c         i   T
c         o   W
c(* W ER UNDERPROGRAM TIL QQ*)
      real A1,B1,C1,D1,P2,Q2,BB,AA,AN,AH,AK,A,B,C,D
      real WR,WI,PI,dum
      integer N
      complex T,W
      
      PI = 4.*ATAN(1.)
      A=ABS(REAL(T))
      B=ABS(AIMAG(T))
      
      
      IF((A.LE.3.9).AND.(B.LE.3.0)) THEN
         A1=COS(2.0*A*B)
         B1=SIN(2.0*A*B)
         C1=EXP(-2.0*B*PI/0.8)-COS(2.0*A*PI/0.8)
         D1=SIN(2.0*A*PI/0.8)
         P2=2.0*EXP(-(A**2.+2.0*B*PI/0.8-B**2.))
         P2=P2*((A1*C1-B1*D1)/(C1**2.+D1**2.))
         Q2=2.0*EXP(-(A**2.+2.0*B*PI/0.8-B**2.))
         Q2=Q2*((A1*D1+B1*C1)/(C1**2.+D1**2.))
         BB=B**2.
         AA=A**2.
         AH=0.8*B/(PI*(AA+BB))
         AK=0.8*A/(PI*(AA+BB))
         
         do 100 N=1,5
            dum=float(N)
            AN=(BB-AA+dum**2.* 0.64)**2. + 4.0*AA*BB
            AH=AH+
     1         2.0*B*0.8/PI*EXP(-dum**2.*0.64)*(BB+AA+0.64*dum**2.)/AN
            AK=AK+
     1         2.0*A*0.8/PI*EXP(-dum**2.*0.64)*(BB+AA-0.64*dum**2.)/AN
 100     continue
 
         IF(B.LT.PI/0.8) THEN
            AH=AH+P2
            AK=AK-Q2
         ENDIF
         
         IF(B.EQ.PI/0.8) THEN
            AH=AH+P2/2.0
            AK=AK-Q2/2.0
         ENDIF
         
         W=CMPLX(AH,AK)
         
      ELSE
      
         W=CMPLX(0.,0.)
         C=0.4613135
         D=0.1901635
         call SUBW(A,B,C,D,W)
         C=0.09999216
         D=1.7844927
         call SUBW(A,B,C,D,W)
         C=0.002883894
         D=5.5253437
         call SUBW(A,B,C,D,W)
         
      ENDIF
      
      
      IF(REAL(T).LT.0.) THEN
         W=CONJG(W)
         A=-A
      ENDIF
      
      IF(AIMAG(T).LT.0.) THEN
         WR=2.0*EXP(B**2.-A**2.)*COS(2.0*A*B)-REAL(W)
         WI=2.0*EXP(B**2.-A**2.)*SIN(2.0*A*B)+AIMAG(W)
         W=CMPLX(WR,WI)
      ENDIF
      
      return
      end
