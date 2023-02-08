      real FUNCTION EFFV(arg1,arg2,arg3)

      real arg1,arg2,arg3      
      real*8 V,Z0,F
      real*8 LKG,VABS,HEFF,V0

      V=dble(arg1)
      Z0=dble(arg2)
      F=dble(arg3)
      IF(Z0.LT.0.0D0) THEN
c          {Lin var for negative z0}
         LKG=V
      ELSE
c          {effective height determination}
         HEFF=10.0D0*(340.0D0/F)
c          {determine v0 in log-formula}
         V0=V/DLOG(10.0D0/Z0)
c          {340/F= 1WVL}
         LKG=10.0D0*V0*(DLOG(HEFF/Z0)-DLOG(340.0D0/F/Z0))/
     1       (HEFF-340.0D0/F)
c          {RESULT IS EQUIV. WIND SPEED AT 10M}
         VABS=DABS(V)
         IF(LKG.GT.3.0D0*VABS) LKG=3.0D0*V
         IF(LKG.LT.-3.0D0*VABS) LKG=3.0D0*V
      ENDIF
      
      EFFV=sngl(LKG)
      
      END
