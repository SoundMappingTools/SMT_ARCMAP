      subroutine SUBW2(A,B,C,D,W)
c        i   A,B,C,D
c        o   W
c(* SUBW2 ER UNDERPROGRAM TIL W OG DERMED TIL QQ*)
      real*8 A,B,C,D,AN,WR,WI
      complex*16 W

      AN=(A**2.0D0- B**2.0D0-D)**2.0D0+(2.0D0*A*B)**2.0D0
      WR=C*B*(A**2.0D0+ B**2.0D0+D)/AN+DREAL(W)
      WI=C*A*(A**2.0D0+ B**2.0D0-D)/AN+DIMAG(W)
      W=DCMPLX(WR,WI)
      
      return
      end
