      subroutine SUBW(A,B,C,D,W)
c        i   A,B,C,D
c        o   W
c(* SUBW ER UNDERPROGRAM TIL W OG DERMED TIL QQ*)
      real A,B,C,D,AN,WR,WI
      complex W

      AN=(A**2.- B**2.-D)**2+(2.0*A*B)**2
      WR=C*B*(A**2.+ B**2.+D)/AN+REAL(W)
      WI=C*A*(A**2.+ B**2.-D)/AN+AIMAG(W)
      W=CMPLX(WR,WI)

      return
      end
