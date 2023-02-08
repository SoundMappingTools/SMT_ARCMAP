      subroutine FRES(Y,FR)
c         i   Y
c         o   FR
c(* BEREGNING AF FRESNELINTEGRAL --BRUGES I VD PROCEDURE*)
c(* FRR ER REALDELEN OG FRI ER IMAGIN’RDELEN AF FRESNELINTEGR.*)
      real Y,C,X,F,G,SI,CO
      complex FR
      
      C=0.797885
      X=C*Y
      F=(1.0+0.926*X)/(2.0+1.792*X+3.104*X**2.)
      G=1.0/(2.0+4.142*X+3.492*X**2.+6.67*X**3.)
      SI=SIN((X/C)**2.)
      CO=COS((X/C)**2.)
      FR=CMPLX((-F*SI+G*CO)/C,(F*CO+G*SI)/C)
      
      return
      end
