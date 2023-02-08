      subroutine DELBAZ2(f,si,z)
c
c      inputs:
c                 f       Input frequency in Hz
c                 si      Ground Flow Resistivity in kN*s/m^4
c      output:    z       Complex reflection cooeficient
c
c
      real*8 f,si
      complex*16 z
      real*8 a,dumr,dumi
c
      a=si/f
      dumr=1.0D0 + 9.08D0 * a**0.75D0
c         {negative x because of exp(jwt) convention}
      dumi=-11.9D0 * a**0.73D0
c
      z=Dcmplx(dumr,dumi)
c
      return
      end