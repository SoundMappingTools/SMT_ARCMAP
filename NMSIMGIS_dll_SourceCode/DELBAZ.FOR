      subroutine DELBAZ(f,si,z)
c
c      inputs:
c                 f       Input frequency in Hz
c                 si      Ground Flow Resistivity in kN*s/m^4
c      output:    z       Complex reflection cooeficient
c
      real f,si
      complex z
      real a,dumr,dumi
c
      a=si/f
      dumr=1.0 + 9.08 * a**0.75
      dumi=11.9 * a**0.73
c
      z=cmplx(dumr,dumi)
c
      return
      end
