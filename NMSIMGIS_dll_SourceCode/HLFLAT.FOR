      subroutine hlflat(x1,z1,x2,z2,x3,z3,hs,bbs)
c
c This routine computes the "hill flat" component of the hill model,
c as defined by Equations 32 through 34 and Figure 4.  It is written
c for the source side of the system.  To obtain the value for the
c receiver side, call this with the various points complemented.c
c
c Modifications 3/16/96 by KJP:
c
c 1. Compute all angles based on the absolute values of the relevant
c    distances.  (Figure 4 indicates magnitudes of angles, and the
c    equations appear to be written for positive angles.)
c
c     Input:
c       x1          X position of source side of hill
c       z1          Z position of source side of hill
c       x2          X position of hill top
c       z2          Z position of hill top
c       x3          X position of receiver side of hill
c       z3          Z position of receiver side of hill
c       hs          Height of source, perpendicular to local surface.
c     Outputs:
c       bbs         "hill flat" component of the hill model
c
c Distance from source base to barrier top
      hhs = sqrt((x3-x1)**2 + (z3-z1)**2)
c Angle from source base to barrier top and angle of source plane
      beta13 = atan2(abs(z3-z1),abs(x3-x1))
      beta12 = atan2(abs(z2-z1),abs(x2-x1))
c Screen height re base-base line
      hscr = hhs*sin(abs(beta13 - beta12))
c Source base to screen base distance, along source baseline
      dscr = hhs*cos(abs(beta13-beta12))
c
      pfac = 57.29577951
      pfac2 = pfac**2
      bbs = pfac2*atan2(hs,dscr)*atan2(hscr,dscr)
      return
      end
