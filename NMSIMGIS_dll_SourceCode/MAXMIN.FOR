      subroutine maxmin(prof,npts,hsrc,hrec,hgts,klocs)
c
c This routine examines profile PROF and returns the heights at several
c critical locations, indices LOCS.  For starters it finds the max and
c min heights of a source-receiver line above the profile.  It returns
c these as indices 1 and 2, respectively, of HGTS and LOCS.  HGTS repre-
c sents the height of the terrain above the line of sight.
c
c Arguments:
c
c  prof(2,npts)  i  x,z profile of terrain cut
c  npts          i  number of points in profile
c  hsrc          i  height of source above local ground
c  hrec          i  height of receiver above local ground
c  hgts(2)       o  highest (1) and lowest (2) profile points,
c                   relative to LOS from source to receiver
c  klocs(2)      o  indices (along prof) of HGTS
c
      real prof(2,npts),hgts(2),klocs(2)
c
c Dumb error kick-out: no array
c
      if(npts.le.1) then
           klocs(1) = 1
           klocs(2) = 1
           hgts(1) = - hsrc
           hgts(2) = - hrec
           return
        endif
c
c Test for the profile having a finite length; will be zero if directly
c overhead
c
      if(abs(prof(1,1) - prof(1,npts)) .lt. 10.) then
           hgts(1) = hsrc
           hgts(2) = hrec
           klocs(1) = 1
           klocs(2) = npts
           return
        endif
c
c Define parameters A and B such that height of the source-receiver line
c is given by z = Ax + B.
c
      x1 = prof(1,1)
      x2 = prof(1,npts)
      z1 = prof(2,1) + hsrc
      z2 = prof(2,npts) + hrec
      a = (z2 - z1)/(x2 - x1)
      b = z1 - a*x1
c
c Look for the max and min
c
      zmax = -10000.
      zmin = 10000.
      do 10, i = 1,npts
      z = a*prof(1,i) + b
      hgt = prof(2,i) - z
      if(hgt.gt.zmax) then
           zmax = hgt
           klocs(1) = i
        endif
c
      if(hgt.lt.zmin) then
           zmin = hgt
           klocs(2) = i
        endif
c
 10   continue
c
      hgts(1) = zmax
      hgts(2) = zmin
c
c      write(*,*)'hgts, locs =',hgts,klocs
c
      return
      end
