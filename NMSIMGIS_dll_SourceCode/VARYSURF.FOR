      subroutine VARYSURF(R,Gndtyp,N,Hsrc,Hrec,attsoft,atthard,Atten)
c
c  This subroutine calculates the acoustic attenuation across any ground
c  with varying surface type (ie: hard & soft ground combined).  The hard
c  and soft attenuations are supplied as inputs; this figgers the regions
c  of influence and does the weighted average.
c
c The theory upon which this subroutine is based is presented in the following:
c   Birger Plovsing,  July, 1994,
c  "Aircraft Sound Propagation over Non-Flat Terrain.  Prediction Algorithms,"
c   AV 7/94 - DELTA Acoustics & Vibration
c   Danish Acoustical Institute
c
c  Programmer:   Juliet Anne Page
c                Wyle Laboratories
c                Arlington, VA 22202
c                (703) 415-4550
c                13 October, 1994
c                Job# USAF-19092-2
c
c  VARIABLE DICTIONARY - Input / Output / iNternal - Integer / Real / Character
c  Variable  ION  IRC   Definition
c
c  Atten      O    R    Attenuation, dB from source to receiver, negative # returned
c  Atthard    N    R    Attenuation, dB Hard ground using LEVEL
c  Attsoft    N    R    Attenuation, dB Soft ground using LEVEL
c  Drefl      N    R    Distance from source to reflection point - method of mirrors
c  Gndtyp(500)I    N    Surface type at the R points, 0:  Soft Ground, ie: Grass
c                                                     1:  Hard Ground, ie: Water
c  Hrec       I    R    Receiver height defined normal to ground, m
c  Hsrc       I    R    Source height defined normal to ground, m
c  i          N    I    Loop counter through R points
c  N          I    I    Number of points along cut, point 1=source, point N=receiver
c  R(500)     I    R    Distance along cut from source to receiver, m, =0 at source
c                       may be non-uniformly spaced, defined along the flat ground
c
      dimension R(n)
      integer Gndtyp(n)
c
c      write(*,*)'N,Hsrc,Hrec,attsoft,atthard,Atten:'
c      write(*,*)N,Hsrc,Hrec,attsoft,atthard,Atten
c
c  Calculate the reflection point:
      drefl = Hsrc*R(N)/(Hsrc+Hrec)
c
c  Calculate boundary of Inactive / Active region on source side (Eqn 47-49)
c  Calculate boundary of Active / Inactive region on receiver side (Eqn 50-52)
      if (Gndtyp(1).eq.0) then
          if (hsrc.le.3.0) Ds = 10.0*Hsrc
          if (hsrc.gt.3.0) Ds = (20.*Hsrc/3. - 10.)*Hsrc
        else
          Ds = 30.*Hsrc
      endif
c
      if (gndtyp(n).eq.0) then
           if (hrec.le.3.0) Dr = 10.0*Hrec
           if (hrec.gt.3.0) Dr = (20.*Hrec/3. - 10.)*Hrec
         else
           dr = 30.*hrec
        endif
c
c  Check that Ds doesn't exceed 70% of Drefl & Dr doesn't exceed 70% of remainder
      ds = amin1(ds,0.7*drefl)
      dr = amin1(dr,0.7*(r(n)-drefl))
c
c Radii of the regions of influence, as measured from the source
c
      srcinf = ds
      recinf = r(n) - dr
c
c  Calculate the DL's for hard & soft surface based on LEVEL
c     write(*,*)'making LEVEL calls in varysurf'
c     call LEVEL(hsrc,hrec,R(N),0,Attsoft)
c     call LEVEL(hsrc,hrec,R(N),1,Atthard)
c     write(*,*)'attsoft,atthard =',attsoft,atthard
c
c  Copy the attenuations to delta L's
      dLsoft = Attsoft
      dLhard = Atthard
c
c  Evaluate via Equation (54)
c
      sum=0.
      denom=0.
      do 110, i=1,n
c  Bug out if not in either domain of influence
      if(r(i).lt.srcinf) go to 110
      if(r(i).gt.recinf) go to 120
c
c  Get weightings based on which side of reflection point
c
      if(r(i).lt.drefl) then
           w1 = hrec
         else
           w1 = hsrc
        endif
      if(gndtyp(i).eq.1) then
           w2 = .5
           deltl = dlhard
         else
           w2 = 1.0
           deltl = dlsoft
        endif
c
c  Add pieces to summation.  No need for a deltaR since the profile is
c  unifomly segmented
c
      sum = sum + w1*w2*deltl
      denom = denom + w1*w2
c      write(*,*)'sum,denom =',sum,denom
c
 110  continue
 120  continue
c
      deltal = sum/denom
      Atten = deltaL
c      write(*,*)'In varysurf: deltal,atten =',deltal,atten
      return
      end
