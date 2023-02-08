      subroutine incidence(xyzsrcs,xyzrecs,zsrcs,zrecs,headings,
     &                           pitchs,roll,distance,thrad,phi,elvang)
      

      implicit none
!
!     This subroutine will calculate the incidence angle from the source
!     to the receiver, computing the Phi and Theta position on the source
!
! KABnote: this assumes the following:
!          xyzsrc() is the ground location under the aircraft (source)
!          xyzrec() is the ground location under the receiver

!          zsrc is the height above ground of the source     
!          zrec is the height above ground of the receiver     
!
      real xyzsrcs(3),xyzrecs(3),zsrcs,zrecs,headings,pitchs,roll,elvang

      real*8 xyzsrc(3),xyzrec(3),zsrc,zrec,heading,pitch

      real*8 rx, ry, rz, UnitX, UnitY, UnitZ
      real*8 coshead, denom, dir, fnum, dphi, dhead   

      real*8 arg
      real*8 degrad, degpi 

      real distance,thrad, phi, elevang
      
      data degrad/.0174532925/
      data degpi/57.2957795131/
 
      pitch=dble(pitchs) * degrad
      heading=dble(90.0d0 - headings) * degrad

      xyzsrc=dble(xyzsrcs)
      xyzrec=dble(xyzrecs)
      zsrc=dble(zsrcs)
      zrec=dble(zrecs)


      distance=sngl(dsqrt((xyzsrc(1)-xyzrec(1))**2+
     &              (xyzsrc(2)-xyzrec(2))**2+
     &              ((xyzsrc(3)+zsrc)-(xyzrec(3)+zrec))**2))

c Create unit vector components from source to receiver
      rx = (xyzrec(1) - xyzsrc(1))/dble(distance)
      ry = (xyzrec(2) - xyzsrc(2))/dble(distance)
      rz = ((xyzrec(3)+zrec)-(xyzsrc(3)+zsrc))/dble(distance)
    
      if(distance.gt.abs((xyzsrc(3)+zsrc)-(xyzrec(3)+zrec))) then
        arg = ((xyzsrc(3)+zsrc)-(xyzrec(3)+zrec))/dble(distance)
      else
        arg = 1.d0
      endif
      arg = max(arg,-1.d0)
      arg = min(arg,1.d0)
      
      elvang = asin(sngl(arg))

c Below are the three components of the AC direction vector

ckab subroutine angles uses tau for the first argument used below  
ckab where tau is defined as elevation angle + angle of attack.
ckab Below we are using pitch since that is provided

ckab subroutine angles uses psi as the second argument below 
ckab where psi is defined as the heading angle + yaw angle.
ckab Below, we are using heading since that is provided
ckab These can be refined to include angle of attack & sideslip, as required. 
ckab Or, it's probably better to deal with AOA & sideslip in the wind frame. 

      UnitX = dcos(pitch)*cos(heading)
      UnitY = dcos(pitch)*sin(heading)
      UnitZ = dsin(pitch)
c
c Use rx,ry,rz from above as the three components of the vector from the AC to receiver       

c
c  The dot product is the cosine of theta.
c  Equivalent to the earth"s longitude, when the aircraft points due North.
      coshead = unitX*rx + unitY*ry + unitZ*rz


      dhead = dacos(dmax1(dmin1(coshead,1.0d0),-1.0d0))

      dir  = UnitY*rx - UnitX*ry
c
c  Calculate Phi.
      denom = dble(distance)*dsin(dhead)

c check this!
      if(dhead .le. 0.0d0) then
         thrad = 0.0
         heading = 0.0
c         heading = 0.0 
         phi = 0.0 + roll
         return
ckab in the real version of this, just enable the return and cut out the following goto

      end if


      
c
c the original equation for fnum uses zsrc-zrec to compute the elevation angle.
c this version assumes the height difference is as shown below. 


      fnum = dble((xyzsrc(3)+zsrc)-(xyzrec(3)+zrec))/dcos(pitch)
     +     + dble(distance)*coshead*dtan(pitch)

c
c  Arc sin argument gives phi which is equivalent to the latitude.
c  90 degrees is pointing directly under the aircraft.
      dphi = dacos(dmin1(1.0d0,dmax1(-1.0d0,fnum/denom)))

c
c  Make phi=0 when observer is directly under flight and
c  positive when observer is on the right side of the path.
c  Z part of cross product gives vector direction
      dphi = dsign(dphi,dir)
c
c  Convert theta and phi to degrees
      thrad = sngl(dhead)
      phi = sngl(degpi*dphi + roll)

      end subroutine