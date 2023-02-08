      subroutine setatm(absatm,freq,nfreq,tempc,rh)

      !DEC$ ATTRIBUTES DLLEXPORT::setatm

      implicit none

      integer,intent(in) :: nfreq
      real,intent(in) :: freq(nfreq), tempc, rh
      real,intent(inout) :: absatm(nfreq)
      
      real acomps(4),con,p,tempk
      real humid,alpha3,i,j

      p=101.325
      tempk = 273.15 + tempc
c
c Alpha3 returns the absorption coefficient in Nepers/m.  Conversion
c factor 2.649 is 8.69 dB/Neper divided by 3.28 feet/meter, so units
c are now dB/ft.
c
c Call alpha3 in "full" mode at the first frequency
      if(rh.gt.0) then
           con=humid(rh,p,tempk)
           absatm(1) = 8.69*alpha3(freq(1),p,tempk,con,acomps)
        endif
c Get the rest of them, calling alpha3 in "shortcut" mode
      do 80 j=2,nfreq
      if(rh.gt.0.) absatm(j) = 8.686*alpha3(freq(j),p,tempk,-1.,acomps)
 80   continue
c
      return
      end
