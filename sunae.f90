! Copyright (C) 2011, Juan Pablo Justiniano  <jpjustiniano@gmail.com>
!
Program TestSunae
implicit none

 Real :: year
 Real :: day
 Real :: hour
 Real :: lat
 Real :: lon
 Real :: az
 Real :: el
 Real :: ha
 Real :: dec
 Real :: soldst

	Year = 2011
	lat = -33.45
	lon= -70.66
	
	write(*,*)
40	write(*,*) 'day of year:  '
	read (*,*) day
50	write(*,*) 'time:       (UTC)  (99 = exit) (999 = change day)'
	read (*,*) hour
	If (hour == 99) goto 99
	If (hour == 999) goto 40
	
	call sunae(year,day,hour,lat,lon,az,el,ha,dec,soldst)
	write (*,*)
	If (el<0.) then
		print *,' Night time!! Remember, UTC decimal hours.!!'
		go to 50
	end if	
	write (*,100) az,el,ha,dec,soldst
	write (*,*)
	goto 50
	
100 Format (' Azimuth:    ',F6.1,/' Elevation:  ', F4.1,/' Hour angle: ',F5.1,/ &
			& ' Declination: ', F4.1,/' Sun dist.:   ',F4.1)

99 stop


 contains
!-------------------------------------------------------------------------------+
!  This subroutine calculates the local azimuth and elevation of the sun at		|
!  a specific location and time using an approximation to equations used   		|
!  to generate tables in The Astronomical Almanac.  				  	   		|
!  Refraction correction is added so sun position is apparent one. 		   		|
!                                                                      			|
!The Astronomical Almanac, U.S. Gov't Printing Office, Washington,D.C. (1985).	|
!                                                                      			|
!    input parameters                                                  			|
!      year=year, e.g., 1986                                           			|
!      day=day of year, e.g., feb 1=32                                 			|
!      hour=hours plus fraction in UT, e.g., 8:30 am eastern daylight  			|
!        time is equal to 8.5 + 5(5 hours west of Greenwich) -1(for    			|
!        daylight savings time correction)                             			|
!      lat=latitude in degrees (north is positive)                     			|
!      long=longitude in degrees (east is positive)                    			|
!                                                                      			|
!    output parameters                                                 			|
!	az= sun azimuth angle (measured east from north, 0 to 360 degs)				|
!	el= sun elevation angle (degs)                                 				|
!	ha= solar hour angle														|
!	dec = declination															|
!	soldst = solar distance                                        				|				
!-------------------------------------------------------------------------------+
 Subroutine sunae(year,day,hour, lat, long,az,el,ha,dec,soldst)  ! Sun's position, Michalsky
 implicit none
 Real, intent(in) :: year
 Real, intent(in) :: day
 Real, intent(in) :: hour
 Real, intent(in) :: lat
 Real, intent(in) :: long
 Real, intent(out) :: az
 Real, intent(out) :: el
 Real, intent(out) :: ha
 Real, intent(out) :: dec
 Real, intent(out) :: soldst
 
 real :: twopi, pi, rad
 real :: delta, leap, jd, time
 real :: mnlong, mnanom, eclong, oblqec, num, den, ra
 real :: gmst, lmst, latrad, elc, refrac, soldia

 data twopi,pi,rad/6.2831853,3.1415927,.017453293/
  
!   get the current julian date (actually add 2,400,000 for jd)
      delta=year-1949.
      leap=aint(delta/4.)
      jd=32916.5+delta*365.+leap+day+hour/24.
!   1st no. is mid. 0 jan 1949 minus 2.4e6; leap=leap days since 1949
!  the last yr of century is not leap yr unless divisible by 400
      if (amod(year,100.).eq.0.0.and.amod(year,400.).ne.0.0) jd=jd-1.

!   calculate ecliptic coordinates
      time=jd-51545.0
!   51545.0 + 2.4e6 = noon 1 jan 2000

!   force mean longitude between 0 and 360 degs
      mnlong=280.460+.9856474*time
      mnlong=mod(mnlong,360.)
      if(mnlong.lt.0.)mnlong=mnlong+360.

!   mean anomaly in radians between 0 and 2*pi
      mnanom=357.528+.9856003*time
      mnanom=mod(mnanom,360.)
      if(mnanom.lt.0.)mnanom=mnanom+360.
      mnanom=mnanom*rad

!   compute the ecliptic longitude and obliquity of ecliptic in radians
      eclong=mnlong+1.915*sin(mnanom)+.020*sin(2.*mnanom)
      eclong=mod(eclong,360.)
      if (eclong.lt.0.) eclong=eclong+360.
      oblqec=23.439-.0000004*time
      eclong=eclong*rad
      oblqec=oblqec*rad

!   calculate right ascension and declination
      num=cos(oblqec)*sin(eclong)
      den=cos(eclong)
      ra=atan(num/den)
!   force ra between 0 and 2*pi
      if (den.lt.0) then
          ra=ra+pi
      elseif (num.lt.0) then
          ra=ra+twopi
      endif

!   dec in radians
      dec=asin(sin(oblqec)*sin(eclong))

!   calculate Greenwich mean sidereal time in hours
      gmst=6.697375+.0657098242*time+hour 

!   hour not changed to sidereal time since 'time' includes the fractional day 
      gmst = mod(gmst,24.)
      if(gmst.lt.0.) gmst=gmst+24.

!   calculate local mean sidereal time in radians 
      lmst=gmst+long/15.
      lmst=mod(lmst,24.)
      if(lmst.lt.0.) lmst=lmst+24.
      lmst=lmst*15.*rad

!   calculate hour angle in radians between -pi and pi
      ha=lmst-ra
      if(ha.lt.-pi) ha=ha+twopi
      if(ha.gt.pi) ha=ha-twopi
!   change latitude to radians
      latrad=lat*rad

!   calculate azimuth and elevation
      el=asin(sin(dec)*sin(latrad)+cos(dec)*cos(latrad)*cos(ha))
      az=asin(-cos(dec)*sin(ha)/cos(el))

!!   this puts azimuth between 0 and 2*pi radians
      if (sin(dec)-sin(el)*sin(latrad).ge.0.) then
		if(sin(az).lt.0.) az=az+twopi
      else
      az=pi-az
      endif
!   if az=90 degs, elcritical=asin(sin(dec)/sin(latrad))
    elc=asin(sin(dec)/sin(latrad))
    if(el.ge.elc)az=pi-az
    if(el.le.elc.and.ha.gt.0.)az=twopi+az

!   calculate refraction correction for US stan. atmosphere
!   need to have el in degs before calculating correction
      el=el/rad

      if(el.ge.19.225) then 
         refrac=.00452*3.51823/tan(el*rad)
      else if (el.gt.-.766.and.el.lt.19.225) then
         refrac=3.51823*(.1594+.0196*el+.00002*el**2)/(1.+.505*el+.0845*el**2)
      else if (el.le.-.766) then
         refrac=0.0
      end if

!   note that 3.51823=1013.25 mb/288 C
      el=el+refrac
!   elevation in degs
!   calculate distance to sun in A.U. & diameter in degs
      soldst=1.00014-.01671*cos(mnanom)-.00014*cos(2.*mnanom)
      soldia=.5332/soldst

!   convert az and lat to degs before returning
    az=az/rad
	ha=ha/rad
	dec=dec/rad

!   mnlong in degs, gmst in hours, jd in days if 2.4e6 added;
!   mnanom,eclong,oblqec,ra,and lmst in radians
End subroutine

end program TestSunae
