      program ptoradian
c-----
c     This program considers at range of teleseismic
c     ray parameters typically observed. It then 
c     converts the ray parameters to angles in radians
c     for use with the Cerveny code. In addiiton
c     it gives the horizontal offset for a given source
c     depth.
c
c     It is assumed tha tht source depth is much greater than 
c     the thickness of the structure being modeled.
c-----
c     define the velocity of the P wave at
c     the source depth
c-----
      vel=8.5
c-----
c     define the source depth deemed sufficient for a plane-wave
c     approximation
c-----
      depth = 16000.0
      call bdoit(0.045,0.055,vel,depth,.true.)
      call bdoit(0.055,0.065,vel,depth,.false.)
      call bdoit(0.065,0.075,vel,depth,.false.)
      end

      subroutine bdoit(plow,phgh,vel,depth,lprinthead)
      real plow, phgh, vel, depth
      logical lprinthead
      real anglow, anghgh, angmid, thetamid
      real xoff
      real pi2

      if(lprinthead)then
        write(6,1)depth,vel
      endif
    1 format('    Mapping of ray parameter range to horizontal offset ',
     1'and ray angles in radians'/
     2'    for source depth of',f10.2,' and velocity of ',f10.3,' km/s'/
     3'  p range (s/km)  Theta   Xoffset            Forward Rays (+x) ',
     4'             Reverse Rays (-x)')

      pi2 = 3.1415927/2.0
      call getang(plow,vel,anglow)
      call getang(phgh,vel,anghgh)
      angmid = 0.5*(anghgh + anglow)
      thetamid = angmid * 180.0 / 3.1415927
c-----
c     the horizontal offset is x=depth*tan(angmid)
c-----
      xoff = depth*tan(angmid)
c-----
c     get the angles in radians for the forward and reverse
c     profiles
c     The Cerveny code measures angles with respect to the
c     horizontal. Upward rays are negative. 
c     Rays in the +x direction, the forward direction, will
c     have angles in the range [0, -pi/2]. 
c     Rays in the -x direction, the reverse direction, will 
c     have angles oin the range [-pi/2, -pi]
c-----
      flow = -pi2 + anghgh
      fhgh = -pi2 + anglow
      rlow = -pi2 - anglow
      rhgh = -pi2 - anghgh
      write(6,2)plow,phgh,thetamid,xoff,flow,fhgh,rlow,rhgh 
    2 format('[',f6.3,' to',f6.3,']',f7.3,f10.3,5x,
     1      2('[',f10.4,'  to',f10.4,']',5x))
      return
      end

      subroutine getang(p,vel,ang)
      real p, vel, ang
      real pi2
      pi2 = 3.1415927/2.0
      ang = asin(p*vel)
      return
      end
     

