 1000       continue
       read(5,*,end=9999)xlat,xlon,xmag,stk,dip,rake
       call dofmplt(xlat,xlon,xmag,dip,stk,rake,3)
       go to 1000
 9999       continue
       end

       subroutine dofmplt(xlat,xlon,xmag,dip,stk,rake,ityp)
c-----
c       ityp       I*4       - 0 plot mechanism
c                       1 plot T axis as a bar
c                       2 Plot P axis as a bar
c                       3 Plot regional max compressive stress as per
c                            Zoback
c-----


       parameter (mp=10,np=10)
       real mom(3,3)
       real*8 xmt(np,np),ev(np),ev1(np)
       logical hemis

       hemis = .true.

       call tensor(stk,dip,rake,1.0,mom)
c-----
c       force symmetry of moment tensor - this is valid for
c       double couples. 
c       If the field of a single couple is desired, do not force
c       this symmetry
c-----
       mom(3,1) = mom(1,3)
       mom(3,2) = mom(2,3)
       mom(2,1) = mom(1,2)

c-----
c       compute eigenvalues and eigenvectors of moment tensor matrix
c       Get index of largest and smallest eigenvalues
c-----
       do 100 i=1,3
              do 110 j=1,3
                     xmt(i,j) = dble(mom(i,j))
  110              continue
  100       continue
       call tred2(xmt,3,np,ev,ev1)
       call tqli(ev,ev1,3,np,xmt)
c-----
c       get the largest and smallest eigenvectors
c-----
       elg = -1.0e+38
       esm =  1.0e+38
       ilg = 1
       ism = 1
C       write(6,*)' EIGENVALUE, AND EIGENVECTOR OF M(i,j)'
C    2       format(' ',e11.4,'  (',e11.4,',',e11.4,',',e11.4,')')
       do 120 j=1,3
C              write(6,2)ev(j),(xmt(i,j),i=1,3)
              if(ev(j).gt.elg)then
                     elg = ev(j)
                     ilg = j
              endif
              if(ev(j).lt.esm)then
                     esm = ev(j)
                     ism = j
              endif
  120       continue
c-----
c       if(ityp.ne.0)just output the mechanism as the P or T axis, but
c       temper everything with the apparent dip.
c-----
       if(ityp.ge.1 .and. ityp.le.2)then
              if(ityp.eq.1)then
c-----
c                     T axis
c-----
                     j = ilg
              else if(ityp.eq.2)then
c-----
c                     P axis
c-----
                     j = ism
              endif
c-----
c       the projection on the horizontal axis is just sqrt(1-Z**2)
c-----
c-----
c       the apparent azimuth is atan2(ev(1),ev(2)
c-----
              call gttrpl(ev(j),xmt(1,j),xmt(2,j),xmt(3,j),
     1                hemis,trn,pln)
              return
       else if(ityp.eq.3)then
c-----
c              we must define identify the B axis
c-----
              ipaxis = ism
              itaxis = ilg
              if(ipaxis.eq.1)then
                     if(itaxis.eq.2)then
                            ibaxis = 3
                     else
                            ibaxis = 2
                     endif
              else if(ipaxis.eq.2)then
                     if(itaxis.eq.1)then
                            ibaxis = 3
                     else
                            ibaxis = 1
                     endif
              else if(ipaxis.eq.3)then
                     if(itaxis.eq.1)then
                            ibaxis = 2
                     else
                            ibaxis = 1
                     endif
              endif
              j = ipaxis
              call gttrpl(ev(j),xmt(1,j),xmt(2,j),xmt(3,j),
     1                hemis,trnp,plnp)
       write(6,'(a,3f8.2,3f6.0,i5,2f7.0)')'P: ',
     1       xlat,xlon,xmag,dip,stk,rake,j,trnp,plnp
              j = itaxis
              call gttrpl(ev(j),xmt(1,j),xmt(2,j),xmt(3,j),
     1                hemis,trnt,plnt)
       write(6,'(a,3f8.2,3f6.0,i5,2f7.0)')'T: ',
     1       xlat,xlon,xmag,dip,stk,rake,j,trnt,plnt
              j = ibaxis
              call gttrpl(ev(j),xmt(1,j),xmt(2,j),xmt(3,j),
     1                hemis,trnb,plnb)
c-----
c              now use the Zoback, 1992) classification
c-----
              call boxzob(xlat,xlon,xmag,xc,yc,radc,trnp,
     1                plnp,trnt,plnt,trnb,plnb)
              return
                            
       endif
       return
       end

       subroutine gttrpl(ev,xmt1,xmt2,xmt3,hemis,tr,pl)
       real*8 ev, xmt1, xmt2, xmt3
       logical hemis
       real tr, pl
       
       real*8 r
c-----
c       plot pole on focal sphere
c-----
c       convert to trend and plunge, lower hemisphere
c-----
       degrad = 3.1415927/180.0
       if(xmt3 .lt. 0.0d+00)then
              xmt1 = - xmt1
              xmt2 = - xmt2
              xmt3 = - xmt3
       endif
       r = xmt1**2 + xmt2**2 + xmt3**2
       r = dsqrt(r)
       xmt1 = xmt1/r
       xmt2 = xmt2/r
       xmt3 = xmt3/r
       pl = 90.0 - sngl(dacos(xmt3))/degrad
       r = xmt1**2 + xmt2**2 
       r = dsqrt(r)
       if(r.le.1.0d-06)then
              tr = 0.0
       else
              tr = sngl(datan2(xmt2,xmt1))/degrad
       endif
       if(.not.hemis)then
              pl =  - pl
              tr = 180 + tr
       endif
       return
       end


       subroutine tensor(stk,dip,slip,xmom,m)
c-----
c       calculate moment tensor for a double couple mechanism
c-----
c       stk       - R*4        - strike, measured clockwise from north
c                       when looking down, e.g., east = 90
c       dip       - R*4       - dip, measured from horizontal when looking
c                       in direction of strike, fault dips down to right
c                       0 <= dip <= 90
c       slip       - R*4       - slip. Measured from horizontal with respect
c                       to strike. If -180 < slip <= 180 is taken as
c                       the convention, then a negative slip, implies 
c                       that the P-wave first motion in the center of the
c                       focal sphere is negative (e.g., dilatation)
c       xmom       - R*4       - moment value
c       m(3,3)       - R*4       - moment tensor
c-----
       integer LOT
       parameter (LOT=6)
       real stk, dip, slip, xmom, m(3,3)
       real degrad, tol, dp, st, sl
       real sinp, cosp, sint, cost, sinp2, cosp2, sinlp, coslp
       real sint2, cost2
       
              degrad=0.0174532925
              tol = 1.0e-7
              dp = degrad*dip
              st = degrad*stk
              sl = degrad*slip
              sinp=sin(dp)
              cosp=cos(dp)
              sinlp=sin(sl)
              coslp=cos(sl)
              sint=sin(st)
              cost=cos(st)
              sinp2=sin(2.*dp)
              cosp2=cos(2.*dp)
              sint2=sin(2.*st)
              cost2=cos(2.*st)
              m(1,1)=(-sinp*coslp*sint2-sinp2*sinlp*sint*sint)*xmom
              m(2,2)=(sinp*coslp*sint2-sinp2*sinlp*cost*cost)*xmom
              m(3,3)=(sinp2*sinlp)*xmom
              m(1,2)=(sinp*coslp*cost2+0.5*sinp2*sinlp*sint2)*xmom
              m(1,3)=(-cosp*coslp*cost-cosp2*sinlp*sint)*xmom
              m(2,3)=(-cosp*coslp*sint+cosp2*sinlp*cost)*xmom
              m(2,1) = m(1,2)
              m(3,1) = m(1,3)
              m(3,2) = m(2,3)
c-----
c       clean up small values
c-----
       xmax= -1.0e+38
       do 4 i=1,3
              do 5 j=1,3
                     if(abs(m(i,j)).gt.xmax)xmax = abs(m(i,j))
    5              continue
    4       continue
       thresh = tol * xmax
       do 6 i=1,3
              do 7 j=1,3
                     if(abs(m(i,j)).lt.thresh) m(i,j) = 0.0
    7              continue
    6       continue

c-----
c       write out the information
c-----
C       write(LOT,*)' DIP       =',dip
C       write(LOT,*)' SLIP      =',slip
C       write(LOT,*)' STRIKE    =',stk
C       write(LOT,*)' MOMENT    =',xmom
C       write(LOT,*)' MOMENT TENSOR:'
C       write(LOT,*)' M(1,1)    =',m(1,1)
C       write(LOT,*)' M(1,2)    =',m(1,2)
C       write(LOT,*)' M(1,3)    =',m(1,3)
C       write(LOT,*)' M(2,1)    =',m(2,1)
C       write(LOT,*)' M(2,2)    =',m(2,2)
C       write(LOT,*)' M(2,3)    =',m(2,3)
C       write(LOT,*)' M(3,1)    =',m(3,1)
C       write(LOT,*)' M(3,2)    =',m(3,2)
C       write(LOT,*)' M(3,3)    =',m(3,3)
       return
       end

      subroutine tred2(a,n,np,d,e)
      parameter (LER=0, LIN=5, LOT=6)
c----                                             
c---- Ref.: Numerical Recipes by Press, Flannery, Teukolsky,
c---- and Vetterling, Cambridge, 1987, 818pp
c---- EQN 11.2.xx p355
c----
c---- Householder reduction of a real, symmetric, N by N matrix A, stored
c---- in an NP by NP physical array. On output, A is replaced by the ortho-
c---- gonal matrix Q effecting the transformation. D returns the diagonal
c---- elements in the tridiagonal matrix, and E the off-diagonal elements
c---- with E(1)=0. 
      implicit real*8 (a-h,o-z)
      dimension a(np,np),d(np),e(np)
      if(n.gt.1) then
      do 18 i=n,2,-1
      l=i-1
      h=0.
      scale=0.
      if(l.gt.1) then
      do 11 k=1,l
      scale=scale+abs(a(i,k))
   11 continue
      if(scale.eq.0.) then
      e(i)=a(i,l)
      else
      do 12 k=1,l
      a(i,k)=a(i,k)/scale
      h=h+a(i,k)**2
   12 continue
      f=a(i,l)
      g=-sign(sqrt(h),f)
      e(i)=scale*g
      h=h-f*g
      a(i,l)=f-g
      f=0.
      do 15 j=1,l
      a(j,i)=a(i,j)/h
      g=0.
      do 13 k=1,j
      g=g+a(j,k)*a(i,k)
   13 continue
      if(l.gt.j) then
      do 14 k=j+1,l
      g=g+a(k,j)*a(i,k)
   14 continue
      endif
      e(j)=g/h
      f=f+e(j)*a(i,j)
   15 continue
      hh=f/(h+h)
      do 17 j=1,l
      f=a(i,j)
      g=e(j)-hh*f
      e(j)=g
      do 16 k=1,j
      a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
   16 continue
   17 continue
      endif
      else
      e(i)=a(i,l)
      endif
      d(i)=h
   18 continue
      endif
      d(1)=0.
      e(1)=0.
      do 23 i=1,n
      l=i-1
      if(d(i).ne.0.) then
      do 21 j=1,l
      g=0.
      do 19 k=1,l
      g=g+a(i,k)*a(k,j)
   19 continue
      do 20 k=1,l
      a(k,j)=a(k,j)-g*a(k,i)
   20 continue
   21 continue
      endif
      d(i)=a(i,i)
      a(i,i)=1.
      if(l.ge.1) then
      do 22 j=1,l
      a(i,j)=0.
      a(j,i)=0.
   22 continue
      endif
   23 continue
      return
      end
 
      subroutine tqli(d,e,n,np,z)
      parameter (LER=0, LIN=5, LOT=6)
c----                                             
c---- Ref.: Numerical Recipes by Press, Flannery, Teukolsky,
c---- and Vetterling, Cambridge, 1987, 818pp
c---- EQN 11.2.xx p362
c----
c---- QL algorithm with implicit shifts, to determine the eigenvalues and
c---- eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
c---- symmetric matrix previously reduced by TRED2. D is a vector of length
c---- NP. On input, its first N elements are the diagonal elements of the
c---- tridiagonal matrix. On output, it returns the eigenvalues. The vector
c---- E inputs the subdiagonal elements of the tridiagonal matrix, with E(1)
c---- arbitrary. On output, E is destroyed. The N by N matrix Z is input as
c---- the identity matrix. If a matrix has been reduced by TRED2, then Z is
c---- input as the matrix output by TRED2. On output, the k-th column of Z
c---- returns the normalized eigenvector corresponding to D(k). 
      implicit real*8 (a-h,o-z)
      dimension z(np,np),d(np),e(np)
      if(n.gt.1) then
      do 11 i=2,n
      e(i-1)=e(i)
   11 continue
      e(n)=0.
      do 15 l=1,n
      iter=0
    1 do 12 m=l,n-1
      dd=abs(d(m))+abs(d(m+1))
      if(abs(e(m))+dd.eq.dd) goto 2
   12 continue
      m=n
    2 if(m.ne.l) then
      if(iter.eq.30) write(LOT,*) 'too many iterations'
      iter =iter+1
      g=(d(l+1)-d(l))/(2.*e(l))
      r=sqrt(g**2+1.)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1.
      c=1.
      p=0.
      do 14 i=m-1,l,-1
      f=s*e(i)
      b=c*e(i)
      if(abs(f).ge.abs(g)) then
      c=g/f
      r=sqrt(c**2+1.)
      e(i+1)=f*r
      s=1./r
      c=c*s
      else
      s=f/g
      r=sqrt(s**2+1.)
      e(i+1)=g*r
      c=1./r
      s=s*c
      endif
      g=d(i+1)-p
      r=(d(i)-g)*s+2.*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
      do 13 k=1,n
      f=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*f
      z(k,i)=c*z(k,i)-s*f
   13 continue
   14 continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.
      goto 1
      endif
   15 continue
      endif
      return
      end

       subroutine boxzob(xlat,xlon,xmag,xc,yc,radc,trnp,
     1                plnp,trnt,plnt,trnb,plnb)
       character regime*2
       integer red, green, blue
c-----
c       first determine the type of mechanism according to Zoback 1992)
c       dir = direction of maximum compressive stress
c-----
       red = 1000
       green = 1050
       blue = 1100
       dir = 0.0
       kolor1 = 1
       kolor2 = 1
       regime = 'XX'
C  /* Regime characterization based on WSM - Zoback */
C  if ( (Ppbt->p.dip >= 52.0) && (Ppbt->t.dip <= 35.0) )  {
C     strcpy( regime, "NF");
C     *azi = (int) ( AZIMUTH(Ppbt->b.strike) + 0.5);
       if(plnp .ge. 52 .and. plnt .le. 35.0)then
              regime = 'NF'
              dir = trnb
              kolor1 = red
              kolor2 = red
C  } else if ( (40.0 <= Ppbt->p.dip) && (Ppbt->p.dip <= 52.0) && (Ppbt->t.dip <= 20.0) ) {
C     strcpy( regime, "NS");
C     *azi = (int) ( AZIMUTH(Ppbt->t.strike) + 90.5);
       else if(plnp .ge. 40 .and. plnp .le. 52.0 .and. plnt.le.20.0)then
              regime = 'NS'
              dir = trnt +90
              kolor1 = red
              kolor2 = green
C  } else if ( (40.0 < Ppbt->p.dip) && (Ppbt->b.dip >= 45.0) && (Ppbt->t.dip <= 20.0) )  {
C     strcpy( regime, "SS");
C     *azi = (int) ( AZIMUTH(Ppbt->t.strike) + 90.5);
       else if(plnp.lt.40. .and. plnb.ge.45 .and. plnt.le.22)then
              regime = 'SS'
              dir = trnt + 90.
              kolor1 = green
              kolor2 = green
C  } else if ( (Ppbt->p.dip <= 20.0) && (Ppbt->b.dip >= 45.0) && (Ppbt->t.dip < 40.0) )  {
C     strcpy( regime, "SS");
C     *azi = (int) ( AZIMUTH(Ppbt->p.strike) + 0.5);
       else if(plnp.le.22. .and. plnb.ge.45 .and. plnt.le.40)then
              regime = 'SS'
              dir = trnp 
              kolor1 = green
              kolor2 = green
C  } else if ( (Ppbt->p.dip <= 20.0) && (40.0 <= Ppbt->t.dip) && (Ppbt->t.dip <= 52.0) )  {
C     strcpy( regime, "TS");
C     *azi = (int) ( AZIMUTH(Ppbt->p.strike) + 0.5);
       else if(plnp.le.20. .and. plnt.ge.40. .and. plnt.le.52.)then
              regime = 'TS'
              dir = trnp
              kolor1 = blue
              kolor2 = green
C  } else if ( (Ppbt->p.dip <= 35.0) && (Ppbt->t.dip >= 52.0) )  {
C     strcpy( regime, "TF");
C     *azi = (int) ( AZIMUTH(Ppbt->p.strike) + 0.5);
       else if(plnp.le.35. .and. plnt.ge.52.)then
              regime = 'TF'
              dir = trnp
              kolor1 = blue
              kolor2 = blue
       endif
       write(6,'(a,1x,a2,2x,4f10.2,2i6)' )'Regime:',
     1                regime,xlat,xlon,xmag,dir,kolor1,kolor2
       
       return
       end

         
