       program strainpulse96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME III                                                     c
c                                                                     c
c      PROGRAM: STRAINPULSE96                                              c
c                                                                     c
c      COPYRIGHT 2021                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       06 SEP 2000 - build in P, S  and S  first arrival times
c       14 JAN 2001 - use new FILE96.2 with sa, sc, sf, sl, sn, sr
c       19 FEB 2001 - fixed obscure error having that messed up
c               spectra in fundamental mode only at long period
c               by forcing kkf = 4 for all runs. Also
c               redefined mode =0 for fundamental in -M mode
c       30 JAN 2002 - caught error in excitr for ZDD had 
c            d1=w1*v1*fact1 should
c               be d1=d1*v1*fact1
c       17 OCT 2002 - Added description of dfile format to usage routine
c       10 SEP 2003 - error in pressure field in fluid due to a
c               point horizontal force
c           BAD  cvp(6) = cvp(6) + dkp(6)*xmom*cmplx(st2,-ct2)
c           GOOD cvp(6) = cvp(6) + dkp(6)*xmom*cmplx(st1,-ct1)
c           also changes notation slightly dk, dkk, dkp in excitr
c               to dkz, dkr , dkp
c       05 FEB 2004 - modified to be slightly more tolerant about DT for
c           user supplied pulse - now issues WARNING and not termination
c           mlarocca@ov.ingv.it
c       07 FEB 2005 - add a -Z flag to indicate that the 
c            internal parabolic 
c           or triangular pulses are to be zero phase 
c       04 AUG 2006 - corrected error in first arrival pick that
c               falsely gave the refraction time instead of the
c               direct time because the refraction arrival was
c               unphysical
c       12 DEC 2006 - set header values evlat, evlon, stlat, stlon to -12345
c           for compatibility of resulting SAC traces
c       26 MAY 2007 - careful reworking to make the slat2d work
c          FIXED exceed array bound problem with -2 flat. Everything
c          now defaults to NSAMP as the maximum of anything
c       21 JAN 2008 - upgraded travel time computation so that all code
c          is the same for CPS isotropic media
c       22 JAN 2008 - rearranged routines, put in new common
c                     isotropic travel time routines from time96
c       25 JAN 2008 - put Radius of Earth into common/earth/radius for
c                      generality
c                   -  define a separate common block for the
c                      S  velocity and density
c                   -  have sphericity correction
c                      work on common blocks instead of procedure call
c                   -  create a default adomod  to fill the S  for a flat model
c                      note the separation of S  is important for wavenumber
c                      integration code
c       08 FEB 2008 -  subtle change in fstarr for source receiver in same layer -
c                      spherical mapping was not done
c       18 FEB 2009 -  caught egregious error in frstar where refdep was not set
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       01 JUN 2013 -  Modified subroutine to prevent NaN for direct ray by replacing
c                      pnew = 0.999d+00 * pupper to
c                      pnew = 0.99d+00 * pupper
c                      also modified subroutine frstar to have the dogeom argument.
c                      We do not want to compute teleseismic geometrical spreading
c       22 JUL 2013    set NLAY=200 from NLAY=100 to be compatible 
c                      model format
c       06 JUL 2017    Corrected the lines in excitr
c            from
c                 v1 = dsqrt(dble(ares(j,1))*dble(arer(j,1)))
c              1      /dsqrt(dble(wvmr(j,1))*rx)
c                 v2 = dsqrt(dble(ares(j,2))*dble(arer(j,2)))
c              1      /dsqrt(dble(wvmr(j,2))*rx)
c                 to
c                 v1 = dsqrt(dble(ares(j,1))*dble(arer(j,1)))
c              1      /dsqrt(dble(wvmrr(j,1))*rx)
c                 v2 = dsqrt(dble(ares(j,2))*dble(arer(j,2)))
c              1      /dsqrt(dble(wvmrr(j,2))*rx)
c                 and  in excitl
c                 from
c                     v1 = dsqrt(dble(ales(j,1))*dble(aler(j,1)))
c              1          /dsqrt(dble(wvml(j,1))*rx)
c                     v2 = dsqrt(dble(ales(j,2))*dble(aler(j,2)))
c              1          /dsqrt(dble(wvml(j,2))*rx)
c                 to
c                     v1 = dsqrt(dble(ales(j,1))*dble(aler(j,1)))
c              1          /dsqrt(dble(wvmlr(j,1))*rx)
c                     v2 = dsqrt(dble(ales(j,2))*dble(aler(j,2))c
c              1          /dsqrt(dble(wvmlr(j,2))*rx)
c            To agree with Levshin's theory that is used with slat2d
c       03 JAN 2021    modified spulse96 to create sac files of strain 
c             for a given source. The array xx(21) is replaced by xx(63)
c             xx i= 1,21 are the components of the Green's functions
c                 =20,42 are the far-field d/dr of the 1,21 values
c                        obtained by multiplying by ( -i k )
c                 =43,63 are d/dz of the Grene's functions at the
c                        receiver
c             
c             One problem with computing strain is that we must be 
c             very careful about the units. Since we are creating the
c             synthetics and the strain from surface-wave modal superposition,
c             we will need terms such as dU/dr and dU/dz at the receiver.
c             
c             In the far-field expansion the dUdr ~ - i k U. To get the
c             dU/dz we need to use the stress definition, e.g.,
c             
c             TT = mu dUT/dz
c             TR = mu ( d UR/dz + k UZ)
c             TZ = (lambda + 2 mu) d UZ/dz - k lambda UR
c
c             or
c             d UT/dz = TT/my
c             d UR/dz = ( TR - k UZ)/mu
c             d UZ/dz = ( TZ + k lambda UZ)/(lambda + 2 mu)
c             
c             Note these terms use the units of the model. Typically the model 
c             is given as KM KM/S GM/CM^3`
c
c             Because of this we will have to be careful about the strains, e.g.,
c             The units will be DIsplacement (which is another question)/KM
c
c             So for a moment of 1.0e+20 the output will be cm so the strain will be
c             multiplied 0.01 / 1000  = 1.0e+6
c
c            Note the original spulse96 was defined to have Uz positive up
c            we must be careful of this sign change for the comptuation of strain
c            wbich uses a coordinate system of  z positive down
c
c            Logical units used for IO
c            1  - Love wave eigenfunciton file
c            2  - Rayleigh wave eigenfunciton file
c-----

c-----
c       This program takes the output of sregn96 and slegn96
c       and makes surface wave time series in the file96(V) format
c
c       PROGRAM CONTROL
c           INPUT  - from command line
c           OUTPUT - standard output
c-----
c       spulse96 -d Distance_File -v  -t -o -p -i -a alpha \
c           -l dur  -D -V -A -F rfile -m mult \
c       
c        -d Distance_File Distance control information
c        -v           Verbose output 
c        -t           Triangular pulse of base 2 dur dt 
c        -p           Parabolic Pulse of base  4 dur dt 
c        -o           Ohnaka pulse with parameter alpha 
c        -i           Dirac Delta function 
c        -a alpha     Shape parameter for Ohnaka pulse 
c        -D           Output is ground displacment 
c        -V           Output is ground velocity (default) 
c        -A           Output is ground acceleration 
c        -F rfile     User supplied pulse 
c        -m mult      Multiplier (default 1.0) 
c        -?           Write this help message 
c        -h           Write this help message 
c               
c-----
        parameter (LER=0, LIN=5, LOT=6)
        character*80 dfile, ryfile, lvfile, rfile
        integer*4 ntau, ipt, idva
        real*4 xmult, alp
        logical dolat, dodble, dolock, dozero, dotest1,dostep

        common/stresstrain/dostrain,dostress,
     1     az, baz, xmt(3,3), xmom,forcex, forcey,forcez,
     1     lrot, caz, saz, c2az, s2az
        logical dostrain, dostress
        real az, baz, xmt, xmom,forcex, forcey,forcez
        real caz, saz, c2az, s2az
        logical lrot
c-----
c       dostrain  - L .true. output strain
c       dostress  - L .true. outpu stress
c-----
        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt
c-----
c       at the receiver position
c       l2mu = rho Vp^2, mu = rho Vs^2
c-----

        character ostr*120

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(ntau,ipt,alp,dfile,idva,
     1      xmult,rfile,ostr,dolat,dodble,nmode,dolock,
     1      dozero,
     1      dostrain,dostress,az,baz,lrot,xmom,xmt,forcex,forcey,forcez,
     1      dotest1,dostep,outfmt)
        WRITE(0,*)'ntau              :',ntau
        WRITE(0,*)'ipt               :',ipt 
        WRITE(0,*)'alp               :',alp 
        WRITE(0,*)'dfile             :',dfile
        WRITE(0,*)'idva              :',idva
        WRITE(0,*)'xmult             :',xmult
        WRITE(0,*)'rfile             :',ntau
        WRITE(0,*)'ostr              :',ostr
        WRITE(0,*)'dolat             :',dolat
        WRITE(0,*)'dodble            :',dodble
        WRITE(0,*)'nmode             :',nmode
        WRITE(0,*)'dolock            :',dolock
        WRITE(0,*)'dozero            :',dozero
        WRITE(0,*)'dostrain          :',dostrain
        WRITE(0,*)'dostress          :',dostress
        WRITE(0,*)'dotest1           :',dotest1
        WRITE(0,*)'az                :',az
        WRITE(0,*)'baz               :',baz
        WRITE(0,*)'lrot              :',lrot
        WRITE(0,*)'xmom              :',xmom
        WRITE(0,*)'xmt               :',xmt(1,1),xmt(1,2),xmt(1,3)
        WRITE(0,*)'                  :',xmt(2,1),xmt(2,2),xmt(2,3)
        WRITE(0,*)'                  :',xmt(3,1),xmt(3,2),xmt(3,3)
        WRITE(0,*)'forcex            :',forcex
        WRITE(0,*)'forcey            :',forcey
        WRITE(0,*)'forcez            :',forcez
        WRITE(0,*)'dostep            :',dostep
        WRITE(0,*)'outfmt            :',outfmt
c-----
c       get the azimuths for the source
c-----
        degrad = 3.1415927/180.0
        caz = cos(degrad*az)
        saz = sin(degrad*az)
        c2az = cos(2*degrad*az)
        s2az = sin(2*degrad*az)
c-----
c       define the eigenfunction files
c-----
        if(dolat)then
            lvfile = 'slatl96.egn'
            ryfile = 'slatr96.egn'
        else
            lvfile = 'slegn96.egn'
            ryfile = 'sregn96.egn'
        endif
c-----
c       process to make time series
c-----
        call process(ntau,ipt,alp,dfile,idva,
     1      xmult,rfile,ryfile,lvfile,ostr,
     2      mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat,
     3      dotest1,dostep,dostrain,dostress)
        end

        subroutine adomod()
c-----
c       just fill the rhosh, bsh and qbsh arrays 
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        do  i=1,mmax
            bsh(i)=b(i)
            qbsh(i)=qb(i)
            rhosh(i) = rho(i) 
        enddo
        return
        end

        subroutine adosph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c-----
c       mmax    I*4 number of layers
c       ipsvsh  I*4     1 - get P time
c                       2 - get S  time
c                       3 - get S  time
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        double precision z0,z1,r0,r1,dr,ar,tmp
        integer ifunc

        common/earth/radius
        real radius

        ar=radius
        dr=0.0d0
        r0=ar + refdep
        d(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(d(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            d(i)=z1-z0
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)

            a(i)=a(i)*tmp
            b(i)=b(i)*tmp
            bsh(i)=b(i)
            qbsh(i)=qb(i)
            rhosph=rho(i)
            rhosh(i) = rhosph * tmp **(-5.0)
            rho(i) = rhosph * tmp **(-2.275)
            r0 = r1
   10   continue
        d(mmax)=0.0
        return
        end

        subroutine deriv(y,n,dt)
c-----
c       Use centered difference to take derivative of time series
c-----
c       y   R*4 - time series to be differentiated and returned as y
c       n   I*4 - length of time series
c       dt  R*4 - sample interval
c-----
        real*4 y(1)
        b0=1.0/dt
        xm1=y(1)
        xm2=0.
        y1=0.
        do 11 i=1,n
            y2=b0*(y(i)-xm1) 
            xm2=xm1
            xm1=y(i)
            y(i)=y2
   11   continue
        return
        end

        subroutine excitl(per,mode,rr,tshft,xx,ipar)
c-----
c     This generates the T component of seismograms
c     for  1: dip-slip source     2: strike-slip source.
c
c     vt(i,j)   i=1,2  real & imag part    j=1,2  source
c---
        integer ipar(10)
        parameter (NGRN=63)
        complex xx(NGRN)
        parameter (NMD=100)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),tt0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr0,si0,xmom0,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)
        common/receiver/l2mu,mu
        real l2mu, mu
c-----
c       at the receiver position
c       l2mu = rho Vp^2, mu = rho Vs^2
c-----

        real*8 rx,tx,rat,wvno,atn1,atn2,fact1,fact2
        real*8 v1,v2,w1,w2,t1,ct1,st1,sr,si,xmom,pi,pi4
        real*8 ww1, ww2
        real*8 dk(3),omega,ddk(3)
        complex*16 cvt(3+3+3)
        integer*4 kkr,kkl,kkf,np0
        complex*16 kfac
        data pi/3.141592653589793d+00/
c
        xx(5) = cmplx(0.0,0.0)
        xx(8) = cmplx(0.0,0.0)
        xx(15) = cmplx(0.0,0.0)
c-----
c       if source or receiver is in fluid, there is no 
c            far-field S  motion
c-----
        if(ipar(2).eq.1 .or. ipar(3).eq.1)return
        xmom=dble(xmom0)/dsqrt(2.d+00*pi)
        pi4=pi/4.d+00
        do 90 j=1,3+3+3
            cvt(j)=dcmplx(0.0d+00,0.0d+00)
   90   continue
c-----
c       initialize processing
c-----

            j1=mods
            j2=modt
            if(j2.gt.mode) j2=mode
            j3=modi

        if(j1.gt.j2) return
        omega=2.d+00*pi/dble(per)
        per1=pers1(1)
        per2=pers2(1)
        n1=0
        n2=0
        ii=0
  100   continue
        ii=ii+1
        call finds(1,ii,j1,j2,j11,j22,j33,n1,n2,n3,mode)
        icnt=0
        do 300 j=j11,j22
            if(ales(j,1).lt.1.e-30.and.ales(j,2).lt.1.e-30 .and.
     1      aler(j,1).lt.1.e-30.and.aler(j,2).lt.1.e-30) then
                icnt=icnt+1
            else
                icnt=0
            endif
            if(icnt.gt.3) go to 310
            if(ales(j,1).lt.1.e-30.or.ales(j,2).lt.1.e-30 .and.
     1      aler(j,1).lt.1.e-30.or.aler(j,2).lt.1.e-30) go to 300
            k=(n3-1)*NMD+j
            if(k.ge.j1.and.k.le.j2) then
                if(mod(k-j1,j3).eq.0) go to 150
            endif
            go to 300
  150       continue
            rat=0.0d+00
            rat=(dble(per)-dble(per1))/(dble(per2)-dble(per1))
            wvno=(dble(wvml(j,1)*per1) + 
     1          dble(wvml(j,2)*per2-wvml(j,1)*per1)*rat)/per
            kfac = dcmplx(0.0d+00,-wvno)
c-----
c       process this distance
c-----
            rx=dble(rr)
            tx=dble(tshft)
            if(zeta.lt.0.0) then
                fact1=1.d+00
                fact2=1.d+00
                if(fact2.eq.1.d+00) go to 160
            endif
            atn1=1.0d+00
            if(zeta.gt.2.0) then
                atn1=zeta/100.0
            else
                if(per.lt.10.0) then
                    atn1=dble(per)**dble(zeta)
                else
                    atn1=dble(10.0)**dble(zeta)
                endif
            endif
            atn2=atn1*dble(gaml(j,2))*rx
            atn1=atn1*dble(gaml(j,1))*rx
            fact1=0.0d+00
            if(atn1.lt.80.0) fact1=1./dexp(atn1)
            fact2=0.0d+00
            if(atn2.lt.80.0) fact2=1./dexp(atn2)
  160       continue
c------
c       DS component.
c------
            
            v1 = dsqrt(dble(ales(j,1))*dble(aler(j,1)))
     1          /dsqrt(dble(wvmlr(j,1))*rx)
            v2 = dsqrt(dble(ales(j,2))*dble(aler(j,2)))
     1          /dsqrt(dble(wvmlr(j,2))*rx)
            w1 = dble(dut(j,1))*v1*fact1
            w2 = dble(dut(j,2))*v2*fact2
            ww1= dble(ut0(j,1))*w1
            ww2= dble(ut0(j,2))*w2
            dk(1)  = ww1+(ww2-ww1)*rat
            ww1= dble(tt0(j,1))*w1
            ww2= dble(tt0(j,2))*w2
            ddk(1) = ww1+(ww2-ww1)*rat
            ddk(1) = ddk(1) /mu
c------
c       SS component.
c------
            w1 = dble(wvmls(j,1))*dble(ut(j,1))
            w1 = w1*v1*fact1
            w2 = dble(wvmls(j,2))*dble(ut(j,2))
            w2 = w2*v2*fact2
            ww1= dble(ut0(j,1))*w1
            ww2= dble(ut0(j,2))*w2
            dk(2)  = ww1+(ww2-ww1)*rat
            ww1= dble(tt0(j,1))*w1
            ww2= dble(tt0(j,2))*w2
            ddk(2) = ww1+(ww2-ww1)*rat
            ddk(2) = ddk(2) /mu
c-----
c       HF Component
c-----
            w1 = dble(ut(j,1))*v1*fact1
            w2 = dble(ut(j,2))*v2*fact2
            ww1= dble(ut0(j,1))*w1
            ww2= dble(ut0(j,2))*w2
            dk(3)  = ww1+(ww2-ww1)*rat
            ww1= dble(tt0(j,1))*w1
            ww2= dble(tt0(j,2))*w2
            ddk(3) = ww1+(ww2-ww1)*rat
            ddk(3) = ddk(3) /mu
c-----
c      account for units up to here dx = km rho=gm/cm^3 vel=km/s
c-----
C           ddk(1) = 1.0e-3 *ddk(1)
C           ddk(2) = 1.0e-3 *ddk(2)
C           ddk(3) = 1.0e-3 *ddk(3)
c-----
c       Introduce time shift and also phase shift from eigenfunctions
c-----
            t1=omega*tx-wvno*rx+pi4
            ct1=dcos(t1)
            st1=dsin(t1)
c-----
c       TDS - requires -i in front of eigenfunction
c          note that these are in order Ut, d Ut/dr and d Ut/dz
c-----
            cvt(1)       = cvt(1)     + dk(1) *xmom*cmplx(st1,-ct1)
            cvt(1+3)     = cvt(1+3)   + dk(1) *xmom*cmplx(st1,-ct1)*kfac
            cvt(1+3+3)   = cvt(1+3+3) + ddk(1)*xmom*cmplx(st1,-ct1)
c-----
c       TSS - eigenfunction excitation is pure real
c          note that these are in order Ut, d Ut/dr and d Ut/dz
c-----
            cvt(2)     = cvt(2)     + dk(2) *xmom*cmplx(ct1,st1)
            cvt(2+3)   = cvt(2+3)   + dk(2) *xmom*cmplx(ct1,st1)*kfac
            cvt(2+3+3) = cvt(2+3+3) + ddk(2)*xmom*cmplx(ct1,st1)
c-----
c       THF - requires -i in front of eigenfunction
c          note that these are in order Ut, d Ut/dr and d Ut/dz
c-----
            cvt(3)     = cvt(3)     + dk(3) *xmom*cmplx(st1,-ct1)
            cvt(3+3)   = cvt(3+3)   + dk(3) *xmom*cmplx(st1,-ct1)*kfac
            cvt(3+3+3) = cvt(3+3+3) + ddk(3)*xmom*cmplx(st1,-ct1)
  300   continue
        if(j2.gt.j33) go to 100
  310   continue
c-----
c       introduce source time function
c-----
        if(ms.eq.1)then
            sr=1.0d+00
            si=0.0d+00
        else
            sr=dble(sr0)
            si=dble(si0)
        endif
        do 400 j=1,3+3+3
            cvt(j) = cvt(j)*dcmplx(sr,si)
  400   continue
        do 501 j=1,3+3+3
            t1 = cdabs(cvt(j))
            if(t1 .le. 1.0d-20)cvt(j) = 0.0d+00
  501   continue
c-----
c       the comman arises from the source excitation
c       5 - TDS
c       8 - TSS
c       15 - THF
c-----
        xx( 5)    = -cvt(1)
        xx( 8)    = -cvt(2)
        xx(15)    = -cvt(3)
        xx( 5+21) = -cvt(1+3)
        xx( 8+21) = -cvt(2+3)
        xx(15+21) = -cvt(3+3)
        xx( 5+42) = -cvt(1+3+3)
        xx( 8+42) = -cvt(2+3+3)
        xx(15+42) = -cvt(3+3+3)
        return
        end

        subroutine excitr(per,mode,rr,tshft,xx,ipar)
c-----
c
c     This generates the Z and R components of
c     seismogram for
c      1: 45-deg dip-slip source  2: strike-slip source
c      3: dip-slip source         4: explosion source.
c
c     number of mode is unlimited.
c
        integer ipar(10)
        parameter (NGRN=63)
        parameter (NMD=100)
        real*8 rat,atn1,atn2,fact1,fact2,v1,v2,w1,w2,u1,u2
        real*8 uu1, ww1, pp1
        real*8 uu2, ww2, pp2
        real*8 d1, d2
        real*8 p1,p2
        real*8 t1,t2,ct1,ct2,st1,st2,sr,si,xmom,rx,tx
        real*8 wvno,omega,pi,pi4,pi34
        real*8  dkz(6), dkr(6), dkp(6)
        real*8 ddkz(6),ddkr(6),ddkp(6)
        complex*16 cvz(6+6+6), cvr(6+6+6), cvp(6+6+6)
        integer*4 kkr,kkl,kkf,np0,np
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/ctrl/   sr0,si0,xmom0,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)

        common/water/rrho

        common/receiver/l2mu,mu
        real l2mu, mu
c-----
c       at the receiver position
c       l2mu = rho Vp^2, mu = rho Vs^2
c-----

        real duzdz, durdz
        complex xx(NGRN)
        complex*16 kfac
        data pi/3.141592653589793d+00/
c------
        xmom=dble(xmom0)/dsqrt(2.d+00*pi)
        pi4=pi/4.d+00
        pi34=3.d+00*pi/4.d+00
        do 91 j=1,6+6+6
            cvz(j)  = dcmplx(0.0d+00,0.0d+00)
            cvr(j)  = dcmplx(0.0d+00,0.0d+00)
            cvp(j)  = dcmplx(0.0d+00,0.0d+00)
   91   continue

        j1=mods
        j2=modt
        if(j2.gt.mode) j2=mode
        j3=modi

        if(j1.gt.j2) return
        omega=2.d+00*pi/dble(per)
        per1=pers1(2)
        per2=pers2(2)
        n1=0
        n2=0
        ii=0
  100   continue
        ii=ii+1
c------
c       find where the data locate.
c------
        call finds(2,ii,j1,j2,j11,j22,j33,n1,n2,n3,mode)
c------
c       put different elements into the final solution 
c            after interpolation.
c------
        icnt=0
        do 400 j=j11,j22
        if(ares(j,1).lt.1.e-30.and.ares(j,2).lt.1.0e-30 .and.
     1  arer(j,1).lt.1.e-30.and.arer(j,2).lt.1.0e-30) then
          icnt=icnt+1
        else
          icnt=0
        endif
        if(icnt.gt.3) go to 410
        if(ares(j,1).lt.1.e-30.or.ares(j,2).lt.1.0e-30 .and.
     1  arer(j,1).lt.1.e-30.or.arer(j,2).lt.1.0e-30) go to 400
        k=(n3-1)*NMD+j
        if(k.ge.j1.and.k.le.j2) then
          if(mod(k-j1,j3).eq.0) go to 150
        endif
        go to 400
  150   continue
        rat=0.0d+00
        rat=(dble(per)-dble(per1))/(dble(per2)-dble(per1))
        wvno=(dble(wvmr(j,1)*per1) + 
     1      dble(wvmr(j,2)*per2-wvmr(j,1)*per1)*rat)/per
        kfac = dcmplx(0.0d+00,-wvno)
c------
c       For different distance using same npt.
c------
        rx=dble(rr)
        tx=dble(tshft)
c
        if(zeta.lt.0.0) then
          fact1=1.d+00
          fact2=1.d+00
          if(fact2.eq.1.d+00) go to 160
        endif
        atn1=1.0d+00
        if(zeta.gt.2.0) then
            atn1=zeta/100.0
        else
            if(per.lt.10.0) then
                atn1=dble(per)**dble(zeta)
            else
                atn1=dble(10.0)**dble(zeta)
            endif
        endif
        atn2=atn1*dble(gamr(j,2))*rx
        atn1=atn1*dble(gamr(j,1))*rx
        fact1=0.0d+00
        if(atn1.lt.80.0) fact1=1./dexp(atn1)
        fact2=0.0d+00
        if(atn2.lt.80.0) fact2=1./dexp(atn2)
  160   continue
c
        v1 = dsqrt(dble(ares(j,1))*dble(arer(j,1)))
     1      /dsqrt(dble(wvmrr(j,1))*rx)
        v2 = dsqrt(dble(ares(j,2))*dble(arer(j,2)))
     1      /dsqrt(dble(wvmrr(j,2))*rx)
c------
c       DD component.
c------
        d1      = 2.0*dble(duz(j,1))+dble(wvmrs(j,1))*dble(ur(j,1))
        d1      = d1*v1*fact1
        d2      = 2.0*dble(duz(j,2))+dble(wvmrs(j,2))*dble(ur(j,2))
        d2      = d2*v2*fact2

        u1      = d1*dble(ur0(j,1))
        u2      = d2*dble(ur0(j,2))
        dkr(1)  = u1+(u2-u1)*rat

        w1      = d1*dble(uz0(j,1))
        w2      = d2*dble(uz0(j,2))
        dkz(1)  = w1+(w2-w1)*rat

        p1      = -d1*dble(tz0(j,1))
        p2      = -d2*dble(tz0(j,2))
        dkp(1)  = p1+(p2-p1)*rat

        uu1     =  d1*dble(tr0(j,1)/mu - wvmrr(j,1)*uz0(j,1))
        uu2     =  d2*dble(tr0(j,2)/mu - wvmrr(j,2)*uz0(j,2))
        ddkr(1) =  uu1+(uu2-uu1)*rat

        ww1     =  d1*
     1      dble((tz0(j,1) + wvmrr(j,1)*(l2mu - mu -mu)*ur0(j,1))/l2mu)
        ww2     =  d2*
     1      dble((tz0(j,2) + wvmrr(j,2)*(l2mu - mu -mu)*ur0(j,2))/l2mu)
        ddkz(1) =  ww1+(ww2-ww1)*rat

c------
c       SS component.
c------
        d1    = dble(wvmrs(j,1))*dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(wvmrs(j,2))*dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(2) = w1+(w2-w1)*rat
        dkr(2) = u1+(u2-u1)*rat
        dkp(2) = p1+(p2-p1)*rat

        uu1     =  d1*dble(tr0(j,1)/mu - wvmrr(j,1)*uz0(j,1))
        uu2     =  d2*dble(tr0(j,2)/mu - wvmrr(j,2)*uz0(j,2))
        ddkr(2) =  uu1+(uu2-uu1)*rat

        ww1     =  d1*
     1      dble((tz0(j,1) + wvmrr(j,1)*(l2mu - mu -mu)*ur0(j,1))/l2mu)
        ww2     =  d2*
     1      dble((tz0(j,2) + wvmrr(j,2)*(l2mu - mu -mu)*ur0(j,2))/l2mu)
        ddkz(2) =  ww1+(ww2-ww1)*rat

c------
c       DS component.
c------
        d1    = dble(wvmrs(j,1))*dble(uz(j,1)) + dble(dur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(wvmrs(j,2))*dble(uz(j,2)) + dble(dur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(3) = w1+(w2-w1)*rat
        dkr(3) = u1+(u2-u1)*rat
        dkp(3) = p1+(p2-p1)*rat

        uu1     =  d1*dble(tr0(j,1)/mu - wvmrr(j,1)*uz0(j,1))
        uu2     =  d2*dble(tr0(j,2)/mu - wvmrr(j,2)*uz0(j,2))
        ddkr(3) =  uu1+(uu2-uu1)*rat

        ww1     =  d1*
     1      dble((tz0(j,1) + wvmrr(j,1)*(l2mu - mu -mu)*ur0(j,1))/l2mu)
        ww2     =  d2*
     1      dble((tz0(j,2) + wvmrr(j,2)*(l2mu - mu -mu)*ur0(j,2))/l2mu)
        ddkz(3) =  ww1+(ww2-ww1)*rat

c------
c       EXPLOSION.
c------
        d1    = dble(duz(j,1))-dble(wvmrs(j,1))*dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(duz(j,2))-dble(wvmrs(j,2))*dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(4) = w1+(w2-w1)*rat
        dkr(4) = u1+(u2-u1)*rat
        dkp(4) = p1+(p2-p1)*rat

        uu1     =  d1*dble(tr0(j,1)/mu - wvmrr(j,1)*uz0(j,1))
        uu2     =  d2*dble(tr0(j,2)/mu - wvmrr(j,2)*uz0(j,2))
        ddkr(4) =  uu1+(uu2-uu1)*rat

        ww1     =  d1*
     1      dble((tz0(j,1) + wvmrr(j,1)*(l2mu - mu -mu)*ur0(j,1))/l2mu)
        ww2     =  d2*
     1      dble((tz0(j,2) + wvmrr(j,2)*(l2mu - mu -mu)*ur0(j,2))/l2mu)
        ddkz(4) =  ww1+(ww2-ww1)*rat

c------
c       VF component.
c------
        d1    = dble(uz(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(uz(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(5) = w1+(w2-w1)*rat
        dkr(5) = u1+(u2-u1)*rat
        dkp(5) = p1+(p2-p1)*rat

        uu1     =  d1*dble(tr0(j,1)/mu - wvmrr(j,1)*uz0(j,1))
        uu2     =  d2*dble(tr0(j,2)/mu - wvmrr(j,2)*uz0(j,2))
        ddkr(5) =  uu1+(uu2-uu1)*rat

        ww1     =  d1*
     1      dble((tz0(j,1) + wvmrr(j,1)*(l2mu - mu -mu)*ur0(j,1))/l2mu)
        ww2     =  d2*
     1      dble((tz0(j,2) + wvmrr(j,2)*(l2mu - mu -mu)*ur0(j,2))/l2mu)
        ddkz(5) =  ww1+(ww2-ww1)*rat

c------
c       HF component.
c------
        d1    = dble(ur(j,1))
        d1    = d1*v1*fact1
        u1    = d1*dble(ur0(j,1))
        w1    = d1*dble(uz0(j,1))
        p1    = -d1*dble(tz0(j,1))
        d2    = dble(ur(j,2))
        d2    = d2*v2*fact2
        u2    = d2*dble(ur0(j,2))
        w2    = d2*dble(uz0(j,2))
        p2    = -d2*dble(tz0(j,2))
        dkz(6) = w1+(w2-w1)*rat
        dkr(6) = u1+(u2-u1)*rat
        dkp(6) = p1+(p2-p1)*rat

        uu1     =  d1*dble(tr0(j,1)/mu - wvmrr(j,1)*uz0(j,1))
        uu2     =  d2*dble(tr0(j,2)/mu - wvmrr(j,2)*uz0(j,2))
        ddkr(6) =  uu1+(uu2-uu1)*rat

        ww1     =  d1*
     1      dble((tz0(j,1) + wvmrr(j,1)*(l2mu - mu -mu)*ur0(j,1))/l2mu)
        ww2     =  d2*
     1      dble((tz0(j,2) + wvmrr(j,2)*(l2mu - mu -mu)*ur0(j,2))/l2mu)
        ddkz(6) =  ww1+(ww2-ww1)*rat

c------
c       1/sqrt(2 pi).  exp(-i pi/4).  exp(-i 3pi/4).
c       A far-field approximation of Hankel function.
c------
c-----
c       introduce time shift and also phase shift from eigenfunction
c-----
        t1 =omega*tx-wvno*rx-pi4
        t2 =omega*tx-wvno*rx-pi34
        ct1=dcos(t1)
        st1=dsin(t1)
        ct2=dcos(t2)
        st2=dsin(t2)
c-----
c       DD - requires real eigenfunction
c-----
        cvz(1)     = cvz(1)     -  dkz(1)*xmom*cmplx(ct1,st1)
        cvr(1)     = cvr(1)     +  dkr(1)*xmom*cmplx(ct2,st2)
        cvp(1)     = cvp(1)     -  dkp(1)*xmom*cmplx(ct1,st1)

        cvz(1+6)   = cvz(1+6)   -  dkz(1)*xmom*cmplx(ct1,st1)*kfac
        cvr(1+6)   = cvr(1+6)   +  dkr(1)*xmom*cmplx(ct2,st2)*kfac
        cvp(1+6)   = cvp(1+6)   -  dkp(1)*xmom*cmplx(ct1,st1)*kfac

        cvz(1+6+6) = cvz(1+6+6) - ddkz(1)*xmom*cmplx(ct1,st1)
        cvr(1+6+6) = cvr(1+6+6) + ddkr(1)*xmom*cmplx(ct2,st2)
        cvp(1+6+6) = cvp(1+6+6) - ddkp(1)*xmom*cmplx(ct1,st1)
c-----
c       DS - requires -i  eigenfunction
c-----
        cvz(3)     = cvz(3)     +  dkz(3)*xmom*cmplx(st1,-ct1)
        cvr(3)     = cvr(3)     -  dkr(3)*xmom*cmplx(st2,-ct2)
        cvp(3)     = cvp(3)     +  dkp(3)*xmom*cmplx(st1,-ct1)

        cvz(3+6)   = cvz(3+6)   +  dkz(3)*xmom*cmplx(st1,-ct1)*kfac
        cvr(3+6)   = cvr(3+6)   -  dkr(3)*xmom*cmplx(st2,-ct2)*kfac
        cvp(3+6)   = cvp(3+6)   +  dkp(3)*xmom*cmplx(st1,-ct1)*kfac

        cvz(3+6+6) = cvz(3+6+6) + ddkz(3)*xmom*cmplx(st1,-ct1)
        cvr(3+6+6) = cvr(3+6+6) - ddkr(3)*xmom*cmplx(st2,-ct2)
        cvp(3+6+6) = cvp(3+6+6) + ddkp(3)*xmom*cmplx(st1,-ct1)
c-----
c       SS - requires real eigenfunction
c-----
        cvz(2)     = cvz(2)     -  dkz(2)*xmom*cmplx(ct1,st1)
        cvr(2)     = cvr(2)     +  dkr(2)*xmom*cmplx(ct2,st2)
        cvp(2)     = cvp(2)     -  dkp(2)*xmom*cmplx(ct1,st1)

        cvz(2+6)   = cvz(2+6)   -  dkz(2)*xmom*cmplx(ct1,st1)*kfac
        cvr(2+6)   = cvr(2+6)   +  dkr(2)*xmom*cmplx(ct2,st2)*kfac
        cvp(2+6)   = cvp(2+6)   -  dkp(2)*xmom*cmplx(ct1,st1)*kfac

        cvz(2+6+6) = cvz(2+6+6) - ddkz(2)*xmom*cmplx(ct1,st1)
        cvr(2+6+6) = cvr(2+6+6) + ddkr(2)*xmom*cmplx(ct2,st2)
        cvp(2+6+6) = cvp(2+6+6) - ddkp(2)*xmom*cmplx(ct1,st1)
c-----
c       EX - requires real eigenfunction
c-----
        cvz(4)     = cvz(4)     -  dkz(4)*xmom*cmplx(ct1,st1)
        cvr(4)     = cvr(4)     +  dkr(4)*xmom*cmplx(ct2,st2)
        cvp(4)     = cvp(4)     -  dkp(4)*xmom*cmplx(ct1,st1)

        cvz(4+6)   = cvz(4+6)   -  dkz(4)*xmom*cmplx(ct1,st1)*kfac
        cvr(4+6)   = cvr(4+6)   +  dkr(4)*xmom*cmplx(ct2,st2)*kfac
        cvp(4+6)   = cvp(4+6)   -  dkp(4)*xmom*cmplx(ct1,st1)*kfac

        cvz(4+6+6) = cvz(4+6+6) - ddkz(4)*xmom*cmplx(ct1,st1)
        cvr(4+6+6) = cvr(4+6+6) + ddkr(4)*xmom*cmplx(ct2,st2)
        cvp(4+6+6) = cvp(4+6+6) - ddkp(4)*xmom*cmplx(ct1,st1)
c-----
c       VF - requires real eigenfunction
c-----
        cvz(5)     = cvz(5)     -  dkz(5)*xmom*cmplx(ct1,st1)
        cvr(5)     = cvr(5)     +  dkr(5)*xmom*cmplx(ct2,st2)
        cvp(5)     = cvp(5)     -  dkp(5)*xmom*cmplx(ct1,st1)

        cvz(5+6)   = cvz(5+6)   -  dkz(5)*xmom*cmplx(ct1,st1)*kfac
        cvr(5+6)   = cvr(5+6)   +  dkr(5)*xmom*cmplx(ct2,st2)*kfac
        cvp(5+6)   = cvp(5+6)   -  dkp(5)*xmom*cmplx(ct1,st1)*kfac

        cvz(5+6+6) = cvz(5+6+6) - ddkz(5)*xmom*cmplx(ct1,st1)
        cvr(5+6+6) = cvr(5+6+6) + ddkr(5)*xmom*cmplx(ct2,st2)
        cvp(5+6+6) = cvp(5+6+6) - ddkp(5)*xmom*cmplx(ct1,st1)
c-----
c       HF - requires -i eigenfunction
c-----
        cvz(6)     = cvz(6)     +  dkz(6)*xmom*cmplx(st1,-ct1)
        cvr(6)     = cvr(6)     -  dkr(6)*xmom*cmplx(st2,-ct2)
        cvp(6)     = cvp(6)     +  dkp(6)*xmom*cmplx(st1,-ct1)

        cvz(6+6)   = cvz(6+6)   +  dkz(6)*xmom*cmplx(st1,-ct1)*kfac
        cvr(6+6)   = cvr(6+6)   -  dkr(6)*xmom*cmplx(st2,-ct2)*kfac
        cvp(6+6)   = cvp(6+6)   +  dkp(6)*xmom*cmplx(st1,-ct1)*kfac

        cvz(6+6+6) = cvz(6+6+6) + ddkz(6)*xmom*cmplx(st1,-ct1)
        cvr(6+6+6) = cvr(6+6+6) - ddkr(6)*xmom*cmplx(st2,-ct2)
        cvp(6+6+6) = cvp(6+6+6) + ddkp(6)*xmom*cmplx(st1,-ct1)
  400   continue
        if(j2.gt.j33) go to 100
c------
c       The standard output is impulse response type. i.e.,
c       ground velocity, no matter what source type is specified.
c------
410   continue
c-----
c       include source spectrum
c-----
        if(ms.eq.1)then
            sr=1.0d+00
            si=0.0d+00
        else
            sr=dble(sr0)
            si=dble(si0)
        endif
        do 500 j=1,6+6+6
            cvz(j) = cvz(j) * dcmplx(sr,si)
            cvr(j) = cvr(j) * dcmplx(sr,si)
            cvp(j) = cvp(j) * dcmplx(sr,si)
  500   continue
        do 601 j=1,6+6+6
            t1 = cdabs(cvz(j))
            if(t1.le.1.0d-20)cvz(j) = dcmplx(0.0d+00,0.0d+00)
            t2 = cdabs(cvr(j))
            if(t2.le.1.0d-20)cvr(j) = dcmplx(0.0d+00,0.0d+00)
            t2 = cdabs(cvp(j))
            if(t2.le.1.0d-20)cvp(j) = dcmplx(0.0d+00,0.0d+00)
  601   continue
        if(ipar(2).eq.0)then
c-----
c       source in solid
c-----
            xx( 1) =     cvz(1)
            xx( 2) =     cvr(1)
            xx( 3) =     cvz(3)
            xx( 4) =     cvr(3)
            xx( 6) =    -cvz(2)
            xx( 7) =    -cvr(2)
            xx( 9) =     cvz(4)
            xx(10) =     cvr(4)
            xx(11) =     cvz(5)
            xx(12) =     cvr(5)
            xx(13) =     cvz(6)
            xx(14) =     cvr(6)
            xx( 1+21) =     cvz(1+6)
            xx( 2+21) =     cvr(1+6)
            xx( 3+21) =     cvz(3+6)
            xx( 4+21) =     cvr(3+6)
            xx( 6+21) =    -cvz(2+6)
            xx( 7+21) =    -cvr(2+6)
            xx( 9+21) =     cvz(4+6)
            xx(10+21) =     cvr(4+6)
            xx(11+21) =     cvz(5+6)
            xx(12+21) =     cvr(5+6)
            xx(13+21) =     cvz(6+6)
            xx(14+21) =     cvr(6+6)
            xx( 1+42) =     cvz(1+12)
            xx( 2+42) =     cvr(1+12)
            xx( 3+42) =     cvz(3+12)
            xx( 4+42) =     cvr(3+12)
            xx( 6+42) =    -cvz(2+12)
            xx( 7+42) =    -cvr(2+12)
            xx( 9+42) =     cvz(4+12)
            xx(10+42) =     cvr(4+12)
            xx(11+42) =     cvz(5+12)
            xx(12+42) =     cvr(5+12)
            xx(13+42) =     cvz(6+12)
            xx(14+42) =     cvr(6+12)
            if(ipar(3).eq.1)then
c-----
c       receiver in fluid
c-----
                xx(16) =     cvp(4)
                xx(17) =     cvp(1)
                xx(18) =     cvp(3)
                xx(19) =    -cvp(2)
                xx(20) =     cvp(5)
                xx(21) =     cvp(6)
                xx(16+21) =     cvp(4+6)
                xx(17+21) =     cvp(1+6)
                xx(18+21) =     cvp(3+6)
                xx(19+21) =    -cvp(2+6)
                xx(20+21) =     cvp(5+6)
                xx(21+21) =     cvp(6+6)
                xx(16+42) =     cvp(4+12)
                xx(17+42) =     cvp(1+12)
                xx(18+42) =     cvp(3+12)
                xx(19+42) =    -cvp(2+12)
                xx(20+42) =     cvp(5+12)
                xx(21+42) =     cvp(6+12)
            else
c-----
c       receiver in solid
c-----
                xx(16) = 0.0
                xx(17) = 0.0
                xx(18) = 0.0
                xx(19) = 0.0
                xx(20) = 0.0
                xx(21) = 0.0
                xx(16+21) = 0.0
                xx(17+21) = 0.0
                xx(18+21) = 0.0
                xx(19+21) = 0.0
                xx(20+21) = 0.0
                xx(21+21) = 0.0
                xx(16+42) = 0.0
                xx(17+42) = 0.0
                xx(18+42) = 0.0
                xx(19+42) = 0.0
                xx(20+42) = 0.0
                xx(21+42) = 0.0
            endif
        else 
c-----
c       source in fluid
c-----
            xx(1) = 0.0
            xx(2) = 0.0
            xx(3) = 0.0
            xx(4) = 0.0
            xx(6) = 0.0
            xx(7) = 0.0
            xx(9) =     cvz(4)
            xx(10)=     cvr(4)
            xx(11)= 0.0
            xx(12)= 0.0
            xx(13)= 0.0
            xx(14)= 0.0
            xx(1+21) = 0.0
            xx(2+21) = 0.0
            xx(3+21) = 0.0
            xx(4+21) = 0.0
            xx(6+21) = 0.0
            xx(7+21) = 0.0
            xx(9+21) =     cvz(4+6)
            xx(10+21)=     cvr(4+6)
            xx(11+21)= 0.0
            xx(12+21)= 0.0
            xx(13+21)= 0.0
            xx(14+21)= 0.0
            xx(1+42) = 0.0
            xx(2+42) = 0.0
            xx(3+42) = 0.0
            xx(4+42) = 0.0
            xx(6+42) = 0.0
            xx(7+42) = 0.0
            xx(9+42) =     cvz(4+6+6)
            xx(10+42)=     cvr(4+6+6)
            xx(11+42)= 0.0
            xx(12+42)= 0.0
            xx(13+42)= 0.0
            xx(14+42)= 0.0
            if(ipar(3).eq.1)then
c-----
c       receiver in fluid
c-----
                xx(16)    = cvp(4)
                xx(17)    = 0.0
                xx(18)    = 0.0
                xx(19)    = 0.0
                xx(20)    = 0.0
                xx(21)    = 0.0
                xx(16+21) = cvp(4+6)
                xx(17+21) = 0.0
                xx(18+21) = 0.0
                xx(19+21) = 0.0
                xx(20+21) = 0.0
                xx(21+21) = 0.0
                xx(16+42) = cvp(4+6+6)
                xx(17+42) = 0.0
                xx(18+42) = 0.0
                xx(19+42) = 0.0
                xx(20+42) = 0.0
                xx(21+42) = 0.0
            else
c-----
c       receiver in solid
c-----
                xx(16)    = 0.0
                xx(17)    = 0.0
                xx(18)    = 0.0
                xx(19)    = 0.0
                xx(20)    = 0.0
                xx(21)    = 0.0
                xx(16+21) = 0.0
                xx(17+21) = 0.0
                xx(18+21) = 0.0
                xx(19+21) = 0.0
                xx(20+21) = 0.0
                xx(21+21) = 0.0
                xx(16+42) = 0.0
                xx(17+42) = 0.0
                xx(18+42) = 0.0
                xx(19+42) = 0.0
                xx(20+42) = 0.0
                xx(21+42) = 0.0
            endif
        endif
            
      return
      end
        subroutine finds(m,ii,j1,j2,j11,j22,j33,n1,n2,n3,mode)
c-----
c       In order to have no limit on the maximum number of modes
c       the modal output is written in groups of 100 modes
c       find the part where the modes>100 are stored.
c-----
        parameter (NMD=100)
        integer*4 kkr,kkl,kkf,np0
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),tt0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *      kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *      ipart1(2),ipart2(2)
        common/where/  pp(2,2)
        data epi/1.0e-8/
c-----
c
c-----
        if(ipart1(m).gt.1) then
            n1=n1+1
            krec0=0
c------
c     If want max mode < 5000, change krec0=40 to krec0=100
c       if(pp(m,2).eq.pers1(m)) krec0=40
c------
            if(pers1(m).ge.(pp(m,2)-epi)
     1          .and.pers1(m).le.(pp(m,2)+epi) )then
                    krec0=40
            endif
            irec=ii+krec0
            if(m.eq.1) then
                read(7,rec=irec) (wvml(k,1),
     1          ales(k,1),aler(k,1),gaml(k,1),ut(k,1),
     2              dut(k,1),ut0(k,1),
     3              tt0(k,1),wvmlr(k,1),wvmls(k,1),
     4              k=1,NMD)
            endif
            if(m.eq.2) then
                read(8,rec=irec) (wvmr(k,1),
     1          ares(k,1),arer(k,1),gamr(k,1),ur(k,1),
     2              dur(k,1),uz(k,1),duz(k,1),
     3              ur0(k,1),tr0(k,1),uz0(k,1),
     4              tz0(k,1),
     5              wvmrs(k,1),wvmrr(k,1),
     6                  k=1,NMD)
            endif
        endif
        if(ipart2(m).gt.1) then
            n2=n2+1
            krec0=0
c------
c     If want max mode < 5000, change krec0=40 to krec0=100
c       if(pp(m,2).eq.pers2(m)) krec0=40
c------
            if(pers2(m).ge.(pp(m,2)-epi)
     1          .and.pers2(m).le.(pp(m,2)+epi) )then
                krec0=40
            endif
            irec=ii+krec0
            if(m.eq.1) then
                read(7,rec=irec) (wvml(k,2),
     1          ales(k,2),aler(k,2),gaml(k,2),ut(k,2),
     2              dut(k,2),ut0(k,2),
     3              tt0(k,2),wvmlr(k,2),wvmls(k,2),
     4              k=1,NMD)
            endif
            if(m.eq.2) then
                read(8,rec=irec) (wvmr(k,2),
     1          ares(k,2),arer(k,2),gamr(k,2),ur(k,2),
     2              dur(k,2),uz(k,2),duz(k,2),
     3              ur0(k,2),tr0(k,2),uz0(k,2),
     4              tz0(k,2),
     5              wvmrs(k,2),wvmrr(k,2),
     6                  k=1,NMD)
            endif
        endif
        if(n1*n2.eq.0) then
            j11=j1
            j22=j2
            j33=j2
            n3=1
        else
            j11=1
c------
c     find the range of mode to be treated.
c------
            j22=NMD
            if(ipart1(m).le.ipart2(m)) then
                if(n1.eq.ipart1(m)) j22=mod(mode,NMD)
                n3=n1
            else
                if(n2.eq.ipart2(m)) j22=mod(mode,NMD)
                n3=n2
            endif
            if(j22.eq.0) j22=NMD
            j33=(n3-1)*NMD+j22
        endif
        return
        end

        subroutine frstar(r,hs,hr,mname,ipsvsh,time,pvel,svel,den,
     1      vsa, vsb, vsr, rayp, geom, tstar, dolock, dogeom)
c-----
c       r   R   Epicentral distance
c       hs  R   Source depth
c       hr  R   Receiver depth
c       mname   Ch*(*)  Name of model file
c       ipsvsh  I*4 1 - get P time
c               2 - get S  time
c               3 - get S  time
c               4 - get pP time
c               5 - get sP time
c       time    R   First arrival time
c       pvel    R   Velocity of P wave at receiver
c       svel    R   Velocity of S wave at receiver
c       den     R   Density at receiver
c       vsa R   P-wave velocity at source
c       vsb R   S-wave velocity at source
c       vsr R   Density at source
c       rayp R   Ray parameter in sec/km
c       geom R   geometrical spreading factor
c       tstar R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c       dogeom L .true. compute geometrical spreading. Only do this for
c                teleseisms
c-----
        real r, hs, hr, time, pvel, svel, vsa, vsb, vsr
        real rayp, geom, tstar
        logical dolock, dogeom
        character mname*(*)
        integer ipsvsh
c-----
c-----
c       internal variables
c-----
        real depths, depthr
        real dphs, dphr, dphref
        integer lmaxs, lmaxr, lmaxref
        real tm, rp(2),ts
        real dpdx, deginrad, sinis, sinir, cosir, cosis, sindeg
        real vs, vr, rs, rr

        
        common/earth/radius
        real radius

        common/depref/refdep
        real refdep

        radius = 6371.
c-----
c       compute the travel time
c-----
        call fstarr(r,time,lmaxs, lmaxr, lmaxref,
     1      hs, hr, ipsvsh,iflsph, rayp,
     2      tstar, dolock,
     3      mname, varec, vbrec, rhorec,
     4      vasrc, vbsrc, rhosrc)

        svel = vbrec
        pvel = varec
        den  = rhorec
        vsa  = vasrc
        vsb  = vbsrc
        vsr  = rhosrc
c-----
C       compute the geometrical spreading
C
C       since a flat earth is always used, we use the flat layered
C       geometrical spreading from Officer
C       The geometrical spreading is dimensionless and gives the decrease in amplitude from
C       a distance of 1 km from the source to the receiver
C                        2                            2
C            ( rhos vs  a sin Is Vs                  d T      )
C       sqrt |  -------------------------------     -----     |
C            |                                         2      |
C            ( rhor vr sin DEG  cos Ir Rs Cos Is     dx       )
C
C       where p (sec/km) = dT/dx
C       a = radius of sphere about source - we use 1.0 km
C       Is= incident angle at source
C       Ir= incident angle at receiver
C       Rs= distance from center of Earth to source
C       Rr= distance from center of Earth to receiver
C       DEG=epicental distance in degrees
C       rhos and vs are the density and wave velocity at the source depth
c       rhor and vr are the density and wave velocity at the receiver depth
c
c       To get the dp/dx, we determine p at different distances, and then form
c       an interpolating polynomial
c
        if(dogeom)then
              call fstarr(r-500,tm,lmaxs, lmaxr, lmaxref,
     1            hs+refdep, hr+refdep, ipsvsh,iflsph, rp(1),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)

              call fstarr(r+500,tm,lmaxs, lmaxr, lmaxref,
     1            hs+refdep, hr+refdep, ipsvsh,iflsph, rp(2),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)


              dpdx = abs(rp(1) - rp(2))/(500.0 - ( - 500.0))
              deginrad = r/radius
              sindeg = sin(deginrad)
  
              if(ipsvsh.eq.1)then
                  vs = vsa
                  vr = pvel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.2)then
                  vs = vsb
                  vr = svel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.3)then
                  vs = vsb
                  vr = svel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.4)then
                  vs = vsa
                  vr = pvel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.5)then
                  vs = vsb
                  vr = pvel
                  rs = vsr
                  rr = den
              endif
              sinis = rayp*vs
              cosis = sqrt(1.0-sinis*sinis)
              sinir = rayp*vr
              cosir = sqrt(1.0-sinir*sinir)

              fac = (rs*vs*sinis*vs*dpdx)/
     1           (rr*vr*sindeg*cosir*(radius-hs)*cosis)
              geom = sqrt(abs(fac))
      else
c-----
c     default to have a defined return value
c-----
              geom = 1.0
      endif
              

c-----
        return
        end

        subroutine fstarr(dist,tfirst,lmaxs,lmaxr,lmaxref,
     1      hs,hr,ipsvsh,iflsph, rayp,
     2      tstar, dolock,
     3      mname, varec, vbrec, rhorec,
     4      vasrc, vbsrc, rhosrc)
c-----
c       given a distance, the source depth, receiver depth,
c       get time of first arrival of P
c-----
c       dist    R   - distance
c       tfirst  R   - first arrival time
c       mmax    I*4 - number of layers in model
c       lmaxs   I*4 - layer index for source
c       lmaxr   I*4 - layer index for receiver
c       lmaxref I*4 - layer index for reference depth,
c                     used only for pP and sS
c       hs      R   - depth of source
c       hs      R   - depth of receiver
c       ipsvsh  I*4 1 - get P time
c               2 - get S  time
c               3 - get S  time
c               4 - get pP time
c               5 - get sP time
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       rayp    R   - ray parameter in sec/km
c       geom R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c       mname   Ch* - name of the velocity model
c       varec   R - P velocity at receiver (untransformed)
c       vbrec   R - S velocity at receiver (untransformed)
c       rhorec  R - Density at receiver (untransformed)
c       vasrc   R - P velocity at source (untransformed)
c       vbsrc   R - S velocity at source (untransformed)
c       rhosrc  R - Density at source (untransformed)
c-----
c       since this routine is to be used for omega-k,
c       we will approximate the direct arrival
c
c       18 JAN 2008 - everything is straightforward. The addition of
c          the request for pP and sP changes the logic in that
c          the direct arrival is ignored, and that the upgoing refraction 
c          from the source is ignored. We handle this by just setting
c          a very large tfirst before trying to do the modified 
c          downward path refraction to avoid another level of
c          if/then/else/endif
c       24 MAR 2008 - modified to place the model read into this
c          routine instead of in frstar
c-----
        real dist, tfirst, dphs, dphr, hs, hr, depths, depthr
        real rayp
        real varec, vbrec, rhorec
        real vasrc, vbsrc, rhosrc
        logical dolock
        integer lmaxs, lmaxr, lmaxref
        character mname*(*)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmmax
        integer mmmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer mmax

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        real v(NL), h(NL), qi(NL)

        real*8 c, s, t, x, p, tint, dxdp, vel, pnew, pupper
        real*8 ts
        real*8 sumx, sumt
        logical ext

        real tds, tdr
        common/earth/radius
        real radius

c-----
c       first read in the model and determine the medium parameters at the
c       source and receiver depths
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        l = lgstr(mname)
        if(.not. ext) then
             write(0,*)'Model:', mname(1:l)
             call usage('Model file does not exist')
        endif

                call getmod(11,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        mmmax = mmax
        if(ierr .lt. 0)return
                call adomod()
c-----
c       insert the source and receiver depths into the model
c       placing the source and receiver on a layer boundary
c-----
        call insert(hs+refdep)
        call insert(hr+refdep)
        call insert(   refdep)

c-----
c       get the layer in which the source lies
c-----
        call srclyr(hs+refdep, lmaxs, dphs)
        call srclyr(hr+refdep, lmaxr, dphr)
        call srclyr(   refdep, lmaxref, dphref)
c-----
c       get the medium parameters at the source and reciever depths
c-----
        varec = a(lmaxr)
        vbrec = b(lmaxr)
        rhorec = rho(lmaxr)
        vasrc = a(lmaxs)
        vbsrc = b(lmaxs)
        rhosrc = rho(lmaxs)
c-----
c       prepare for the computations
c-----
        depths = hs + refdep
        depthr = hr + refdep

c-----
c       set up default
c-----
        tfirst = 1.0e+30
c-----
c       special case for locked mode
c-----
        if(dolock)then
            mmax = mmmax -1
        else
            mmax = mmmax
        endif

c-----
c       get specifics about upward and downward distances
c       with a layer. We need his to define ray paths
c       We will also use the fact that the source/receiver are
c       on layer boundaries
c
c       lmn = layer number of shallowest of source/receiver
c       lmx = layer number of deepest    of source/receiver
c-----
        lmn = min(lmaxs,lmaxr)
        lmx = max(lmaxs,lmaxr)

c-----
c       perform spherical -> flat earth model transformation
c-----
        if(iflsph.ne.0)then
            call adosph()
            tds = radius*alog(radius/(radius-hs))
            tdr = radius*alog(radius/(radius-hr))
        else
            tds = depths
            tdr = depthr
        endif
c-----
c       now fill in velocity array according to desired first arrival
c       for S  there can be no water layer
c       for S  can be a water layer
c       Also define the Q for the T* analysis. Note we define
c        eventually q = 1/Q based on whether the given Q > or < 1
c-----
        do 100 i=1,mmax
            if(ipsvsh.eq.1)then
                v(i) = a(i)
                qi(i) = qa(i)
            else if(ipsvsh.eq.2)then
                v(i) = b(i)
                qi(i) = qb(i)
                if(b(i).le.0.001)then
                    v(i) = a(i)
                    qi(i) = qa(i)
                endif
            else if(ipsvsh.eq.3)then
                v(i) = bsh(i)
                qi(i) = qbsh(i)
            else if(ipsvsh.eq.4)then
                v(i) = a(i)
                qi(i) = qa(i)
            else if(ipsvsh.eq.5)then
                v(i) = a(i)
                qi(i) = qa(i)
            endif
            if(qi(i) .gt. 1.0)then
                qi(i) = 1.0 / qi(i)
            endif
            h(i) = d(i)
 100    continue
c-----
c       For the computations we look at four cases
c       1) direct path between source and receiver 
c          a) source and receiver in the same layer
c          b) source and receiver in different layers
c       2) refracted arrivals       
c          a) path is downward from source and then up to
c             receiver
c          b) path is upward from the source and then down to
c             receiver
c          This recognized the possibility that velocity does
c          not increase uniformly with depth
c-----
                    
c-----
c       direct arrival source/receiver at same layer
c-----
        if(v(lmaxs).eq.0.0)return
        if(v(lmaxr).eq.0.0)return
        if(lmaxs .eq. lmaxr)then
            tfirst = sqrt(dist**2 + abs(tds - tdr)**2)/
     1          v(lmaxs)
            rayp = (dist/sqrt(dist**2 + abs(tds - tdr)**2))/
     1          v(lmaxs)
            tstar = tfirst*qi(lmaxs)
        else
c-----
c       direct arrival source/receiver in different layers
c-----
c       Newton Iteration for direct arrival source/receiver at
c           different depths
c           
c           x = SUM h  tan theta
c                    i          i
c
c           t = SUM h  / v  cos theta
c                    i    i          i
c                                                          2 2
c       where sin theta  = p V  , cos theta  = sqrt ( 1 - p V )
c                      i      i                              i
c       and p is the ray parameter bounded by [0, 1/V  ] where V
c                                                    sr         sr
c       is the wave velocity at the starting point of the ray. 
c       Since the ray must also reach the receiver, we consider
c       that possibility too. The upper bound is MIN ( 1/Vs, 1/Vr)
c       Also we test for a real ray path, between source and receiver
c
c       Because source/receiver at top of a layer boundary, we have
c
c           -----------X----------
c           h(lmn)      \
c           ----------------------
c                      ....
c           ----------------------
c           h(lmx-1)        \
c                            \
c           ------------------X---
c            
c-----
            ps = 1.0/v(lmaxs)
            pr = 1.0/v(lmaxr)
            if(ps.lt.pr)then
                pupper = ps
            else
                pupper = pr
            endif
            do 1000 l=lmn,lmx
                if(v(l).eq.0.0)return
                p = 1.0/v(l)
                if(p.lt.pupper)pupper = p
 1000       continue
            p = 0.5d+00  * pupper
            do 1111 iter=1,10
                x = 0.0d+00
                t = 0.0d+00
                ts = 0.0d+00
                tint = 0.0d+00
                dxdp = 0.0d+00
                do 1500 l=lmn,lmx - 1
                    vel = dble(v(l))
                    s = p*vel
                    c = dsqrt(1.0d+00 - s*s)
                    t = t + dble(h(l)) /(vel*c)
                    x = x + dble(h(l)) * s / c
                    dxdp  = dxdp + dble(h(l)) *
     1                  vel/(c*c*c)
                    tint = tint + dble(h(l)) * c / vel
                    ts = ts + qi(l) * dble(h(l))/(c*vel)
                   

 1500           continue
                pnew = p - (x-dble(dist))/dxdp
c-----
c       safety - we must have a real ray, with upper bound
c       of  min[ 1/v(src), 1/v(rec)]
c-----  
                if(pnew .gt. pupper)then
                    if(iter.lt.10)then
                        pnew = 0.99d+00 * pupper
                    else
c-----
c       this is propably working like a refraction, so stop iterations
c-----  
                        t = tint + p * (dist)
                        go to 1112
                    endif
                endif
                p = pnew
 1111       continue
 1112       continue
            tfirst = t
            rayp = p
            tstar = ts
        endif
c-----
c       now proceed through the possible refracted arrivals
c       considering first upward rays from the source
c-----  
        if(lmn.gt.1)then
        do 3020 m=1,lmn-1
c-----
c       m is the refracting layer
c
c       get velocity of refractor
c-----
            vel = v(m)
            if(v(m).eq.0.0)return
            p = 1.0/vel
c-----
c
c           --------------------------------
c           h(1)
c           --------------------------------
c                      ....
c           --------------------------------
c           h(m)
c           ----------------...-------------
c           h(m+1)         /   \
c           --------------------------------
c                         /     \
c                      ....
c           --------------------------------
c           h(lmn-1)              \
c           -----------------------X--------
c               
c           h(lmn)     /    
c           --------------------------------
c                      ....
c           --------------------------------
c           h(lmx-1) /
c           --------X-----------------------
c
c       safety check, velocity at source or receiver must be less than
c       refraction velocity
c-----
        if(v(lmn).ge.vel)go to 3020
        if(v(lmx).ge.vel)go to 3020
c-----
c       single leg
c-----
        sumt = 0.0
        sumx = 0.0
        ts = 0.0
            do 3021 l=1,lmx-1,lmn
                if(v(l).gt.vel)go to 3020
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + h(l)*cs/v(l)
                sumx = sumx + h(l)*p/cs
                ts = ts + qi(l)*h(l)/(cs * v(l))
 3021       continue
            do 3022 l=m+1,lmn-1
                if(v(l).gt.vel)go to 3020
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + 2.0*h(l)*cs/v(l)
                sumx = sumx + 2.0*h(l)*p/cs
                ts = ts + 2.0*qi(l)*h(l)/(cs * v(l))
 3022       continue
            tint = sumt
            tt = tint + dist / vel
            ts = ts + qi(m)*(dist-sumx)/v(m)
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                  tfirst = tt
                  rayp = p
                 tstar = ts
            endif
 3020       continue
        endif
c-----
c       For the special case of the depth phases, ignore previous
c       first arrival times
c-----
        if(ipsvsh.eq.4 .or. ipsvsh.eq.5)then
             tfirst = 1.0e+30
        endif
c-----
c       now proceed through the possible refracted arrivals
c       considering first downward rays from the source
c
c       We start considering the deepest point since we place
c       a source/receiver position just below a layer boundary
c       and thus should consider a horizontal ray
c
c       The refraction is accepted only if the desired distance >
c       first refraction from the source - this puts physics in the problem
c           
c           x = SUM h  tan theta
c                    i          i
c
c           t = SUM h  cos theta / V
c                    i          i   i
c                                                          2 2
c       where sin theta  = p V  , cos theta  = sqrt ( 1 - p V )
c                      i      i                              i
c       For the T* computation we need to follow the path, e.g.,
c       SUM h qi / ( cos theta  / V ) + qi (dist -  SUM h tan theta / V )/V
c            i  i             i    i      i              i         i   i   r
c-----  
        do 2020 m=lmx+1, mmax
c-----
c       m is the refracting layer
c
c       get velocity of refractor
c-----
            vel = v(m)
            if(v(m).eq.0.0)return
            p = 1.0/vel
c-----
c
c           -----------X--------------------
c           h(lmn)      \
c           --------------------------------
c                      ....
c           --------------------------------
c           h(lmx-1)        \             
c                            \           
c           ------------------X--------X----
c           h(lmx)             \       /
c           --------------------\-----/-----
c                      ....      \   /
c           ----------------------...-------
c           h(m)
c
c-----
c       safety check, velocity at source or receiver must be less than
c       refraction velocity
c-----
        if(v(lmn).ge.vel)go to 2020
        if(v(lmx).ge.vel)go to 2020
c-----
c       single leg
c-----
        sumx = 0.0
        sumt = 0.0
        ts = 0.0
c-----
c       special case for depth phases
c-----
            if(ipsvsh.eq.4)then
c-----
c               pP
c-----
                  do  l=lmaxref,lmaxs - 1
                      if(a(l).gt.vel)go to 2020
                      cs = sqrt(abs(1.0 - p*p*a(l)*a(l)))
                      sumt = sumt + 2.*h(l)*cs/a(l)
                      sumx = sumx + 2.*h(l)*p*a(l)/cs
                      if(qa(l).gt.1.0)qa(l) = 1.0/qa(l)
                      ts = ts + 2.*qa(l)*h(l)/(cs * a(l))
                  enddo
            else if(ipsvsh.eq.5)then
c-----
c               sP
c-----
                  do  l=lmaxref,lmaxs - 1
                      if(a(l).gt.vel)go to 2020
                      if(b(l).gt.vel)go to 2020
                      csa = sqrt(abs(1.0 - p*p*a(l)*a(l)))
                      csb = sqrt(abs(1.0 - p*p*b(l)*b(l)))
                      sumt = sumt + h(l)*csa/a(l)
     1                        +h(l)*csb/b(l)
                      sumx = sumx + 2.*h(l)*p*a(l)/csa
                      if(qa(l).gt.1.0)qa(l) = 1.0/qa(l)
                      if(qb(l).gt.1.0)qb(l) = 1.0/qb(l)
                      ts = ts + qa(l)*h(l)/(csa * a(l))
     1                        + qb(l)*h(l)/(csb * b(l))
                  enddo
            endif
c-----
c       continue
c-----
            do 2021 l=lmn,lmx - 1
                if(v(l).gt.vel)go to 2020
                if(v(l).eq.0.0)return
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + h(l)*cs/v(l)
                sumx = sumx + h(l)*p*v(l)/cs
                ts = ts + qi(l)*h(l)/(cs * v(l))
 2021       continue
c-----
c       double leg
c-----
            do 2022 l=lmx,m-1
                if(v(l).gt.vel)go to 2020
                if(v(l).eq.0.0)return
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + 2.0*h(l)*cs/v(l)
                sumx = sumx + 2.0*h(l)*p*v(l)/cs
                ts = ts + 2.*qi(l)*h(l)/(cs * v(l))
 2022       continue
            tint = sumt
            tt = tint + dist / vel
            ts = ts + qi(m)*(dist-sumx)/vel
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                 tfirst = tt
                 rayp = p
                 tstar = ts
            endif
 2020       continue
             if(tfirst .eq. 1.0e+30)then
                tfirst = -12345.
                rayp   = -12345.
                tstar  = -12345.
             endif
        return
        end

        subroutine gtemod(depths,depthr,nper,mname,ipar,fpar)
c------
c       read in the earth model. consistency check between files.
c-----
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NLAY=200)
        real*4 d(NLAY),a(NLAY),b(NLAY),rho(NLAY),qa1(NLAY),qb1(NLAY)
        integer*4 kkr,kkl,kkf,mmax,nperl,mperl,nperr,mperr,nper,np0
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
        common/water/rrho
c-----
        dphsr=0.0
        dphre=0.0
        mperr = 0
        mperl = 0
        ijk=1
        if(kkl.eq.1) then
            rewind 1
            call getmdl(1,mmax,d,a,b,rho,qa1,qb1,nperl,dphsrl,dphrel,
     2              mname,ipar,fpar)
            depths=dphsrl
            depthr=dphrel
            nper=nperl
        endif
        if(kkr.eq.1) then
            rewind 2
            call getmdl(2,mmax,d,a,b,rho,qa1,qb1,nperr,dphsrr,dphrer,
     2              mname,ipar,fpar)
            depths=dphsrr
            depthr=dphrer
            nper=nperr
            deplw = 0.0
            depup = 0.0
            do 1000 i=1,mmax
                depup = deplw + d(i)
                if(depthr.ge.deplw .and. depthr.lt.depup)then
                    rrho = rho(i)
                endif
                deplw = depup
 1000       continue
        endif
c       write(LER,8)
c       write(LER,10)
C       write(LER,20) (i,d(i),a(i),b(i),rho(i),qa1(i),qb1(i),i=1,mmax)
        bmax=b(mmax)
c       write(LER,30) dphsr,dphre,nperl,mperl,nperr,mperr
        return
        end

        subroutine insert(dph)
        implicit none
        real dph
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        real dp, dphh, hsave, dep
        integer m, ls
c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        if(dph.le.0)then
            d(1) = d(1) - dph
            return
        else if(dph.ge.0)then
            dep = 0.0 
            dp = 0.0 
            dphh = -1.0
            do 100 m = 1,mmax 
                dp = dp + d(m) 
                dphh = dp - dph 
                if(m.eq.mmax)then
                    if(d(mmax).le.0.0 .or. dphh.lt.0.0)then
                        d(mmax) = (dph - dp)
                    endif
                endif
                dep = dep + d(m) 
                dphh = dep - dph 
                ls = m 
                if(dphh.ge.0.0) go to 101 
  100       continue 
  101       continue 
        endif
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c           at surface and internally
c       However do put in a zero thickness layer 
c           at the base if necessary
c-----
        if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            return
        else
c-----
c           adjust layering
c-----
             do 102 m = mmax,ls,-1
                d(m+1) = d(m)
                a(m+1) = a(m)
                b(m+1) = b(m)
                rho(m+1) = rho(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                etap(m+1) = etap(m)
                etas(m+1) = etas(m)
                frefp(m+1) = frefp(m)
                frefs(m+1) = frefs(m)
                bsh(m+1) = b(m)
                qbsh(m+1) = qb(m)
                rhosh(m+1) = rho(m)
  102       continue
            hsave = d(ls)
            d(ls) = hsave - dphh
            d(ls+1) = dphh
            ls = ls + 1
            mmax = mmax + 1
            if(d(mmax).lt.0.0)d(mmax)=0.0
        endif
        return
        end

        subroutine integ(x,nt,dt,dodc)
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       Form the integral x(t) dt numerically
c       x   R*4 - array to be integrated
c       nt  I*4 - number of data points in time series
c       dt  R*4 - sampling interval
c       dodc    L*4 - .true. remve x(1) before summing
c-----
        
        real x(nt)
        integer*4 nt
        real*4 dt
        logical dodc
        sum=0.0
        if(dodc)then
            do 100 k = 1,nt
                sum = sum + dt*(x(k) - x(1))
                x(k)=sum
  100       continue
        else
            do 200 k = 1,nt
                sum = sum + dt*x(k)
                x(k)=sum
  200       continue
        endif
        return
        end

        subroutine maksyn(ntau,ipt,alp,idva,
     1      xmult,rfile,ryfile,lvfile,dist,dt,npts,
     2      tfirst,ostr,
     2      mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat,
     4      dotest1,dostep,dostrain,dostress)
        
c-----
c                                                                       
c     This program generates sixteen basic types of synthetic          
c     spectrum after combining the eigenfunctions from                 
c     sregn96 and/or slegn96  with the source spectrum.              
c                                                                    
c     The output is a time series in file96 format
c                                                                 
c     When number of points <= NSAMP, do FFT (call sub four) here,
c     otherwise, it will need 'supfft96' (through system call). 
c                                                              
c     The maximum dimension for time histories is LIMITED  by dimension
c     The maximum mode number at any period is NOT LIMITED.  
c                                                           
c     The Q is independent of frequency                   
c                                                        
c     The program uses a large number of files in order to function 
c     The file Logical Unit Numbers and their use are as follow:   
c                                                                 
c     1       Love Eigenfunction Source Depth                    
c     2       Rayleigh Eigenfunction Source Depth               
c     3       Distance control file                            
c     4       Temporary source pulse file used in pulud       
c     5       Standard Input                                 
c     6       Standard Output                               
c     7       Temporary Love Mode File                     
c     8       Temporary Rayleigh Mode File                
c     9       Internal Spectra File and then Output time series File 
c                                                                   
c     The program may seem complicated, but the computational structure
c     is designed for efficiency and relies upon the order of the output
c      of sdisp96, sregn96 and slegn96
c                                                                     
c     For each frequency, the 16 Green s functions are computed
c     at each distance.
c-----
        parameter (LER=0, LIN=5, LOT=6)
        character*80 ryfile, lvfile, rfile
        integer*4 ntau, ipt, idva,  npts
        real*4 xmult, alp, dist, dt, tfirst
        character ostr*(*)
        logical dodble, dolock, dozero, dolat, dotest1,dostep
        logical dostrain, dostress
C-----
c       MSDOS, and SUN F77
c-----
        integer kerr, system
C-----
c       END MSDOS
c-----

        integer NSAMP2
        parameter(NSAMP2=16384)

        logical ext
        parameter (NGRN=63)
        integer*4 kkr,kkl,kkf,nper
        integer*4 np0
        integer*4 isign
        integer NSAMP
        parameter (NSAMP=16384)
          common/srctim/ src(NSAMP)
        parameter (NMD=100)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),tt0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     1                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     1                 ipart1(2),ipart2(2)

        common/receiver/l2mu,mu
        real l2mu, mu
c-----
c       at the receiver position
c       l2mu = rho Vp^2, mu = rho Vs^2
c-----

        character mname*80
        integer ipar(10)
        real*4 fpar(10)

        integer iopen0
        save iopen0
        data iopen0/0/
        itrigl = 0
        itrigr = 0
c-----
c       check for the existence of the Rayleigh and Love wave
c       eigenfunction files
c-----
        inquire(file=lvfile,exist=ext)
        if(ext .and. lvfile.ne.' ')then
            kkl=1
            open(1,file=lvfile,status='old',form='unformatted',
     1          access='sequential')
        else
            kkl=0
        endif
        inquire(file=ryfile,exist=ext)
        if(ext .and. ryfile.ne.' ')then
            kkr=1
            open(2,file=ryfile,status='old',form='unformatted',
     1          access='sequential')
        else
            kkr=0
        endif
        
c-----
c       enter source and other parameters:
c-----
c       Mode selection
c       kkf =   1   all modes
c           2   fundamental only
c           3   higher modes only
c           4   range
c                   mods = first
c                   modt = last
c                   modi = increment
c-----
        if(nmode.gt.0)then
            kkf = 4
            mods = nmode
            modt = nmode
            modi = 1
        else if(nmode.eq.-1)then
            kkf = 4
            mods = 1
            modt = 10000
            modi = 1
        else if(nmode.eq.-2)then
            kkf = 4
            mods = 1
            modt = 1
            modi = 1
        else if(nmode.eq.-3)then
            kkf = 4
            mods = 2
            modt = 10000
            modi = 1
        else
            kkf = 4
            mods = 1
            modt = 10000
            modi = 1
        endif
        xmom = xmult
c-----
c       define the source spectrum
c-----
        call srcpul(ntau,ipt,alp,dt,npts,rfile,duration)
        if(.not. dozero)then
            duration = 0.0
        endif
        zeta = 0.0
c-----
c       check for power of 2
c-----
        call npow2(npts)
        if(dodble)npts = npts * 2
        if(npts.gt.NSAMP)npts=NSAMP
c-----
c       read the eigenfunction headers which contain
c       informaiton on source and reiver depths
c-----
        call gtemod(depths,depthr,nper,mname,ipar,fpar)
c-----
c       get elocityes and densities at source and receiver depths
c       and then compute lambda + 2mu and mu at the receiver depth
c-----
            call frstar(dist,depths,depthr,mname,1,tp ,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            mu  = den*svel*svel
            l2mu = den*pvel*pvel
c       here den = GM/CM^3 VEL=KM/S 
c-----
c-----
c       process the spectra for this distance
c-----
        np0=npts
c-----
c       output the file16 header for this distance
c-----
        df=1./(float(np0)*dt)
c------
c       generate the source spectrum.
c       Due to limitation of available space, when np0>NSAMP
c       use 'supfft' to do FFT.
c------
            do 330 i=1,NSAMP
                ssrc(i) = cmplx(0.0,0.0)
  330       continue
            do 340 i=1,NSAMP
                ssrc(i)=cmplx(src(i),0.0)
  340       continue
            call zfour(ssrc,np0,-1,dt,df)
c------
c       open temporary output file for spectra
c-----
        open(10,file='keep.pl0',status='unknown',form='unformatted')
        iopen0=1
        rewind 10

c------
c       the main call to process the spectra for this distance.
c------
        call gtegn(dist,tfirst,ipar)
c------
c       output.
c------
        if(itrigl.eq.1) close(7,status='delete')
        if(itrigr.eq.1) close(8,status='delete')
        if(kkl.eq.1) close(1)
        if(kkr.eq.1) close(2)
        if(dotest1)then
             call outtest1(ntau,idva,dist,dt,npts,
     1           tfirst-0.5*duration,depthr,depths,kkl,kkr,ostr,
     2           mname,ipar,fpar,dodble,dolock,dolat,dostep)
        else
             call output(ntau,idva,dist,dt,npts,
     1           tfirst-0.5*duration,depthr,depths,kkl,kkr,ostr,
     2           mname,ipar,fpar,dodble,dolock,dolat,dostep,dostrain)
             if(dostress)then
             call outstress(ntau,idva,dist,dt,npts,
     1           tfirst-0.5*duration,depthr,depths,kkl,kkr,ostr,
     2           mname,ipar,fpar,dodble,dolock,dolat,dostep)
             endif
        endif
        if(iopen0.eq.1) close(10,status='delete')
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsmp, npts
        nsmp = npts
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsmp)goto 1000
        return
        end

        subroutine process(ntau,ipt,alp,dfile,idva,
     1      xmult,rfile,ryfile,lvfile,ostr,
     2      mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat,
     3      dotest1,dostep,dostrain,dostress)
c-----
c       make Green s functions for a given source pulse and
c       distance
c-----
        parameter (LER=0, LIN=5, LOT=6)
        character*80 dfile, ryfile, lvfile, rfile
        integer*4 ntau, ipt, idva,   nmode
        real*4 xmult, alp
        logical dolock, dozero,dolat, dotest1,dostep
        logical dostrain, dostress

        real*4 dist, dt, tfirst
        integer npts
        logical ext
        character ostr*(*)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
        logical dodble
c-----
c       open file containing distance information
c-----
        inquire(file=dfile, exist=ext)
        if(.not. ext)then
            call usage('Distance Information File does not exist')
        endif
        open(3,file=dfile, status='old',form='formatted',
     1      access='sequential')
        rewind 3
c-----
c       read the data from the distance file, and process each distance
c       separately
c-----
 1000   continue
            read(3,*,err=9999,end=9999)dist,dt,npts,t0,vred
c-----
c       dist    R*4 - desired distance from the source
c               cannot be zero
c       dt  R*4 - desired sampling interval
c       npts    R*4 - desired number of points in the time series
c       t0,vred R*4 - control for the time of the first sample
c               if vred > 0, then the
c                   tfirst = t0 + dist/vred
c               else
c                   tfirst = t0
c-----
            dist = abs(dist)
            if(vred.gt.0.0)then
                tfirst = t0 + dist/vred
            else
                tfirst = t0
            endif
            call maksyn(ntau,ipt,alp,idva,
     1          xmult,rfile,ryfile,lvfile,dist,dt,npts,
     2          tfirst,ostr,
     3          mname,ipar,fpar,dodble,nmode,dolock,dozero,dolat,
     4          dotest1,dostep,dostrain,dostress)
        go to 1000
 9999   continue
        close(3)
        return
        end

        subroutine puldd(dt,n,l)
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/d(NSAMP)
c-----
c       Dirac Delta Pulse
c-----
        do 100 i=1,n
            d(i) = 0.0
  100   continue
        d(1) = 1.0/dt
        l = 2
        return
        end

        subroutine pulod(alp,dt,nt,l)
c-----
c       unit area far field displacement Ohnaka pulse
c           Harkrider (1976) Geophys J. 47, p 97.
c-----
c       alp R*4 - shape parameter, corner frequency
c               fc = alp / 2 pi
c       dt  R*4 - sample rate
c       nt  I*4 - number of points in time series
c       l   I*4 - pulse d(i) = 0 for i >=l
c-----
        real*4 alp, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/d(NSAMP)
        ltest = 0
        al2=alp*alp
        do 100 i=1,nt
            t=(i-1)*dt
            d(i)=0.0
            arg= alp*t
            if(arg.le.25.0)then
                d(i)= al2*t*exp(-arg)
            else
                ltest = ltest +1
                if(ltest.eq.1)l = i
            endif
  100   continue
        return
        end

        subroutine pulpd(tau,dt,nt,l)
c-----
c       unit area far field displacement parabolic pulse
c-----
c       tau R*4 - duration parameter, total duration = 4 tau
c       dt  R*4 - sample rate
c       nt  I*4 - number of points for time series
c       l   I*4 - pulse d(i) = 0 for i >=l
c-----
        real*4 tau, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=16384)
        common /srctim/ d(NSAMP)
        ltest = 0
        tl = tau
        t1 = 0.01*dt
        t2 = t1 + tau
        t3 = t2 + tau
        t4 = t3 + tau
        t5 = t4 + tau
        do 100 i = 1,nt
            y=(i-1)*dt
            z = y - t1
            d(i) = 0.0
            if(y.ge.t1 .and. y.lt.t2)then
                d(i) = 0.5*(z/tl)**2
            else if(y.ge.t2 .and. y.le.t3)then
                d(i)= -0.5*(z/tl)**2 + 2.*(z/tl) - 1
            else if(y.ge.t3 .and. y.le.t4)then
                d(i)= -0.5*(z/tl)**2 + 2.*(z/tl) - 1.
            else if(y.ge.t4 .and. y.le.t5)then
                d(i)= 0.5*(z/tl)**2 - 4.*(z/tl) + 8.
            else
                ltest = ltest + 1
                if(ltest.eq.1) l = i
            endif
  100   continue
c-----
c       pulse normalized so first integral has area of unity
c-----
        do 200 i = 1,nt
            d(i) = d(i)/(2.*tl)
  200   continue
        return
        end

        subroutine pultd(tau,dt,nt,l)
c-----
c       unit area far field displacement triangular pulse
c-----
c       tau R*4 - duration parameter, total duration = 2 tau
c       dt  R*4 - sample rate
c       nt  I*4 - number of points for time series
c       l   I*4 - pulse d(i) = 0 for i >=l
c-----
        real*4 tau, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=16384)
        common /srctim/ d(NSAMP)
        ltest = 0
        fac = 1./tau
        t1 = tau
        t2 = tau + tau
        do 100 i=1,nt
            t = (i-1)*dt
            d(i)=0.0
            if(t.le.t1)then
                z = t - 0.0
                d(i) = z*fac
            elseif(t.gt.t1.and.t.le.t2)then
                z = t - t1
                d(i)= fac - z*fac
            elseif(t.gt.t2)then
                d(i)=0.0
                ltest = ltest + 1
                if(ltest.eq.1)l = i
            endif
  100   continue
        return
        end

        subroutine pulud(rfile,n,tau,dtt)
        character rfile*(*)
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/d(NSAMP)
        do 100 i=1,n
            d(i) = 0.0
  100   continue
        open(4,file=rfile,status='unknown',form='formatted',
     1      access='sequential')
        rewind 4
        read(4,*)np,dtt
        read(4,*)(d(i),i=1,np)
        close(4)
        tau = np*dtt
        return
        end

        subroutine readin(igo,jj,ifunc,mode,ipart,per)
c-----
c       read in the eigenfunction files generated by
c       reigen85 and leigen85
c
c       use the definitions of the eigenfunction and stress to get
c       the necessary vertical derivatives
c
c       Instead of reading only a single mode at a time, groups of
c       NMD are read in - the purpose of this is to 
c            significantly reduce
c       IO time. If not for this the code would be much simpler
c-----
c       igo I*4 - 0
c                 1
c                 2
c       jj  I*4 -
c       ifunc   I*4 - ifunc for L and R
c       mode    I*4 -
c       ipart   I*4 -
c       per R*4 -
c-----
        parameter (NMD=100)
        integer*4 kkr,kkl,kkf,ifunc(2),mode,np0
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),tt0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *                 kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *                 ipart1(2),ipart2(2)
        common/where/  pp(2,2)
        parameter (NLAY=200)
        real*4 d(NLAY),a(NLAY),b(NLAY),rho(NLAY),qa1(NLAY),qb1(NLAY)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       set fpar(1) = refdep
c       set ipar(1) = 1 if medium is spherical
c       set ipar(2) = 1 if source is in fluid
c       set ipar(3) = 1 if receiver is in fluid
c       set ipar(4) = 1 if eigenfunctions are output with -DER flag
c       set ipar(5) = 1 if dc/dh are output with -DER flag
c       set ipar(6) = 1 if dc/da are output with -DER flag
c       set ipar(7) = 1 if dc/db are output with -DER flag
c       set ipar(8) = 1 if dc/dr are output with -DER flag
c-----
c       initialize
c-----
        if(igo.eq.0)then
            do 90 k=1,2
                ipart1(k)=0
                ipart2(k)=0
   90       continue
c------
c       skip the first several records, which contain the model.
c       This is required only for proper input positioning
c------
            if(kkl.eq.1) then
                rewind 1
            call getmdl(1,mmax,d,a,b,rho,qa1,qb1,nperl,dphsrl,dphrel,
     2              mname,ipar,fpar)
            endif
            if(kkr.eq.1) then
                rewind 2
            call getmdl(2,mmax,d,a,b,rho,qa1,qb1,nperr,dphsrl,dphrel,
     2              mname,ipar,fpar)
            endif
        else
c-----
c       start reading the eigenfunction file 
c       The structure of the file, per read, is
c       earth model
c       number_periods, source_depth
c            | ifunc,mode,per [ifunc=1 for Lv, 
c            2 for RY, -1 for end of file]
c           p|m|        [ mode - number of modes ]
c           e|o|        [ per period of this entry, < 0 
c            end of per set ]
c           r|d| eigenfunctions for this mode and period
c-----
            icode=igo
            call gethed(icode,ifunc(icode),mode,per,ierr)
            if(ifunc(icode).le.0) return
c------
c     Use a buffer with NMD data length to handle large mode number.
c     IF f77 allows bigger dimension, the change in the parameter NMD
c     number will make this program faster.
c
c       Current Love IO has  7 variables, need recl= 7*4*NMD
c       Current Rayl IO has 11 variables, need recl=11*4*NMD
c------
            ipart=(mode-1)/NMD+1
            if(mode.le.0) return
            if(ipart.gt.1) then
                if(itrigl.eq.0.and.jj.eq.1) then
                    locl=0
                    itrigl=1
                    open(7,file='keep.pl7',
     1                  status='unknown',
     2                  form='unformatted',recl=4096
     3                  ,access='direct')
                endif
                if(itrigr.eq.0.and.jj.eq.2) then
                    locr=0
                    itrigr=1
                    open(8,file='keep.pl8',
     1                  status='unknown',
     2                  form='unformatted',recl=6144
     3                  ,access='direct')
                endif
            endif
c------
c     use double buffer. This is a difficult control.
c------
            if(itrigl.eq.1.and.jj.eq.1.and.ipart.gt.1.and.igo.le.2)
     1          locl=locl+1
            if(itrigr.eq.1.and.jj.eq.2.and.ipart.gt.1.and.igo.le.2)
     1          locr=locr+1
c
            do 300 i=1,ipart
                j2=NMD
                if(i.eq.ipart) j2=mod(mode,NMD)
                if(j2.eq.0) j2=NMD
c
                do 200 j=1,j2
                    if(igo.eq.1)then
c-----
c                       READ LOVE
c-----
                        wvml(j,1) = wvml(j,2) 
                        wvmls(j,1) = wvmls(j,2) 
                        wvmlr(j,1) = wvmlr(j,2) 
                        ales(j,1)  = ales(j,2) 
                        aler(j,1)  = aler(j,2) 
                        gaml(j,1)  = gaml(j,2) 
                        ut(j,1)    = ut(j,2) 
                        ut0(j,1)   = ut0(j,2) 
                        dut(j,1)   = dut(j,2) 
                        tt0(j,1)   = tt0(j,2) 
                        call getegn(icode,ifunc(icode),1,wvno,u,gamm,
     1                      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2                      rut,rtt,ruz,rtz,rare,wvnrec,rur0,
     3                      sumkr,sumgr,sumgv,ierr)
                        wvml(j,2) = wvno
c-----
c                       Get Love - 1 = normal, 3 = lat2d
c                       actually this is all the same
c-----
                        wvmlr(j,2) = wvnrec
                        wvmls(j,2) = wvnsrc
                        ales(j,2) = sare
                        aler(j,2) = rare

                        gaml(j,2)  = gamm
                        ut(j,2)    = sur
                        dut(j,2)   = sdur
                        ut0(j,2)   = rut
                        tt0(j,2)  =  rtt
                    else if(igo.eq.2)then
c-----
c                       READ RAYLEIGH
c-----
                        wvmr(j,1)  = wvmr(j,2)
                        wvmrs(j,1) = wvmrs(j,2)
                        wvmrr(j,1) = wvmrr(j,2)
                        ares(j,1)  = ares(j,2)
                        arer(j,1)  = arer(j,2)
                        gamr(j,1)  = gamr(j,2)
                        ur(j,1)    = ur(j,2)
                        dur(j,1)   = dur(j,2)
                        uz(j,1)    = uz(j,2)
                        duz(j,1)   = duz(j,2)
                        ur0(j,1)   = ur0(j,2)
                        tr0(j,1)   = tr0(j,2)
                        uz0(j,1)   = uz0(j,2)
                        tz0(j,1)   = tz0(j,2)
                        call getegn(icode,ifunc(icode),1,wvno,u,gamm,
     1                      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2                      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3                      sumkr,sumgr,sumgv,ierr)
CRBH               if(per.gt.9.5.and.per.lt.10.5)then
CRBH                  WRITE(6,*)'icode    :',icode
CRBH                  WRITE(6,*)'ifunc    :',ifunc(icode)
CRBH                  WRITE(6,*)'per      :',per
CRBH                  WRITE(6,*)'wvno     :',wvno
CRBH                  WRITE(6,*)'u        :',u   
CRBH                  WRITE(6,*)'sur      :',sur 
CRBH                  WRITE(6,*)'sdur     :',sdur 
CRBH                  WRITE(6,*)'suz      :',suz 
CRBH                  WRITE(6,*)'sduz     :',sduz 
CRBH                  WRITE(6,*)'sare     :',sare 
CRBH                  WRITE(6,*)'wvnsrc   :',wvnsrc
CRBH                  WRITE(6,*)'sur0     :',sur0 
CRBH                  WRITE(6,*)'-----'
CRBH                  WRITE(6,*)'rur      :',rur 
CRBH                  WRITE(6,*)'rtr      :',rtr 
CRBH                  WRITE(6,*)'ruz      :',ruz 
CRBH                  WRITE(6,*)'rtz      :',rtz 
CRBH                  WRITE(6,*)'rare     :',rare 
CRBH                  WRITE(6,*)'wvnrec   :',wvnrec
CRBH                  WRITE(6,*)'rur0     :',rur0 
CRBH               endif
                        wvmr(j,2) = wvno
c-----
c                       Get Rayleigh - 2 = normal, 4 = lat2d
c                       actually this is all the same
c-----
                        wvmrr(j,2) = wvnrec
                        wvmrs(j,2) = wvnsrc
                        arer(j,2)  = rare
                        ares(j,2)  = sare
                        
                        gamr(j,2)  = gamm
                        ur(j,2)    = sur
                        dur(j,2)   = sdur
                        uz(j,2)    = suz
                        duz(j,2)   = sduz

                        ur0(j,2)   = rur
                        tr0(j,2)   = rtr
                        uz0(j,2)   = ruz
                        tz0(j,2)   = rtz
                    endif
  200           continue
                if(ipart.gt.1) then
                    if(jj.eq.1) then
                        irec=i
                        if(locl.ge.3) locl=1
                        if(locl.eq.2) irec=irec+40
                        pp(1,locl)=per
                        write(7,rec=irec) (wvml(k,2),
     1              ales(k,2),aler(k,2),gaml(k,2),ut(k,2),
     2                      dut(k,2),ut0(k,2),
     3                      tt0(k,2),wvmlr(k,2),
     4                      wvmls(k,2),
     5                      k=1,NMD)
                    endif
                    if(jj.eq.2) then
                        irec=i
                        if(locr.ge.3) locr=1
                        if(locr.eq.2) irec=irec+40
                        pp(2,locr)=per
                        write(8,rec=irec) (wvmr(k,2),
     1              ares(k,2),arer(k,2),gamr(k,2),ur(k,2),
     2                  dur(k,2),uz(k,2),duz(k,2),
     3                  ur0(k,2),tr0(k,2),uz0(k,2),
     4                  tz0(k,2),
     5                  wvmrs(k,2),wvmrr(k,2),
     6                      k=1,NMD)
                    endif
                endif
c-----
c               DONE WITH WAVE TYPE
c-----
  300       continue
        endif
        return
        end
 
        subroutine srclyr(depth,lmax,dph)
        implicit none
        real depth, dph
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer m
        real dep
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c            = 0 is the free surface 
c       depth = source depth 
c       dph = height of  source above lmax + 1 interface 
c       lmax = 0 is the free surface 
c-----
        if(depth.le.0.0)then
            lmax = 1
            dph = 0.0
        else
            dep = 0.0 
            do 100 m = 2,mmax
                dep = dep + d(m-1) 
                dph = dep - depth 
                lmax = m 
                if(abs(dph).lt. 0.0001*d(m-1) .or.
     1              abs(dph).lt.1.0e-6)go to 101
  100       continue 
  101   continue 
        endif
        return 
        end 


        subroutine srcpul(ntau,ipt,alp,dt,npts,rfile,duration)
c-----
c       define the source time function
c-----
c       ntau    I*4 - pulse duration factor for parabolic and triangular
c               pulses
c       ipt I*4 - pulse type
c               0 - triangular
c               1 - parabolic
c               2 - Ohnaka
c               3 - Dirac delta function
c               4 = user pulse in file rfile
c       alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c       dt  R*4 - sample rate for time series
c       npts    I*4 - number of samples in time series
c       rfile   C*80- name of user provided pulse 
c               note dt must be the same as specified in dfile
c       duration R*4 - duration of triangular/parabolic source pulse
c-----
        parameter (LIN=5, LOT=6, LER=0)
        integer*4 ntau, ipt, npts
        real*4 alp, dt
        character*80 rfile
        integer*4 kkr,kkl,kkf,np0
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     1      kkf,itrigl,itrigr,mods,modt,modi
        integer NSAMP
        parameter (NSAMP=16384)
        common/srctim/ src(NSAMP)
        do 50 i=1,NSAMP
            src(i)=0.0
   50   continue
        if(ipt.eq.0)then
c-----
c           triangluar
c-----
            tau = ntau * dt
            nt = 2*ntau + 1
            call pultd(tau,dt,nt,l)
            duration = 2.*ntau*dt
        else if(ipt.eq.1)then
c-----
c           parabolic
c-----
            tau = ntau * dt
            nt = 4*ntau + 1
            call pulpd(tau,dt,nt,l)
            duration = 4.*ntau * dt
        else if(ipt.eq.2)then
c-----
c           Ohnaka
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call  pulod(alp,dt,nt,l)
            duration = 0.0
        else if(ipt.eq.3)then
c-----
c           Dirac
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call puldd(dt,nt,l)
            duration = 0.0
        else if(ipt.eq.4)then
c-----
c           User specified
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call pulud(rfile,nt,tau,ddt)
                if(ddt.ne.dt)then
                    write(LER,*)'Warning: RFILE dt',ddt,
     1              ' is not same as dt',dt,
     1                  ' specified for synthetic'
                endif

            ntau = tau/dt
            duration = 0.0
        endif
        return
        end
 
        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c       If the string str contains an E or e, then
c       we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c       read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

      subroutine domom(xx,uz,ur,ut,duzdz,durdz,dutdz,
     1      duzdr,durdr,dutdr,duzdt,durdt,dutdt)
c-----
c     apply the moment tensor/force to this set of Green's functions
c     at one frequency
c-----
        integer NGRN
        parameter (NGRN=63)
        complex xx(NGRN)
        complex uz,ur,ut,duzdz,durdz,dutdz,
     1      duzdr,durdr,dutdr,duzdt,durdt,dutdt

        common/stresstrain/dostrain,dostress,
     1     az, baz, oxmt(3,3), xmom,oforcex, oforcey,oforcez,
     1     lrot, caz, saz, c2az, s2az
        logical dostrain, dostress
        real az, baz, oxmt, xmom,oforcex, oforcey,oforcez
        real caz, saz, c2az, s2az

        real xmt(3,3)
        real forcex, forcey, forcez
        integer i,j
c-----
c       assume that the model in KM KM/S GM/CM^3
c       and that the force is in dyne and the moment in dyne-cm
c
c-----
        do i=1,3
           do j=1,3
              xmt(i,j) = oxmt(i,j)/1.0e+20
           enddo
        enddo
        forcex = oforcex/1.0e+15
        forcey = oforcey/1.0e+15
        forcez = oforcez/1.0e+15

        uz = forcex*caz*xx(13) + forcey*saz*xx(13) + forcez*xx(11)
     1     + xmt(1,1)*( c2az*xx(6)/2. - xx(1)/6. + xx(9)/3.)
     1     + xmt(2,2)*(-c2az*xx(6)/2. - xx(1)/6. + xx(9)/3.)
     1     + xmt(3,3)*( xx(1)/3. + xx(9)/3.)
     1     + xmt(1,2)*s2az*xx(6)
     1     + xmt(1,3)* caz*xx(3)
     1     + xmt(2,3)* saz*xx(3)
        ur = forcex*caz*xx(14) + forcey*saz*xx(14) + forcez*xx(12)
     1     + xmt(1,1)*( c2az*xx(7)/2. - xx(2)/6. + xx(10)/3.)
     1     + xmt(2,2)*(-c2az*xx(7)/2. - xx(2)/6. + xx(10)/3.)
     1     + xmt(3,3)*( xx(2)/3. + xx(10)/3.)
     1     + xmt(1,2)*s2az*xx(7)
     1     + xmt(1,3)* caz*xx(4)
     1     + xmt(2,3)* saz*xx(4)
c check the force since this is old?
c-----
c       note this does not follow text because of the way that TSS
c       is defined. The difference is a minus sign
c-----
        ut = ( forcex*saz-forcey*caz)*xx(15)
     1     + xmt(1,1)*s2az*xx(8)/2.
     1     - xmt(2,2)*s2az*xx(8)/2.
     1     - xmt(1,2)*c2a*xx(8)/2.
     1     + xmt(1,3)*saz*xx(5)
     1     - xmt(2,3)*caz*xx(5)

        duzdr = forcex*caz*xx(13+21) + forcey*saz*xx(13+21) 
     1     + forcez*xx(11+21)
     1     + xmt(1,1)*( c2az*xx(6+21)/2. - xx(1+21)/6. + xx(9+21)/3.)
     1     + xmt(2,2)*(-c2az*xx(6+21)/2. - xx(1+21)/6. + xx(9+21)/3.)
     1     + xmt(3,3)*( xx(1+21)/3. + xx(9+21)/3.)
     1     + xmt(1,2)*s2az*xx(6+21)
     1     + xmt(1,3)* caz*xx(3+21)
     1     + xmt(2,3)* saz*xx(3+21)
        durdr = forcex*caz*xx(14+21) + forcey*saz*xx(14+21) 
     1     + forcez*xx(12+21)
     1     + xmt(1,1)*( c2az*xx(7+21)/2. - xx(2+21)/6. + xx(10+21)/3.)
     1     + xmt(2,2)*(-c2az*xx(7+21)/2. - xx(2+21)/6. + xx(10+21)/3.)
     1     + xmt(3,3)*( xx(2+21)/3. + xx(10+21)/3.)
     1     + xmt(1,2)*s2az*xx(7+21)
     1     + xmt(1,3)*xx(4+21)*caz
     1     + xmt(2,3)*xx(4+21)*saz
        dutdr = ( forcex*saz-forcey*caz)*xx(15+21)
     1     + xmt(1,1)*s2az*xx(8+21)/2.
     1     - xmt(2,2)*s2az*xx(8+21)/2.
     1     - xmt(1,2)*c2a*xx(8+21)/2.
     1     + xmt(1,3)*saz*xx(5+21)
     1     - xmt(2,3)*caz*xx(5+21)

        duzdz = forcex*caz*xx(13+42) + forcey*saz*xx(13+42) 
     1     + forcez*xx(11+42)
     1     + xmt(1,1)*(c2az*xx(6+42)/2. - xx(1+42)/6. + xx(9+42)/3.)
     1     + xmt(2,2)*(-c2az*xx(6+42)/2. - xx(1+42)/6. + xx(9+42)/3.)
     1     + xmt(3,3)*( xx(1+42)/3. + xx(9+42)/3.)
     1     + xmt(1,2)*s2az*xx(6+42)
     1     + xmt(1,3)*xx(3+42)*caz
     1     + xmt(2,3)*xx(3+42)*saz
        durdz = forcex*caz*xx(14+42) + forcey*saz*xx(14+42) 
     1     + forcez*xx(12+42)
     1     + xmt(1,1)*(c2az*xx(7+42)/2. - xx(2+42)/6. + xx(10+42)/3.)
     1     + xmt(2,2)*(-c2az*xx(7+42)/2. - xx(2+42)/6. + xx(10+42)/3.)
     1     + xmt(3,3)*( xx(2+42)/3. + xx(10+42)/3.)
     1     + xmt(1,2)*s2az*xx(7+42)
     1     + xmt(1,3)*xx(4+42)*caz
     1     + xmt(2,3)*xx(4+42)*saz
        dutdz = ( forcex*saz-forcey*caz)*xx(15+42)
     1     + xmt(1,1)*s2az*xx(8+42)/2.
     1     - xmt(2,2)*s2az*xx(8+42)/2.
     1     - xmt(1,2)*c2az*xx(8+42)/2.
     1     + xmt(1,3)*saz*xx(5+42)
     1     - xmt(2,3)*caz*xx(5+42)

        duzdt = -forcex*saz*xx(13) + forcey*caz*xx(13) 
     1     + xmt(1,1)*(-2*s2az*xx(6)/2. )
     1     + xmt(2,2)*( 2*s2az*xx(6)/2. )
     1     + xmt(1,2)*2*c2az*xx(6)
     1     - xmt(1,3)*xx(3)*saz
     1     + xmt(2,3)*xx(3)*caz
        durdt = -forcex*saz*xx(14) + forcey*caz*xx(14) 
     1     + xmt(1,1)*(-2*s2az*xx(7)/2. )
     1     + xmt(2,2)*( 2*s2az*xx(7)/2. )
     1     + xmt(1,2)*2*c2az*xx(7)
     1     - xmt(1,3)*xx(4)*saz
     1     + xmt(2,3)*xx(4)*caz
        dutdt = ( forcex*caz+forcey*saz)*xx(15)
     1     + xmt(1,1)*2*c2az*xx(8)/2.
     1     - xmt(2,2)*2*c2az*xx(8)/2.
     1     + xmt(1,2)*2*s2az*xx(8)/2.
     1     + xmt(1,3)*caz*xx(5)
     1     + xmt(2,3)*saz*xx(5)

c-----
c     The CPS convention is the Z+ = up
c     For computation of strain we wand Z+ down
c     this means that all dUz must have sign changed
c-----
        duzdz = - duzdz
        duzdr = - duzdr
        duzdt = - duzdt



        return
        end
        subroutine gcmdln(ntau,ipt,alp,dfile,idva,
     1      xmult,rfile,ostr,dolat,dodble,nmode,
     2      dolock,dozero,
     1      dostrain,dostress,az,baz,lrot,xmom,xmt,fx,fy,fz,
     1      dotest1,dostep,outfmt)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c
c       ntau    I*4 - pulse duration factor for parabolic and triangular
c               pulses
c       ipt I*4 - pulse type
c               0 - triangular
c               1 - parabolic
c               2 - Ohnaka
c               3 - Dirac delta function
c               4 = user pulse in file rfile
c       alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c       dfile   C*80- name of distance file
c       idva    I*4 - time history type
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       xmult   R*4 - moment scaling factor
c       rfile   C*80- name of user provided pulse 
c               note dt must be the same as specified in dfile
c       dolat   L   .true. do laterally varying
c       nmode   I*4 - mode number to plot if >= 0
c               -  0 do fundamental only
c               - -1 do all
c               - -2 do fundamental only
c               - -3 do all higher
c       dolock  L   .true. locked mode used so do not use
c                   bottom layer for first arrival computation
c               .false.  (default)
c       dozero  L   .true. for triangular and parabolic pulses, center
c               at zero lag
c       dostrain  - L .true. output strain
c       dostress  - L .true. outpu stress
c       dostep  L   .true.  source time function is step-like
c                         whose derivative is the triangular,parabolic
c                         Ohnaka or impulse
c                   .false. source time function is impulsive.
c                         The impulse is approximated by the triangular,
c                         parabolic, Ohnaka or impulse function
c                         When -STEP -D are given this is the Green's 
c                         function.
c-----
        character*(*) dfile, rfile
        integer*4 ntau, ipt, idva
        real*4 xmult, alp
        logical  dolat, dodble, dolock, dozero, dotest1,dostep
        logical dostrain, dostress
        real az, baz, xmt(3,3), xmom,forcex, forcey,forcez
        logical lrot
        integer outfmt
        character ostr*(*)

        integer*4 mnmarg
        character*80 name
        integer i,j
        logical lmij, lex, lsdr, lforce
        integer lgstr
c-----
c       initialization
c-----
        ntau = -1
        ipt = -1
        alp = -1.0
        idva = 1
        rfile = ' '
        xmult = 1.0
        dfile = ' '
        ostr = ' '
        dolat = .false.
        dodble = .false.
        nmode = -1
        dolock = .false.
        dozero = .false.
        lrot = .true.
        xmom = 1.0
        stk = 0
        dip = 0
        rake = 0
        xmom = 1.0
        lsdr   = .false.
        lmij   = .false.
        lex    = .false.
        lforce = .false.
        dostrain =.false.
        dostress = .false.
        dotest1 = .false.
        dostep = .true.
        outfmt = 1
        do i=1,3
           do j=1,3
              xmt(i,j) = 0.0
           enddo
        enddo
c-----
c     process the command line arguments
c----
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-l')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')ntau
            elseif(name(1:2).eq.'-d')then
                i=i+1
                call mgtarg(i,dfile)
            else if(name(1:2).eq.'-a')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')alp
            else if(name(1:2).eq.'-t')then
                ipt = 0
            else if(name(1:2).eq.'-p')then
                ipt = 1
            else if(name(1:2).eq.'-o')then
                ipt = 2
            else if(name(1:2).eq.'-i')then
                ipt = 3
            else if(name(1:2).eq.'-2')then
                dodble = .true.
            else if(name(1:2).eq.'-Z')then
                dozero = .true.
            else if(name(1:2).eq.'-F' .and. name(1:3).ne.'-FU'
     1         .and. name(1:3).ne.'-FX'
     1         .and. name(1:3).ne.'-FY'
     1         .and. name(1:3).ne.'-FZ'
     1         .and. name(1:3).ne.'-FM'
     1          )then
                ipt = 4
                        i=i+1
                        call mgtarg(i,rfile)
            else if(name(1:4).eq.'-FMT')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10.0)')outfmt
                if(outfmt.lt.1 .or. outfmt.gt.5)then
                   outfmt = 4
                endif
            else if(name(1:2).eq.'-D' .and.
     1          name(1:3).ne.'-DI')then
                idva = 0
            else if(name(1:2).eq.'-V')then
                idva = 1
            else if(name(1:2).eq.'-A' .and. name(1:4).ne.'-ALL'
     1          .and. name(1:3).ne.'-AZ')then
                idva = 2
            else if(name(1:2).eq.'-m')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')xmult
            else if(name(1:4).eq.'-LAT')then
                dolat = .true.
            else if(name(1:4).eq.'-LOC')then
                dolock = .true.
            else if(name(1:2).eq.'-M' .and.
     1               name(1:3).ne.'-MW' .and.
     1               name(1:3).ne.'-M0')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nmode
                if(nmode.ge.0)nmode = nmode + 1
            else if(name(1:5).eq.'-STEP')then
                dostep = .true.
            else if(name(1:4).eq.'-IMP')then
                dostep = .false.
            else if(name(1:3).eq.'-FU')then
                nmode = -2 
            else if(name(1:3).eq.'-HI')then
                nmode = -3 
            else if(name(1:4).eq.'-DIP')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,dip)
                lsdr = .true.
            else if(name(1:4).eq.'-STK')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,strike)
                lsdr = .true.
            else if(name(1:5).eq.'-RAKE')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,rake)
                lsdr = .true.
            else if(name(1:3).eq.'-M0' .or. name(1:3).eq.'-MO')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
            else if(name(1:3).eq.'-MW' .or. name(1:3).eq.'-Mw')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
                xmom = 10.**(1.5*xmom + 16.10)
            else if(name(1:3).eq.'-AZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,az)
            else if(name(1:4).eq.'-BAZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,baz)
            else if(name(1:3).eq.'-EX')then
                isds = 2
            else if(name(1:3).eq.'-RO')then
                lrot = .true.
            else if(name(1:3).eq.'-xx' .or. name(1:3).eq.'-XX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,1))
                lmij = .true.
            else if(name(1:3).eq.'-yy' .or. name(1:3).eq.'-YY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(2,2))
                lmij = .true.
            else if(name(1:3).eq.'-zz' .or. name(1:3).eq.'-ZZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(3,3))
                lmij = .true.
            else if(name(1:3).eq.'-xy' .or. name(1:3).eq.'-XY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,2))
                xmt(2,1) = xmt(1,2)
                lmij = .true.
            else if(name(1:3).eq.'-xz' .or. name(1:3).eq.'-XZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(1,3))
                xmt(3,1) = xmt(1,3)
                lmij = .true.
            else if(name(1:3).eq.'-yz' .or. name(1:3).eq.'-YZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmt(2,3))
                xmt(3,2) = xmt(2,3)
                lmij = .true.
            else if(name(1:3).eq.'-fx' .or. name(1:3).eq.'-FX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fx)
                lforce = .true.
            else if(name(1:3).eq.'-fy' .or. name(1:3).eq.'-FY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fy)
                lforce = .true.
            else if(name(1:3).eq.'-fz' .or. name(1:3).eq.'-FZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fz)
                lforce = .true.
            else if(name(1:7).eq.'-STRAIN')then
                dostrain = .true.
            else if(name(1:7).eq.'-STRESS')then
                dostress =.true.
            else if(name(1:6).eq.'-TEST1')then
                dotest1 = .true.
            else if(name(1:2).eq.'-?')then
                call usage('Help')
            else if(name(1:2).eq.'-h')then
                call usage('Help')
            endif
        go to 11
   13   continue
c-----
c     do some error checking, e.g., we cannot generate velocity
c     for triangular pulse
c-----
        if(ipt.ge.0 .and. ipt .le.1 .and . ntau.le.0)ntau = 1
        if(ipt.eq.2 .and. alp .le.0.0)alp = 1.0
        if(dfile .eq. ' ')call usage('No distance control data file')
        if(ipt.eq.2 .and. alp.lt.0.0)
     1      call usage('No alpha for Ohnaka pulse')
        if(ipt.lt.0)
     1      call usage('No pulse shape defined')
        lostr = len(ostr)
        ostr = 'strainspulse96'
        do 14 i=1,nmarg
            call getarg(i,name)
            L = lgstr(name)
            lo = lgstr(ostr)
            lomax = lo + 1 + L
            if(lomax .lt. lostr)then
                ostr = ostr(1:lo)//' '//name(1:L)
            endif
   14   continue
c-----
c       Priorties for the source specification
c       Mij > Strike, Dip, Slip >  Explosion > Force
c       This means that if Mij  is specified, the others are ignored
c       If Mij not given, but strike,dip,slip are this overrides the -EX
c
c       not if anything if Mij turn off the force
c-----
        if(.not.lmij)then
           if(lsdr)then
              call tensor(stk,dip,rake,xmom,xmt,.true. )
              WRITE(6,*)' Mij not given, Strike,Dip,Rake, no force'
              fx = 0.0
              fy = 0.0
              fz = 0.0
           else
             if(lex)then
                call tensor(stk,dip,rake,xmom,xmt,.false.)
                WRITE(6,*)' Mij not given, Explosion, no force'
                fx = 0.0
                fy = 0.0
                fz = 0.0
             else
                if(lforce)then
                   WRITE(6,*)'Force only'
                endif
             endif
           endif
        else
c-----
c          turn off parameters that are not related to MT source
c-----     
           WRITE(6,*)' Mij given no force'
           fx = 0.0
           fy = 0.0
           fz = 0.0
        endif
 


        return
        end

        subroutine gtegn(rr,tshft,ipar)
c-----
c       set the end values for interpolation.
c-----
        integer ipar(10)
        parameter (LIN=5, LOT=6, LER=0)
        parameter (NGRN=63)
        parameter (NMD=100)
        integer*4 itrig(2)
        integer*4 kkr,kkl,kkf,ifunc(2),mode2,np0
        integer*4 npx,isign
        integer NSAMP
        parameter (NSAMP=16384)
        common/srcspc/ssrc(NSAMP)
        complex ssrc
        common/rayl/
     1      ur(NMD,2),dur(NMD,2),uz(NMD,2),duz(NMD,2),
     2      ur0(NMD,2),uz0(NMD,2),tr0(NMD,2),tz0(NMD,2),
     3      wvmr(NMD,2),ares(NMD,2),arer(NMD,2),gamr(NMD,2),
     4      wvmrr(NMD,2), wvmrs(NMD,2)
        common/love/
     1      ut(NMD,2),dut(NMD,2),ut0(NMD,2),tt0(NMD,2),
     2      wvml(NMD,2),ales(NMD,2),aler(NMD,2),gaml(NMD,2),
     3      wvmlr(NMD,2), wvmls(NMD,2)
        common/ctrl/   sr,si,xmom,zeta,ms,np0,kkr,kkl,
     *      kkf,itrigl,itrigr,mods,modt,modi
        common/resp/   df,bmax,modes1(2),modes2(2),pers1(2),pers2(2),
     *      ipart1(2),ipart2(2)
        complex xx(NGRN)
        np2=np0/2+1
c------
c     connect with the input data files.
c------
        call readin(0,jj,ifunc,mode2,ipart,per)
        do 150 jj=1,2
            j=jj
            if(jj.eq.1.and.kkl.ne.1) go to 150
            if(jj.eq.2.and.kkr.ne.1) go to 150
            itrig(j)=0
            igo=j
            call readin(igo,j,ifunc,mode2,ipar2,per2)
            modes2(j)=mode2
            ipart2(j)=ipar2
            pers2(j)=per2
  150   continue
c------
c     start the job.
c------
        do 600 ii=np2,1,-1
            do 200 i=1,NGRN
                xx(i)=cmplx(0.0,0.0)
  200       continue
            per=-1.0
            if(ms.ne.1) then
                    sr = real(ssrc(ii))
                    si = aimag(ssrc(ii))
            endif
            if(ii.eq.1) then
                per=1.e+5
                if(ii.eq.1) go to 320
            endif
            per=1./((ii-1)*df)
            if(ii.eq.2) go to 320
            dper=1./((ii-2)*df)-per
  320       continue
            per0=per-0.005*dper
            if(ii.eq.np2) per0=per+0.005*dper
c
            do 500 jj=1,2
                j=jj
                if(jj.eq.1.and.kkl.ne.1) go to 500
                if(jj.eq.2.and.kkr.ne.1) go to 500
c------
c     find the data at two end points. 
c            Interpolation is done between them.
c------
  330           if(per0.lt.pers2(j)) go to 400
                if(itrig(j).eq.2) go to 500
                pers1(j)=pers2(j)
                modes1(j)=modes2(j)
                ipart1(j)=ipart2(j)
                igo=j
                call readin(igo,j,ifunc,mode2,ipar2,per2)
                modes2(j)=mode2
                ipart2(j)=ipar2
                pers2(j)=per2
                itrig(j)=1
                if(ifunc(igo).le.0) itrig(j)=2
                go to 330
  400           if(itrig(j).eq.0) go to 500
                mode0=modes2(j)
                if(mode0.gt.modes1(j)) mode0=modes1(j)
                if(jj.eq.1) then
                     call excitl(per,mode0,rr,tshft,xx,ipar)
                else if(jj.eq.2) then
                     call excitr(per,mode0,rr,tshft,xx,ipar)
                endif
  500       continue
c------
c     output
c     16 complex values correspond to ten Green s functions:
c       01: ZDD     02: RDD     17: PDD
c       03: ZDS     04: RDS     18: PDS
c       05: TDS     06: ZSS     19: PSS
c       07: RSS     08: TSS     20: PVF
c       09: ZEX     10: REX     21: PHF
c       11: ZVF     12: RVF
c       13: ZHF     14: RHF
c       15: THF     16: PEX
c------
            write(10) (xx(i),i=1,NGRN)
  600   continue
        return
C       write(LOT,*) 'Period not agree between source and receiver'
C       write(LOT,*) 'files at period=',per2,perr2
C       stop
        end



        subroutine tensor(stk,dip,rake,xmom,m,lsdr)
c-----
c       calculate moment tensor for a double couple mechanism
c       or for an explostion
c-----
c       stk - R*4   - strike, measured clockwise from north
c                 when looking down, e.g., east = 90
c       dip - R*4   - dip, measured from horizontal when looking
c                 in direction of strike, fault dips down to right
c                 0 <= dip <= 90
c       rake    - R*4   - rake. Measured from horizontal with respect
c                 to strike. If -180 < rake <= 180 is taken as
c                 the convention, then a negative rake, implies 
c                 that the P-wave first motion in the center of the
c                 focal sphere is negative (e.g., dilatation)
c       xmom    - R*4   - moment value
c       m(3,3)  - R*4   - moment tensor
c       lsdr    - L     - .true. get moment tensor for strike,dip,rake
c                         else get moment tensor for isotropic explosion
c-----
        integer LOT
        parameter (LOT=6)
        real*4 stk, dip, rake, xmom, m(3,3)
        logical lsdr

        if(lsdr)then

            degrad=0.0174532925
            tol = 1.0e-7
            dp = degrad*dip
            st = degrad*stk
            sl = degrad*rake
            sind=sin(dp)
            cosd=cos(dp)
            sinr=sin(sl)
            cosr=cos(sl)
            sins=sin(st)
            coss=cos(st)
            sin2d=sin(2.*dp)
            cos2d=cos(2.*dp)
            sin2s=sin(2.*st)
            cos2s=cos(2.*st)
            m(1,1)=(-sind*cosr*sin2s-sin2d*sinr*sins*sins)*xmom
            m(2,2)=(sind*cosr*sin2s-sin2d*sinr*coss*coss)*xmom
            m(3,3)=(sin2d*sinr)*xmom
            m(1,2)=(sind*cosr*cos2s+0.5*sin2d*sinr*sin2s)*xmom
            m(1,3)=(-cosd*cosr*coss-cos2d*sinr*sins)*xmom
            m(2,3)=(-cosd*cosr*sins+cos2d*sinr*coss)*xmom
            m(2,1) = m(1,2)
            m(3,1) = m(1,3)
            m(3,2) = m(2,3)
c-----
c           clean up small values
c-----
            xmax=-1.0e+37
            do  i=1,3
                do  j=1,3
                    if(abs(m(i,j)).gt.xmax)xmax = abs(m(i,j))
                enddo
            enddo
            thresh = tol * xmax
            do  i=1,3
                do  j=1,3
                    if(abs(m(i,j)).lt.thresh) m(i,j) = 0.0
                enddo
            enddo

c-----
c           write out the information
c-----
            write(LOT,*)' STRIKE    =',stk
            write(LOT,*)' DIP       =',dip
            write(LOT,*)' RAKE      =',rake
            write(LOT,*)' MOMENT    =',xmom

        else
            m(1,1)=xmom
            m(2,2)=xmom
            m(3,3)=xmom
            m(1,2)=0.0
            m(1,3)=0.0
            m(2,3)=0.0
            m(2,1) = m(1,2)
            m(3,1) = m(1,3)
            m(3,2) = m(2,3)
            write(LOT,*)' ISOTROPIC SOURCE'
            write(LOT,*)' MOMENT    =',xmom
        endif
      
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'spulse96:',str
        write(LER,*)'USAGE: ',
     1  'strainspulse96 -d Distance_File  [ -t -o -p -i ] [-a alpha]',
     2  ' -l L [ -D -V -A]  [-F rfile ] [ -m mult] [-STEP|-IMP]',
     3  '  [-FUND] [-HIGH] [-Z] ',
     4  ' [-LAT] [-2] [ -M mode ] [-LOCK] -FMT ifmt '
        write(LER,*)
     1   ' [-M0 moment ] [-MW mw] [-STK stk -DIP dip -RAKE rake]',
     2   ' [-FX fx -FY fy -FZ fz] ',
     3   ' [-XX Mxx ... -ZZ Mzz] [-?] [-h]'
        write(LER,*)'TIME FUNCTION SPECIFICATION'
        write(LER,*)
     1  ' -t           Triangular pulse of base 2 L dt'
        write(LER,*)
     1  ' -p           Parabolic Pulse of base  4 L dt'
        write(LER,*)
     1  '              -p -l 1 recommended'
        write(LER,*)
     1  ' -l L         (default 1 )duration control parameter'
        write(LER,*)
     1  ' -o           Ohnaka pulse with parameter alpha'
        write(LER,*)
     1  ' -i           Dirac Delta function'
        write(LER,*)
     1  ' -a alpha     Shape parameter for Ohnaka pulse'
        write(LER,*)
     1  ' -F rfile     User supplied pulse'
        write(LER,*)
     1  ' -m mult      Multiplier (default 1.0)'
        write(LER,*)
     1  ' -STEP        (default)'
        write(LER,*)
     1  ' -IMP         '
        write(LER,*)
     1  '              By default the source time function is '
        write(LER,*)
     1  '              steplike. -IMP forces impulse like. -D -IMP ',
     1  'is Green s function'
        write(LER,*)'OUTPUT FILE NAME'
        write(LER,'(a/a/a/a/a/a/a/a/a)')
     1  '  The format for the name of the binary output attempts to',
     1  '  give information on epicentral distance (km), ',
     1  '  source depth (km), and receiver depth(km). The options are',
     1  '  -FMT 1      DDDDDd_HHHh_ZZZz.cmp',
     1  '              e.g. 005001_1234_0045.Uz',
     1  '  -FMT 2      DDDDDddd_HHHhhh_ZZZzzz.cmp',
     1  '              e.g. 00500123_123456_004578.Erf',
     1  '  -FMT 3      DDDDDdHHHh.grn(default)',
     1  '              e.g. 0050010041.ZVF'
        write(LER,'(a/a/a/a/a/a/a/a/a)')
     1  '  -FMT 4      DDDDdHHHh.grn',
     1  '              e.g. 050010045.Srz',
     1  '  -FMT 5      DDDdddHhhh.grn',
     1  '              e.g. 5001234578.Err',
     1  '  where D is for epicentral distance, H source depth, and',
     1  '  Z receiver depth. The lower case indicates the digits ',
     1  '  to the right of the decimal place.  The examples above',
     1  '  are for an epicentral distance is 500.123 km, source',
     1  '  depth 123.456 km and receiver depth 4.578 km.'

        WRITE(LER,*)'OUTPUT TIMESERIES FOR SOURCE as Ur, Ut, Uz',
     1  ' components with strain, stress optional'
        write(LER,*)
     1  ' -D           Output is ground displacement        (m)'
        write(LER,*)
     1  ' -V           Output is ground velocity (default) (m/s)'
        write(LER,*)
     1  ' -A           Output is ground acceleration       (m/s^2)'
        write(LER,'(a/a)')
     1  '  -STRESS (default .false. ) output stress for mechanism',
     1  '    units are Pa, with suffix Srr, Srf, Srz, Stt, Sfz, Szz'
        write(LER,'(a/a)')
     1  '  -STRAIN (default .false. ) output strain for mechanism',
     1  '       with suffix, Err, Erf, Erz, Eff, Efz, Ezz'
        write(LER,'(a/a/a/a/a)')
     1  '  -TEST1  (default .false.) output CPS Green functions ,e.g.,',
     1  '       ZDS RDS ... RHF THF for use with moment tensor codes',
     1  '       and gsac MT command. This is equivalent to ',
     1  '       spulse96 -V -p -l 1 | f96tosac -G if -FMT 4 is used ',
     1  '       with strainspulse96'

        write(LER,*)'COMPUTATIONS'
        write(LER,*)
     1  ' -d Distance_File {required}    Distance control file '
        write(LER,*)
     1  '   This contains one of more lines with following entries'
        write(LER,*)
     1  '   DIST(km) DT(sec) NPTS T0(sec) VRED(km/s) '
        write(LER,*)
     2  '             first time point is T0 + DIST/VRED'
        write(LER,*)
     3  '             VRED=0 means do not use reduced travel time, e.g.'
        write(LER,*)
     1  '             500.0 0.25 512 -23.33 6.0'
        write(LER,*)
     1  '             500.0 0.25 512  60    0.0 '
        write(LER,*)
     1  '             both have first sample at travel time of 60s'
        write(LER,*)
     1  ' -LAT        (default false) Laterally varying eigenfunctions'
        write(LER,*)
     1  ' -2          (default false) Use double length  internally'
        write(LER,*)
     1  ' -M  nmode   (default all) mode number to ',
     2  'compute[0=fund,1=1st]'
        write(LER,*)
     1  ' -Z          (default false) zero phase ',
     2  'triangular/parabolic pulse'
        write(LER,*)
     1  ' -FUND       (default all) fundamental modes only  '
        write(LER,*)
     1  ' -HIGH       (default all) all higher modes only  '
        write(LER,*)
     1  ' -LOCK       (default false) locked mode used  '
c------
c      strain and stress
c------
        write(LER,*)'SOURCE MECHANISM SPECIFICATION'
        write(LER,*)
     1  ' -DIP dip               dip of fault plane'
        write(LER,*)
     1  ' -STK Strike            strike of fault plane'
        write(LER,*)
     1  ' -RAKE Rake              slip angle on fault plane'
        write(LER,*)
     1  ' -M0 Moment (def=1.0) Seismic moment in units of dyne-cm'
        write(LER,*)
     1  ' -MW mw            Moment Magnitude  '
        write(LER,*)
     1  '            moment (dyne-cm) from log10 Mom = 16.10 + 1.5 Mw'
        write(LER,*)
     1  '              For strike,dip,rake source mw or Moment must',
     1  ' be specified'
        write(LER,*)
     1  ' -EX                  Explosion'
        write(LER,*)
     1  ' -AZ Az                Source to Station Azimuth'
        write(LER,*)
     1  ' -BAZ Baz               Station to Source azimuth'
        write(LER,*)
     1  ' -fx FX -fy Fy -fZ fz  Point force amplitudes ',
     2  ' (N,E,down) in  dynes'
        write(LER,*)
     1  ' -XX Mxx -YY Myy -ZZ Mzz  Moment tensor elements in units of'
        write(LER,*)
     1  ' -XY Mxy -XZ Mxz -YZ Myz    dyne-cm'
        write(LER,*)
     1  ' The moment tensor coordinates are typically X = north',
     2  ' Y = east and Z = down'
        write(LER,*)
     1  ' If by accident more than one source specification is used,'
        write(LER,*)
     1  ' the hierarchy is Mij > Strike,dip,rake > Explosion > Force'
        write(LER,*)
     1  '--------------------------------------------------------------'
        write(LER,*)
     1  'NOTE: The output units are related tot he model specification.'
        write(LER,*)
     1  'To have the desired units the model must be in KM, KM/S',
     1  '  and GM/CM^3'
        write(LER,*)
     1  '--------------------------------------------------------------'

        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end

        subroutine output(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,kkl,kkr,ostr,
     2      mname,ipar,fpar,dodble,dolock,dolat,dostep,dostrain)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       kkl I*4 - = 1 indicates Love wave eigenfunctions
c       kkr I*4 - = 1 indicates Rayleigh wave eigenfunctions
c       ostr    C*(*)   - command line string
c       mname   C*(*)   - name of model file
c       ipar    I*4 - arrau of integer command parameters from sdisp96
c       fpar    R*4 - arrau of float command parameters from sdisp96
c       dodble  L   - if .true. double length of time series for FFT
c       dolock  L   - if .true. do not use bottom layer for 
c                   first arrival computation
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green s function
c               21 - Green s function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
c
c-----

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      tp, tsv, tsh,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        real*4 tp, tsv, tsh
        real*4 sa, sc, sf, sl, sn, sr

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax

        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=63)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr
        character ostr*(*)

        integer kerr, system

        integer NSAMP2
        parameter(NSAMP2=16384)
        common/dd/datc
        complex datc(NSAMP2)
        common/xxx/x
        real*4 x(NSAMP2)
        complex xx(NGRN)
        complex z
         integer iszrt(21)
        character*8 ostnm(2)
        character*8 ocmpnm(9)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
        logical dodble, dolock, dolat, dostep, dostrain

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dufdr,duzdf,durdf,dufdf

        complex err,erf,erz,eff,efz,ezz
        complex xout(9)

        character ofile*80
        integer ls
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm/'Uz      ','Ur      ','Ut      ',
     1       'Err     ', 'Erf     ',
     1       'Erz     ', 'Eff     ', 'Efz     ', 'Ezz     '/


        tau = ntau * dt
        los = lgstr(ostr)
        cpulse = ostr(1:los)
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
c-----
c       get first arrival time - NOTE for implementation of
c-----
            call frstar(dist,depths,depthr,mname,1,tp ,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            call frstar(dist,depths,depthr,mname,2,tsv,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            call frstar(dist,depths,depthr,mname,3,tsh,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
c-----
c       get the source A, C, F, L, N and density values
c-----
            sa = vsr*vsa*vsa
            sc = vsr*vsa*vsa
            sl = vsr*vsb*vsb
            sn = vsr*vsb*vsb
            sf = sa - 2.0*sl
            sr = vsr
c-----
c       lat2d96 we cannot easily predict the first arrival time
c       of a 2-D model. So we set this value to the default for
c       sac, e.g., -12345
c       JUST IGNORE PREVIOUSE - however we must preserve the sa sc cl sn sf sr
c-----
            if(dolat)then
                tp    = -12345.
                tsv   = -12345.
                tsh   = -12345.
            endif
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c       there will be a total of 9, e.g.,
c       Uz  Ur Ut Err Ert Erz Ett Etz Ezz
c-----
c-----
c       jout   1  Uz
c              2  Ur
c              3  Ut
c              4  Err
c              5  Ert
c              6  Erz
c              7  Ett
c              8  Etz
c              9  Ezz
c----
c       This is a two pass process
c       Starting with  ur uz tz durdz and other derivatives on LUN=10
c       compute ur uz ut and strains an then output them in LUN=9
c       
c-----
        open(9,file='file9.tmp',access='sequential', 
     1   form='unformatted', status='unknown')
        rewind 9
        rewind 10
        do  k=np2,1,-1
              read(10) (xx(i),i=1,NGRN)
              call domom(xx,uz,ur,ut,duzdz,durdz,dufdz,
     1            duzdr,durdr,dufdr,duzdf,durdf,dufdf)
c-----
c     dist is in km conert to meters
c     the following is the strain formula
c     be careful of units especially the eff
c     Some dreivatives are dU/dr or dU/dz obviously
c     we must convert from km to m
c-----
              r = dist
              err = durdr
              eff = (dufdf + ur)/r
              ezz = duzdz
              erf = 0.5*(durdf/r  + dufdr - uf/r)
              erz = 0.5*(duzdr + durdz)
              efz = 0.5*(dufdz + duzdf/r)
c-----
c     CONVERT FROM CM TO M
c-----
             xout(1) = uz/1.0e+02
             xout(2) = ur/1.0e+02
             xout(3) = ut/1.0e+02
c-----
c     CONVERT FROM KM TO M for distance and CM to M for displacement
c-----
             xout(4) = err/1.0e+05
             xout(5) = erf/1.0e+05
             xout(6) = erz/1.0e+05
             xout(7) = eff/1.0e+05
             xout(8) = efz/1.0e+05
             xout(9) = ezz/1.0e+05
              write(9)(xout(i),i=1,9)
         enddo
c-----
c     now make 9 passes to make the synthetics
c     for ur uz ut err ert erz ett etz ezz
c-----
      if(dostrain)then
         jkup = 9
      else
         jkup = 3
      endif
         rewind 9
         do 1300 jk=1,jkup
                rewind 9
                do k=np2,1,-1
                    read(9) (xout(i),i=1,9)
                    datc(k)=xout(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                if(ipar(1).gt.0)then
                    rad = distdg*3.1415927/180.0
                    fac = sqrt(dist/(6371.0*sin(rad)) )
                else
                    fac = 1.0
                endif
                    call zfour(datc,nptin,+1,dt,df)
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
 1303               continue
                if(dodble)then
                    npts = nptin / 2
                else
                    npts = nptin
                endif
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                IF(JK.LE.3)THEN
                    if(dostep)then
                       if(idva.eq.0)then
                           call integ(x,npts,dt,.false.)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                       endif
                    else
                       if(idva.eq.1)then
                           call deriv(x,npts,dt)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                           call deriv(x,npts,dt)
                       endif
                    endif
                else
                    if(dostep)then
                        call integ(x,npts,dt,.false.)
                    endif
                endif
c-----
c           vertical
c-----
                if(jk.eq.1)then
                    cmpinc = -90.0 
                    cmpaz  =   0.0
c-----
c           radial
c-----
                else if(jk.eq.2)then
                    cmpinc = 0.0 
                    cmpaz  =   0.0
c-----
c           transverse
c-----
                else if(jk.eq.5)then
                    cmpinc = 0.0 
                    cmpaz  =  90.0
                else
                    cmpinc = -12345.
                    cmpaz  = -12345.
                endif
                stcomp = 'SYN     '
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
c-----
c
c-----
                stcomp = ocmpnm(jk)(1:3)
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,i4.4,'.',a3)
                 endif
                 ls = lgstr(ofile)
         WRITE(0,*)'ofile:',ofile(1:ls)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setkhv('KT0     ','S       ',nerr)
                    call setkhv('KT1     ','S       ',nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',b + indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',b + indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
c----
c       convert from SEED to SAC component incidence
c-----
                cmpinc = cmpinc + 90.0
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','SYN     ' ,nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',stcomp ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                    call bwsac(12,npts,ofile(1:ls),x)
 1300   continue
        close(9,status='delete')
        return
        end

        subroutine outtest1(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,kkl,kkr,ostr,
     2      mname,ipar,fpar,dodble,dolock,dolat,dostep)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       kkl I*4 - = 1 indicates Love wave eigenfunctions
c       kkr I*4 - = 1 indicates Rayleigh wave eigenfunctions
c       ostr    C*(*)   - command line string
c       mname   C*(*)   - name of model file
c       ipar    I*4 - arrau of integer command parameters from sdisp96
c       fpar    R*4 - arrau of float command parameters from sdisp96
c       dodble  L   - if .true. double length of time series for FFT
c       dolock  L   - if .true. do not use bottom layer for 
c                   first arrival computation
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green s function
c               21 - Green s function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
c
c-----

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      tp, tsv, tsh,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        real*4 tp, tsv, tsh
        real*4 sa, sc, sf, sl, sn, sr

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax

        common/receiver/l2mu,mu
        real l2mu, mu

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=63)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr
        character ostr*(*)

        common/grnfmt/outfmt
        integer outfmt

        integer kerr, system

        integer NSAMP2
        parameter(NSAMP2=16384)
        common/dd/datc
        complex datc(NSAMP2)
        common/xxx/x
        real*4 x(NSAMP2)
        complex xx(NGRN)
        complex z
         integer iszrt(21)
        character*8 ostnm(2)
        character*8 ocmpnm1(15)
        character*3 ocmpnm2(3)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
        logical dodble, dolock, dolat, dostep

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dutdr,duzdt,durdt,dufdt

        complex err,ert,erz,ett,etz,ezz
        complex xout(9)

        character ofile*80
        character oofile*80
        integer ls, lgstr
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm1/'ZDD     ','RDD     ',
     1       'ZDS     ','RDS     ','TDS     ',
     1       'ZSS     ','RSS     ','TSS     ',
     1       'ZEX     ','REX     ',
     1       'ZVF     ','RVF     ',
     1       'ZHF     ','RHF     ','THF     '/
        data ocmpnm2/'   ','_dr','_dz'/


        tau = ntau * dt
        los = lgstr(ostr)
        cpulse = ostr(1:los)
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
c-----
c       get first arrival time - NOTE for implementation of
c-----
            call frstar(dist,depths,depthr,mname,1,tp ,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            call frstar(dist,depths,depthr,mname,2,tsv,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            call frstar(dist,depths,depthr,mname,3,tsh,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
c-----
c       get the source A, C, F, L, N and density values
c-----
            sa = vsr*vsa*vsa
            sc = vsr*vsa*vsa
            sl = vsr*vsb*vsb
            sn = vsr*vsb*vsb
            sf = sa - 2.0*sl
            sr = vsr
c-----
c       lat2d96 we cannot easily predict the first arrival time
c       of a 2-D model. So we set this value to the default for
c       sac, e.g., -12345
c       JUST IGNORE PREVIOUSE - however we must preserve the sa sc cl sn sf sr
c-----
            if(dolat)then
                tp    = -12345.
                tsv   = -12345.
                tsh   = -12345.
            endif
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c       there will be a total of 9, e.g.,
c       Uz  Ur Ut Err Ert Erz Ett Etz Ezz
c-----
c-----
c       jout   1  Uz
c              2  Ur
c              3  Ut
c              4  Err
c              5  Ert
c              6  Erz
c              7  Ett
c              8  Etz
c              9  Ezz
c----
c-----
c     now make 9 passes to make the synthetics
c     for ur uz ut err ert erz ett etz ezz
c-----
         do 1300 jk1=1,15
              do 1301 jk2=1,3
                jk = (jk2-1)*21 + jk1
                rewind 10
                do k=np2,1,-1
                    read(10) (xx(i),i=1,63)
                    datc(k)=xx(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                if(ipar(1).gt.0)then
                    rad = distdg*3.1415927/180.0
                    fac = sqrt(dist/(6371.0*sin(rad)) )
                else
                    fac = 1.0
                endif
                    call zfour(datc,nptin,+1,dt,df)
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
 1303               continue
                if(dodble)then
                    npts = nptin / 2
                else
                    npts = nptin
                endif
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
C               if(jk.le.3)then
                IF(JK2.LE.1)THEN
                    if(idva.eq.0)then
                        call integ(x,npts,dt,.false.)
                    else if(idva.eq.2)then
                        call deriv(x,npts,dt)
                    endif
                else
                    call integ(x,npts,dt,.false.)
                endif
                    cmpinc = -12345.
                    cmpaz  = -12345.
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
c-----
c
c-----
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,'_',i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,'_',i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,'_',i4.4,'.',a3)
                 endif
                ls = lgstr(ofile)
                oofile=ofile(1:ls)//ocmpnm1(jk1)(1:3)//ocmpnm2(jk2)(1:3)
                
         WRITE(0,*)'oofile:',oofile(1:ls)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setkhv('KT0     ','S       ',nerr)
                    call setkhv('KT1     ','S       ',nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',b + indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',b + indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
c----
c       convert from SEED to SAC component incidence
c-----
                cmpinc = cmpinc + 90.0
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ',ocmpnm1(jk1),nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',ocmpnm2(jk2) ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                ls = lgstr(oofile)
                    call bwsac(12,npts,oofile(1:ls),x)
 1301   continue
 1300   continue
        return
        end


        subroutine outstress(ntau,idva,dist,dt,nptin,
     1      tfirst,depthr,depths,kkl,kkr,ostr,
     2      mname,ipar,fpar,dodble,dolock,dolat,dostep)
c-----
c       write the Green s functions for this distance
c-----
c       ntau    I*4 - pulse duration parameter
c       idva    I*4 - displacmeent, velocity or acceleration
c       dist    R*4 - distance
c       dt  R*4 - sample rate
c       nptin   I*4 - number of points in time series
c       tfirst  R*4 - time of first sample
c       depthr  R*4 - receiver depth
c       depths  R*4 - source depth
c       kkl I*4 - = 1 indicates Love wave eigenfunctions
c       kkr I*4 - = 1 indicates Rayleigh wave eigenfunctions
c       ostr    C*(*)   - command line string
c       mname   C*(*)   - name of model file
c       ipar    I*4 - arrau of integer command parameters from sdisp96
c       fpar    R*4 - arrau of float command parameters from sdisp96
c       dodble  L   - if .true. double length of time series for FFT
c       dolock  L   - if .true. do not use bottom layer for 
c                   first arrival computation
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green s function
c               21 - Green s function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
c
c-----

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      tp, tsv, tsh,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        real*4 tp, tsv, tsh
        real*4 sa, sc, sf, sl, sn, sr

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec

        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax

        common/receiver/l2mu,mu
        real l2mu, mu
        common/grnfmt/outfmt
        integer outfmt

c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NGRN
        parameter (NGRN=63)
        integer*4 ntau, idva,   nptin
        real*4 dt, tfirst, depths, depthr
        character ostr*(*)

        integer kerr, system

        integer NSAMP2
        parameter(NSAMP2=16384)
        common/dd/datc
        complex datc(NSAMP2)
        common/xxx/x
        real*4 x(NSAMP2)
        complex xx(NGRN)
        complex z
         integer iszrt(21)
        character*8 ostnm(2)
        character*8 ocmpnm(9)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
        logical dodble, dolock, dolat, dostep

        real svel, pvel, vsa, vsb, vsr, geom, rayp, tstar, den

        real caz, saz, c2az, s2az

        complex uz,ur,ut,duzdz,durdz,dufdz,
     1      duzdr,durdr,dufdr,duzdf,durdf,dufdf

        complex srr,srf,srz,sff,sfz,szz
        complex err,erf,erz,eff,efz,ezz
        complex del
        complex xout(9)

        character ofile*80
        integer ls, lgstr
        

        data ostnm/'MT      ', 'F0      '/
        data ocmpnm/'Uz      ','Ur      ','Ut      ',
     1       'Srr     ', 'Srf     ',
     1       'Srz     ', 'Sff     ', 'Sfz     ', 'Szz     '/


        tau = ntau * dt
        ls = lgstr(ostr)
        cpulse = ostr(1:los)
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
c-----
c       get first arrival time - NOTE for implementation of
c-----
            call frstar(dist,depths,depthr,mname,1,tp ,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            call frstar(dist,depths,depthr,mname,2,tsv,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
            call frstar(dist,depths,depthr,mname,3,tsh,
     1         pvel,svel,den,vsa,vsb,vsr,rayp,geom,tstar,dolock,.false.)
c-----
c       get the source A, C, F, L, N and density values
c-----
            sa = vsr*vsa*vsa
            sc = vsr*vsa*vsa
            sl = vsr*vsb*vsb
            sn = vsr*vsb*vsb
            sf = sa - 2.0*sl
            sr = vsr
c-----
c       lat2d96 we cannot easily predict the first arrival time
c       of a 2-D model. So we set this value to the default for
c       sac, e.g., -12345
c       JUST IGNORE PREVIOUSE - however we must preserve the sa sc cl sn sf sr
c-----
            if(dolat)then
                tp    = -12345.
                tsv   = -12345.
                tsh   = -12345.
            endif
c-----
c       Now that all spectra for a given distance are in
c       make time history and output in desired order
c-----
        np2 = nptin/2 + 1
        df = 1.0/(nptin*dt)
c-----
c       perform inverse FFT
c       If the number of data points is too large, then use supfft
c-----
c       here we will form the synthetics
c       there will be a total of 9, e.g.,
c       Uz  Ur Ut Err Ert Erz Ett Etz Ezz
c-----
c-----
c       jout   1  Uz
c              2  Ur
c              3  Ut
c              4  Err
c              5  Ert
c              6  Erz
c              7  Ett
c              8  Etz
c              9  Ezz
c----
c       This is a two pass process
c       Starting with  ur uz tz durdz and other derivatives on LUN=10
c       compute ur uz ut and strains an then output them in LUN=9
c       
c-----
        open(9,file='file9.tmp',access='sequential', 
     1   form='unformatted', status='unknown')
        rewind 9
        rewind 10
        do  k=np2,1,-1
              read(10) (xx(i),i=1,NGRN)
              call domom(xx,uz,ur,ut,duzdz,durdz,dufdz,
     1            duzdr,durdr,dufdr,duzdf,durdf,dufdf)
c-----
c     dist is in km conert to meters
c     the following is the strain formula
c     be careful of units especially the eff
c     Some dreivatives are dU/dr or dU/dz obviously
c     we must convert from km to m
c-----
              r = dist
              err = durdr
              eff = (dufdf + ur)/r
              ezz = duzdz
              erf = 0.5*(durdf/r  + dufdr - uf/r)
              erz = 0.5*(duzdr + durdz)
              efz = 0.5*(dufdz + duzdf/r)

              del = err + eff + ezz
              srz = 2*mu*erz
              sfz = 2*mu*efz
              srf = 2*mu*erf
              srr = (l2mu -mu -mu)*del + 2*mu*err
              sff = (l2mu -mu -mu)*del + 2*mu*eff
              szz = (l2mu -mu -mu)*del + 2*mu*ezz
c-----
c     CONVERT FROM CM TO M
c-----
             xout(1) = uz/1.0e+02
             xout(2) = ur/1.0e+02
             xout(3) = ut/1.0e+02
c-----
c     CONVERT FROM GM/CM_^3 and KM/S to KG/M^3 and M/S
c     err = err/1.0E+05
c     srr = srr * 1.0e+09
c-----
             xout(4) = srr*1.0e+04
             xout(5) = srf*1.0e+04
             xout(6) = srz*1.0e+04
             xout(7) = sff*1.0e+04
             xout(8) = sfz*1.0e+04
             xout(9) = szz*1.0e+04
              write(9)(xout(i),i=1,9)
         enddo
c-----
c     now make 9 passes to make the synthetics
c     for ur uz ut err ert erz ett etz ezz
c-----
         rewind 9
         do 1300 jk=1,9
                rewind 9
                do k=np2,1,-1
                    read(9) (xout(i),i=1,9)
                    datc(k)=xout(jk)
                    if(k.gt.1)then
                        datc(nptin+2-k)=conjg(datc(k))
                    endif
                enddo
c-----
c       perform inverse FFT
c-----
c       also take the opportunity to perform the correct
c       spherical spreading for an earth with radius of 6371.0 km
c-----
                if(ipar(1).gt.0)then
                    rad = distdg*3.1415927/180.0
                    fac = sqrt(dist/(6371.0*sin(rad)) )
                else
                    fac = 1.0
                endif
                    call zfour(datc,nptin,+1,dt,df)
                    do 1303 k=1,nptin
                        x(k) = real(datc(k))*fac
 1303               continue
                if(dodble)then
                    npts = nptin / 2
                else
                    npts = nptin
                endif
c-----
c           convert velocity to displacement or acceleration
c           if required. Since we are using surface waves, the
c           DC component must be zero!
c
c           note the special case for the pressure field
c-----
                IF(JK.LE.3)THEN
                    if(dostep)then
                       if(idva.eq.0)then
                           call integ(x,npts,dt,.false.)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                       endif
                    else
                       if(idva.eq.1)then
                           call deriv(x,npts,dt)
                       else if(idva.eq.2)then
                           call deriv(x,npts,dt)
                           call deriv(x,npts,dt)
                       endif
                    endif
                else
                    if(dostep)then
                        call integ(x,npts,dt,.false.)
                    endif
                endif
c-----
c           vertical
c-----
                if(jk.eq.1)then
                    cmpinc = -90.0 
                    cmpaz  =   0.0
c-----
c           radial
c-----
                else if(jk.eq.2)then
                    cmpinc = 0.0 
                    cmpaz  =   0.0
c-----
c           transverse
c-----
                else if(jk.eq.5)then
                    cmpinc = 0.0 
                    cmpaz  =  90.0
                else
                    cmpinc = -12345.
                    cmpaz  = -12345.
                endif
                stcomp = 'SYN     '
    
                cmpdt  = dt
                ksyear = 1970
                ksdoy  = 1
                ksmon  = 1
                ksday  = j
                kshour = 0
                ksmin  = 0
                ssec =  0.0

                o = 0
c-----
c
c-----
                stcomp = ocmpnm(jk)(1:3)
c-----
c     create the output file name
c-----
                 if(outfmt.eq.1)then
c                   DDDDDd_HHHh_ZZZz.cmp
                    ievdep = 10.*depths
                    istel  = 10.*depthr
                    idist  = 10.0*dist
                    write(ofile,11)idist,ievdep,istel,stcomp
   11               format(i6.6,'_',i4.4,'_',i4.4,'.',a3)
                 else if (outfmt.eq.2)then
c                   DDDDDddd_HHHhhh_ZZZzzz.cmp
                    ievdep = 1000.*depths
                    istel  = 1000.*depthr
                    idist  = 1000.*dist
                    write(ofile,12)idist,ievdep,istel,stcomp
   12               format(i7.7,'_',i6.6,'_',i6.6,'.',a3)
                 else if (outfmt.eq.3)then
c                   DDDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,13)idist,ievdep,stcomp
   13               format(i6.6,'_',i4.4,'.',a3)
                 else if (outfmt.eq.4)then
c                   DDDDdHHHh.grn
                    ievdep = 10.*depths
                    idist  = 10.0*dist
                    write(ofile,14)idist,ievdep,stcomp
   14               format(i5.5,'_',i4.4,'.',a3)
                 else if (outfmt.eq.5)then
c                   DDDdddHhhh.grn
                    ievdep = 1000.*depths
                    idist  = 1000.*dist
                    write(ofile,15)idist,ievdep,stcomp
   15               format(i6.6,'_',i4.4,'.',a3)
                 endif
                ls = lgstr(ofile)
                WRITE(0,*)'ofile:',ofile(1:ls)
                call scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
                call newhdr()
                call setfhv('DELTA   ',cmpdt,nerr)
                    bb = tfirst
                    ee = bb + (npts -1 ) * cmpdt
                    call setfhv('A       ',tp    ,nerr)
                    call setkhv('KA      ','P       ',nerr)
                    call setfhv('T0      ',tsv   ,nerr)
                    call setkhv('KT0     ','S       ',nerr)
                    call setkhv('KT1     ','S       ',nerr)
                    call setfhv('T1      ',tsh   ,nerr)
                call setfhv('B       ',bb    ,nerr)
                call setfhv('E       ',ee    ,nerr)
                call setfhv('O       ',o     ,nerr)
                call setkhv('KO      ','O       ',nerr)
                call setfhv('STLA    ',stlat ,nerr)
                call setfhv('STLO    ',stlon ,nerr)
                call setfhv('TIMMAX  ',b + indmax*cmpdt,nerr)
                call setfhv('TIMMIN  ',b + indmin*cmpdt,nerr)
                call setfhv('DEPMIN  ',depmin,nerr)
                call setfhv('DEPMAX  ',depmax,nerr)
c-----
c               our elevation is KM, SAC is set to KILOMETERS
c-----
                call setfhv('STEL    ',depthr,nerr)
                call setfhv('EVLA    ',evlat ,nerr)
                call setfhv('EVLO    ',evlon ,nerr)
c-----
c               our depth is KM, SAC is set to KILOMETERS
c-----
                call setfhv('EVDP    ',depths ,nerr)
                call setfhv('DIST    ',dist,nerr)
                call setfhv('AZ      ',evstaz,nerr)
                call setfhv('BAZ     ',stevaz,nerr)
                call setfhv('GCARC   ',distdg,nerr)
                call setfhv('CMPAZ   ',cmpaz ,nerr)
c----
c       convert from SEED to SAC component incidence
c-----
                cmpinc = cmpinc + 90.0
                call setfhv('CMPINC  ',cmpinc,nerr)
                call setfhv('DEPMEN  ',depmen,nerr)
c-----
c       set integer header value
c-----
                call setnhv('NZYEAR  ',ksyear,nerr)
                call setnhv('NZJDAY  ',ksdoy ,nerr)
                call setnhv('NZHOUR  ',kshour,nerr)
                call setnhv('NZMIN   ',ksmin ,nerr)
                isec = ssec
                smsec = (ssec - isec) * 1000.0
                call setnhv('NZSEC   ',isec  ,nerr)
                jsmsec = smsec
                call setnhv('NZMSEC  ',jsmsec,nerr)
                call setnhv('NPTS    ',npts, nerr)
c-----
c       set logical header value
c-----
                call setlhv('LEVEN   ',.true.,nerr)
                call setlhv('LPSPOL  ',.true.,nerr)
                call setlhv('LOVROK  ',.true.,nerr)
                call setlhv('LCALDA  ',.false.,nerr)
c-----
c       set character header value
c-----
                call setkhv('KSTNM    ','SYN     ' ,nerr)
                if(iobsyn.eq.2)then
                    call setkhv('KEVNM   ','SYNTHETI',nerr)
                    call setkhv('KEVNMC  ','C       ',nerr)
                else
                    call setkhv('KEVNM   ','OBSERVED',nerr)
                    call setkhv('KEVNMC  ','        ',nerr)
                endif
                call setkhv('KCMPNM   ',stcomp ,nerr)
c-----
c       set enumerated header value
c-----
                call setihv('IFTYPE   ','ITIME   ',nerr)
                call setihv('IZTYPE   ','IB      ',nerr)
                    call bwsac(12,npts,ofile(1:ls),x)
 1300   continue
        close(9,status='delete')
        return
        end
