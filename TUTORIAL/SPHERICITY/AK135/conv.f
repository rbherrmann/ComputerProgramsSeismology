c-----
c       convert minos layer cake into a model96 format
c       note minos assumes a cubic spline interpolation
c       to interpolate between values
c-----=
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/depref/refdep
        real refdep

        real depth(NL), density(NL), pvelocity(NL), svelocity(NL)
        real qkappa(NL), qmu(NL)

        character instr*80

        character mname*80
        character title*80
        title = 'AK135-F model converted to layers'
        ibunit = 0
        iiso = 0
        iflsph = 1
        idimen = 1
        icnvel = 0
        refdep = 0.0

c-----
c        skip 4 lines
c-----
        do i=1,4
                read(5,'(a)')instr
        enddo
        i = 0
 1000        continue
        read(5,*,end=2000)dep,dens,pvel,svel,qk,qm
        i = i + 1
        inmax = i
        depth(i) = dep
        density(i) = dens
        pvelocity(i) = pvel
        svelocity(i) = svel
        qkappa(i) = qk
        qmu(i) = qm
        go to 1000
 2000        continue
c-----
c        output for minos
c-----
        open(1,file='tak135-f.txt',status='unknown',form=
     1                'formatted',access='sequential')
        rewind 1
        write(1,'(a,a)')' MODEL: tak135-f ','
     1    http://wwwrses.anu.edu.au/seismology/ak135/ak135f.html'
c-----
c       ifanis: 1 for an anisotropic (transversely isotropic) model
c       tref  : reference perdion (seconds) for causal attenuation
c       ifdeck: 1 model is given in tabular form
c-----
        ifanis = 1
        tref = 1.0
        ifdeck = 1
        write(1,'(i5,f10.2,i5)')ifanis,tref,ifdeck
c-----
c       BEWARD THIS IS HARD WIRED
c       nic: index of solid side of the inner cord boundary (ICB)
c       noc: index of the fluid side of the mantle core boundary (MCB)
c-----
        nic=25
        noc=71
        write(1,'(3i5)')inmax, nic, noc
c-----
c      output the isotropic model in TI format
c      r rho vpv vsv qkappa qshear vph vsh eta
c       units are meters fro length, m/s for velocity, kg/m^3 for density
c-----
        do i=inmax,1,-1
             rhom = density(i)*1000.0
             rm = (6371.0-depth(i))*1000.0
             vpvm = pvelocity(i)*1000.0
             vsvm = svelocity(i)*1000.0
             vphm = pvelocity(i)*1000.0
             vshm = svelocity(i)*1000.0
             etam = 1.0
        write(1,105)rm,rhom,vpvm,vsvm,qkappa(i),
     1       qmu(i),vphm,vshm,etam
  105 format(f8.0, 3f9.2,2f9.1,2f9.2, f9.5)
        enddo
        close (1)
c-----
c        now create the layers
c-----
        mmax = 0.0
        do i=1,inmax -1 
                dd = depth(i+1) - depth(i)
c-----
c                minos uses two values at the same depth to
c                indicate a discontinuity
c-----
                if(dd.gt.0.0)then
c-----
c        make the velocity in each layer such that the
c        vertical travel time is the same - this means
c        that we assume that there is a linear gradient
c        for simplicity we assume that the slowness is
c        linear - thus by the trapezoidal rule
c
c        1/vavg = 0.5 (1/v1 + 1/v2)
c or
c        vavg = 2 v1 v2/(v1 + v2)
c do the same for t*
c for density assume the same mass in the column
c
c        Ben-Menahem and Singh 10.203
c
c        QP^-1 = LQ^-1 + (1-L)QK^-1
c
c        where L = (4/3)(VS/VP)^2

c        Dahlen and Tromp 9.59 9.60
c
c----- 
                        mmax = mmax + 1
                        d(mmax) = depth(i+1) - depth(i)
                        vavg = 0.5*(1./pvelocity(i) + 1./pvelocity(i+1))
                        a(mmax) = 1./vavg
                        if(svelocity(i) .eq.0.0)then
                                b(mmax) = 0.0
                        else
                        vavg = 0.5*(1./svelocity(i) + 1./svelocity(i+1))
                        b(mmax) = 1./vavg
                        endif
                        fac = (4./3.)*(b(mmax)/a(mmax))**2
                        rho(mmax) = 0.5*(density(i)+density(i+1))
                        if(qmu(i).eq.0.0)then
                                qb(mmax) = 0.0
                                qsavginv = 0.0
                        else
                                qsavginv = 0.5*(1./qmu(i) + 1./qmu(i+1))
                                qb(mmax) = 1./qsavginv
                        endif
                        qpavginv = 0.5*(1./qkappa(i) + 1./qkappa(i+1))
                        qainv = fac*qsavginv + (1.0-fac)*qpavginv
                        qa(mmax) = 1.0/qainv
                        etap(mmax) = 0.0
                        etas(mmax) = 0.0
                        frefp(mmax) = 1.0
                        frefs(mmax) = 1.0
                endif

                
        enddo
c-----
c       do not output all layers since the CPS code cannot handle a
c       solid - fluid - solid model, so drop the inner core
c-----
        do i=1,mmax-1
          if(b(i+1).gt.0.0 .and. b(i).eq.0.0)then
             nmax = i
          endif
        enddo
        call putmod(1,'tak135sph.mod',nmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
        end

        subroutine putmod(wlun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,lverby)
c-----
c       CHANGES
c       03 MAY 2002 permit write to standard output
c-----
c       General purpose model input
c       This model specification is designed to be as 
c           general as possible
c
c       Input lines
c       Line 01: MODEL
c       Line 02: Model Name
c       Line 03: ISOTROPIC or ANISOTROPIC or 
c           TRANSVERSELY ANISOTROPIC
c       Line 04: Model Units, First character is length (k for kilometer
c           second is mass (g for gm/cc), third is time (s for time)
c       Line 05: FLAT EARTH or SPHERICAL EARTH
c       Line 06: 1-D, 2-D or 3-D
c       Line 07: CONSTANT VELOCITY
c       Line 08: open for future use
c       Line 09: open for future use
c       Line 10: open for future use
c       Line 11: open for future use
c       Lines 12-end:   These are specific to the model
c           For ISOTROPIC the entries are
c           Layer Thickness, P-velocity, S-velocity, Density, Qp, Qs,
c           Eta-P, Eta S (Eta is frequency dependence), 
c           FreqRefP, FreqRefP
c-----
cMODEL
cTEST MODEL.01
cISOTROPIC
cKGS
cFLAT EARTH
c1-D
cCONSTANT VELOCITY
cLINE08
cLINE09
cLINE10
cLINE11
c H  VP  VS   RHO   QP  QS   ETAP   ETAS REFP  REFS
c1.0    5.0 3.0 2.5 0.0 0.0 0.0 0.0 1.0 1.0
c2.0    5.1 3.1 2.6 0.0 0.0 0.0 0.0 1.0 1.0
c7.0    6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
c10.0   6.5 3.8 2.9 0.0 0.0 0.0 0.0 1.0 1.0
c20.0   7.0 4.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0
c40.0   8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0
c-----
c-----
c       wlun    I*4 - logical unit for writing model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       title   C*(*)   - title of the model file
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       idimen  I*4 - 1 1-D
c               - 2 2-D
c               - 3 3-D
c       icnvel  I*4 - 0 constant velocity
c                 1 variable velocity
c       lverby  L   - .false. quiet output
c------
        
        implicit none

        character mname*(*)
        character title*(*)
        integer wlun
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        logical lverby
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer lun
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/depref/refdep
        real refdep

        integer lgstr, lt
        integer j
        real curdep

        logical ext
        character cmnt*110
        cmnt(  1: 11) = '      H(KM)    '
        cmnt( 12: 22) = '   VP(KM/S)'
        cmnt( 23: 33) = '   VS(KM/S)'
        cmnt( 34: 44) = ' RHO(GM/CC)'
        cmnt( 45: 55) = '         QP'
        cmnt( 56: 66) = '         QS'
        cmnt( 67: 77) = '       ETAP'
        cmnt( 78: 88) = '       ETAS'
        cmnt( 89: 99) = '      FREFP'
        cmnt(100:110) = '      FREFS'

        lt = lgstr(title)
c-----
c       test to see if the file exists
c-----
        if(MNAME(1:6).eq.'stdout' .or. mname(1:6).eq.'STDOUT')then
c-----
c           do not open anything, use standard output
c-----
            lun = LOT
        else
            inquire(file=mname,exist=ext)
            if(ext .and.  lverby)then
                write(LER,*)'Overwriting Existing model File'
            endif
            lun = wlun
c-----
c           open the file
c-----
            open(lun,file=mname,status='unknown',form='formatted',
     1          access='sequential')
            rewind lun
        endif
c-----
c       verify the file type
c-----
c-----
c       LINE 01
c-----
        write(lun,'(a)')'MODEL.01'
c-----
c       LINE 02
c-----
        write(lun,'(a)')title(1:lt)
c-----
c       LINE 03
c-----
        if(iiso.eq.0)then
            write(lun,'(a)')'ISOTROPIC'
        else if(iiso.eq.1)then
            write(lun,'(a)')'TRANSVERSE ANISOTROPIC'
        else if(iiso.eq.2)then
            write(lun,'(a)')'ANISOTROPIC'
        endif
c-----
c       LINE 04
c-----
        write(lun,'(a)')'KGS'
c-----
c       LINE 05
c-----
        if(iflsph.eq.0)then
            write(lun,'(a)')'FLAT EARTH'
        else if(iflsph.eq.1)then
            write(lun,'(a)')'SPHERICAL EARTH'
        endif
c-----
c       LINE 06
c-----
        if(idimen.eq.1)then
            write(lun,'(a)')'1-D'
        else if(idimen.eq.2)then
            write(lun,'(a)')'2-D'
        else if(idimen.eq.3)then
            write(lun,'(a)')'3-D'
        endif
c-----
c       LINE 07
c-----
        if(icnvel.eq.0)then
            write(lun,'(a)')'CONSTANT VELOCITY'
        else if(icnvel.eq.1)then
            write(lun,'(a)')'VARIABLE VELOCITY'
        endif
c-----
c       put lines 8 through 11
c-----
        write(lun,'(a)')'LINE08'
        write(lun,'(a)')'LINE09'
        write(lun,'(a)')'LINE10'
        write(lun,'(a)')'LINE11'
c-----
c       put model specifically for 1-D flat isotropic
c-----
c-----
c       put comment line
c-----
        write(lun,'(a)')cmnt(1:110)
        curdep = 0.0
        
        do 1000 j=1,mmax
            curdep = curdep + abs(d(j))
            if(curdep .le. refdep)d(j) = - d(j)
            write(lun,'(4f11.4,6g11.3)')d(j),a(j),b(j),
     1          rho(j),qa(j),qb(j),etap(j),etas(j),
     2          frefp(j),frefs(j)
 1000   continue
        if(lun.ne.LOT)close (lun)
        return
        end
        function lgstr(str)
c-----
c       function to find the length of a string
c       this will only be used with file system path names
c       thus the first blank 
c       indicates the end of the string
c-----
        implicit none
        character*(*) str
        integer*4 lgstr
        integer n, i
        n = len(str)
        lgstr = 1
        do 1000 i=n,1,-1
            lgstr = i
            if(str(i:i).ne.' ')goto 100
 1000   continue
  100   continue
        return
        end
