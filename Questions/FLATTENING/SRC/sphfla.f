
        parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)
        integer NP
        parameter (NP=512)

c-----
c     LIN - unit for FORTRAN read from terminal
c     LOT - unit for FORTRAN write to terminal
c     LER - unit for FORTRAN error output to terminal
c     NL  - layers in model
c     NP  - number of unique periods
c-----
        real*4 d(NL),a(NL),b(NL),rho(NL),rtp(NL),dtp(NL),t(NP),btp(NL)
        real*4  qbinv(NL), qainv(NL)
        common/modl/ d,a,b,rho,rtp,dtp,btp


c-----
c     common for iget
c-----
        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rrho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        common/modtit/title
        character title*80

        common/para/ mmax,llw,twopi
        double precision twopi

        call getmod(2,'sph.mod',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        write(0,*)mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr
        nsph = iflsph
        write(6,*)nsph,iflsph
        nsph = 1
c-----
c     save current values
c-----
        do  i=1,mmax
            if(qb(i).gt.0.0)then
                qbinv(i) = 1.0/qb(i)
            else
                qbinv(i) =     qb(i)
            endif
            if(qa(i).gt.0.0)then
                qainv(i) = 1.0/qa(i)
            else
                qainv(i) =     qa(i)
            endif
            b(i) = vb(i)
            a(i) = va(i)
            d(i) = dl(i)
            rho(i) = rrho(i)
        enddo

c---- initialize
        if(nsph.eq.1) call sphere(0,0)
c---- set up density for Rayl
        if(nsph.eq.1) call sphere(2,1)
        iflsph = 0
c-----
c       fill up the output arrays for the new values
c-----
        do  i=1,mmax
            qb(i) = qbinv(i)
            qa(i) = qainv(i)
            vb(i) = b(i)
            va(i) = a(i)
            dl(i) = d(i)
            rrho(i) = rho(i)
        enddo


        call putmod(2,'flat.mod',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,.false.)
	end
        subroutine sphere(ifunc,iflag)
c-----
c     Transform spherical earth to flat earth
c
c     Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free
c     mode computations, in  Methods in Computational Physics, 
c         Volume 11,
c     Seismology: Surface Waves and Earth Oscillations,  
c         B. A. Bolt (ed),
c     Academic Press, New York
c
c     Love Wave Equations  44, 45 , 41 pp 112-113
c     Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
c     slowness instead of using average velocity for the layer
c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
c
c     ifunc   I*4 1 - Love Wave
c                 2 - Rayleigh Wave
c     iflag   I*4 0 - Initialize
c                 1 - Make model  for Love or Rayleigh Wave
c-----
        parameter(NL=200,NP=512)
        real*4 d(NL),a(NL),b(NL),rho(NL),rtp(NL),dtp(NL),btp(NL)
        common/modl/ d,a,b,rho,rtp,dtp,btp
        common/para/ mmax,llw,twopi
        double precision z0,z1,r0,r1,dr,ar,tmp,twopi
        save dhalf
        ar=6370.0d0
        dr=0.0d0
        r0=ar
        d(mmax)=1.0
        if(iflag.eq.0) then
            do 5 i=1,mmax
                dtp(i)=d(i)
                rtp(i)=rho(i)
    5       continue
            do 10 i=1,mmax
         aold =a(i)
         bold = b(i)
         dold = d(i)
                dr=dr+dble(d(i))
                r1=ar-dr
                z0=ar*dlog(ar/r0)
                z1=ar*dlog(ar/r1)
                d(i)=z1-z0
c-----
c               use layer midpoint
c-----
                TMP=(ar+ar)/(r0+r1)
                a(i)=a(i)*tmp
                b(i)=b(i)*tmp
                btp(i)=tmp
                r0=r1
            WRITE(0,*)i,dold,aold,bold,d(i),a(i),b(i)
   10       continue
            dhalf = d(mmax)
        else
            d(mmax) = dhalf
            do 30 i=1,mmax
                if(ifunc.eq.1)then
                     rho(i)=rtp(i)*btp(i)**(-5)
                else if(ifunc.eq.2)then
                     rho(i)=rtp(i)*btp(i)**(-2.275)
                endif
   30       continue
        endif
        d(mmax)=0.0
        return
        end
        subroutine getmod(rlun,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,listmd)
c-----
c       HISTORY
c
c       09 08 2000  gave ierr an initial default value for g77
c       01 13 2001  put in close(lun) if file is not model file
c       03 MAY 2002     Modify to permit read from standard input
c       06 JUL 2005 moved inquire to permit use of STDIN
c
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
c       Line 04: Model Units, First character 
c           is length (k for kilometer
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
c       rlun    I*4 - logical unit for reading model file. This
c                 unit is released after the use of this routine
c       mname   C*(*)   - model name - if this is stdin or 
c           STDIN just read
c                 from standard input
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
c       ierr    I*4 - 0 model file correctly read in
c               - -1 file does not exist
c               - -2 file is not a model file
c                 -3 error in the model file
c       listmd  L   - .true. list the model
c------

        implicit none
        character mname*(*), title*(*)
        integer rlun
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        integer*4 ierr
        character string*80
        logical listmd
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        logical ext
        character ftype*80
        integer lun, j, i, irefdp

c-----
c       test to see if the file exists
c-----
        ierr = 0
c-----
c       test for input
c-----
        if(MNAME(1:5).eq.'stdin' .or. mname(1:5).eq.'STDIN')then
c-----
c           do not open anything, use standard output
c-----
            lun = LIN
        else
            lun = rlun
            inquire(file=mname,exist=ext)
            if(.not.ext)then
                ierr = -1
                write(LER,*)'Model file does not exist'
                return
            endif
c-----
c           open the file
c-----
            open(lun,file=mname,status='old',form='formatted',
     1          access='sequential')
            rewind lun
        endif
c-----
c       verify the file type
c-----
c-----
c       LINE 01
c-----
        read(lun,'(a)')ftype
        if(ftype(1:5).ne.'model' .and. ftype(1:5).ne.'MODEL')then
            ierr = -2
            write(LER,*)'Model file is not in model format'
            close(lun)
            return
        endif
c-----
c       LINE 02
c-----
        read(lun,'(a)')title
c-----
c       LINE 03
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'ISO' .or. string(1:3).eq.'iso')then
            iiso = 0
        else if(string(1:3).eq.'TRA' .or. string(1:3).eq.'tra')then
            iiso = 1
        else if(string(1:3).eq.'ANI' .or. string(1:3).eq.'ani')then
            iiso = 2
        endif
c-----
c       LINE 04
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'KGS' .or. string(1:3).eq.'kgs')then
            iunit = 0
        endif
c-----
c       LINE 05
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'FLA' .or. string(1:3).eq.'fla')then
            iflsph = 0
        else if(string(1:3).eq.'SPH' .or. string(1:3).eq.'sph')then
            iflsph = 1
        endif
c-----
c       LINE 06
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'1-d' .or. string(1:3).eq.'1-D')then
            idimen = 1
        else if(string(1:3).eq.'2-d' .or. string(1:3).eq.'2-D')then
            idimen = 2
        else if(string(1:3).eq.'3-d' .or. string(1:3).eq.'3-D')then
            idimen = 3
        endif
c-----
c       LINE 07
c-----
        read(lun,'(a)')string
        if(string(1:3).eq.'CON' .or. string(1:3).eq.'con')then
            icnvel = 0
        else if(string(1:3).eq.'VAR' .or. string(1:3).eq.'var')then
            icnvel = 1
        endif
c-----
c       get lines 8 through 11
c-----
        do 900 i=8,11
            read(lun,'(a)')string
  900   continue
c-----
c       get model specifically for 1-D flat isotropic
c-----
c-----
c       get comment line
c-----
        read(lun,'(a)')string
        mmax = 0
        refdep = 0.0
        irefdp = 0
        if(iiso.eq.0)then
 1000       continue
            j = mmax +1
                read(lun,*,err=9000,end=9000)d(j),a(j),b(j),
     1              rho(j),qa(j),qb(j),etap(j),etas(j),
     2              frefp(j),frefs(j)
                if(d(j).lt.0.0)then
                    d(j) = -d(j)
                    refdep = refdep + d(j)
                    irefdp = j
                endif
            mmax = j
            go to 1000
 9000       continue
        endif
    1   format(' LAYER             H      P-VEL     S-VEL   DENSITY  ')
    2   format(' ',i5,5x,4f10.3)
    3   format(' ','-SURFACE ','- - - - - ','- - - - - ',
     1      '- - - - - ','- - - - - -')
        if(mmax.gt.0)then
            if(listmd)then
            ierr = 0
            write(LOT,1)
            do 2000 i=1,mmax
                write(LOT,2)
     1              i,d(i),a(i),b(i),rho(i)
                if(i.eq.irefdp)write(LOT,3)
 2000       continue
            endif
        else 
            ierr = -3
            write(LER,*)'Error in model file'
        endif
        if(lun.ne.LIN)close (lun)
        return
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
        real curdep, dout

        logical ext
        character cmnt*110
        cmnt(  1: 11) = '      H(KM)    '
        cmnt( 12: 22) = '   VP(KM/S)'
        cmnt( 23: 33) = '   VS(KM/S)'
        cmnt( 34: 44) = ' RHO(GM/CC)'
        cmnt( 45: 55) = '     QP    '
        cmnt( 56: 66) = '     QS    '
        cmnt( 67: 77) = '   ETAP    '
        cmnt( 78: 88) = '   ETAS    '
        cmnt( 89: 99) = '  FREFP    '
        cmnt(100:110) = '  FREFS    '

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
            if(curdep .le. refdep)then
                dout = - d(j)
            else
                dout = d(j)
            endif

            write(lun,'(4f11.4,6g11.3)')dout,a(j),b(j),
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
