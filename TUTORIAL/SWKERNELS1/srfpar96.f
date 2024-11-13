       program srfpar96
c-----
c       CHANGES
c
c-----
        integer LER, LOT, LIN
        parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
      integer nf10(NL)
c-----
c     get the control parameters from tmpsrfi.00
c-----
        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     4              sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c       iprog   = inversion type. Logical or of  rftn 2, 
c                 surf 1 - if does
c               not match, terminate this run
c       itot    = total number of iterations [37]
c       nf1     = 1 estimated stdev computed from residuals
c                 0 no scaling by residuals
c       nf2     = TOTAL number of Love Wave Gamma Modes
c                 0 DO NOT PROCESS Love Wave Gamma Data for Q
c       nf34    = TOTAL number of Love Wave modes to process 
c                 (from C and U)
c       nf5     = TOTAL number of Rayleigh Wave Gamma Modes
c                 0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
c       nf67    = TOTAL number of Rayleigh Wave modes to process 
c                 (from C and U)
c       nf10    = Input Format (from model file)
c                 0 - Inversion a, rho fixed
c                 1 - Inversion Poisson Ratio Fixed, 
c                 Rho computed from Vp
c       nfilt   = smoothing parameter 
c                 0  No model Weight  No smoothing
c                 1  Model Weight     No smoothing
c                 2  No model weight  Smoothing
c                 3  Model weight     Smoothing
c       nup     = state
c               =1 partials computed 
c               =0 surfinv run =1 before
c               =2 on update of invcsl>=1 or 3
c       dlam    = damping [32]
c       qaqb    = 2.25 default  [34]
c       wref    = reference frequency 1.0 is default [33]
c       invcsl  = 0 acausal, 1 uncoupled causel, 2 coupled causal [35]
c       invdep  = 0 last inversion was for depth
c               = 1 last inversion was for velocity and Q inverse
c       lstinv  = 2,3,4,5 depending on the last inversion
c               invdep = 1 for 2,3,4 and 0 for 5
c       twnmin  These give the receiver function window
c       twnmax
c       iter    Current iteration
c       nurftn  Number of receiver functions to be read
c       invdep  0 invert for layer thickness
c           1 invert for velocity
c       pval    0 invert for RFTN, 1 invert for SURF
c-----
      call mjsamat(wref,invcsl,invdep)
      end
      subroutine mjsamat(wref,invcsl,invdep)
c-----
c     a general purpose reformatting file
c-----
c
c     list partial derivatives or computed dispersion values
c
      parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        character*50 names
        character ostr*12
        save iwant, ihave


        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numr, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm

        character*1 lorr(2), porg(3)
        data lorr/'L','R'/,porg/'C','U','G'/
c-----
c       get earth and q model to set up causal partials
c---- 
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
c-----
c       save current values
c-----
        do 39 i=1,mmax
            if(qb(i).gt.1.0)then
                qbinv(i) = 1.0/qb(i)
            else
                qbinv(i) =     qb(i)
            endif
            if(qa(i).gt.1.0)then
                qainv(i) = 1.0/qa(i)
            else
                qainv(i) =     qa(i)
            endif
            b(i) = vb(i)
            a(i) = va(i)
            d(i) = dl(i)
            r(i) = rho(i)
   39   continue
c-----
c       get inversion control
c-----
        open(1,file='tmpsrfi.04',form='unformatted',access='sequential')
        rewind 1
        read(1) nd,m
        m2 = m + m
        read(1)(dd(i),i=1,m2)
        read(1)(wc(i),i=1,m2)
c-----
c       open dispersion file
c-----
        open(2,file='tmpsrfi.08',form='unformatted',access='sequential')
        rewind 2

        id = 16
CRBH        if(hd.eq.16 )then
CRBH            if(invdep.eq.1)then
CRBH            write(LOT,*)'Period, Partial with beta, Residual:'
CRBH            else if(invdep.eq.0)then
CRBH            write(LOT,*)'Period, Partial with depth, Residual:'
CRBH            endif
CRBH            write(LOT,*)m,'  Layers'
CRBH        endif
c-----
c       read in output of disper, derivl, derivr (unit 1) tmpsrfi.04
c       read in observed data (unit 2) tmpsrfi.08
c
c       Note there is no one-to-one correspondence in the between
c       observed and theoretical. There can be more theoretical than
c       observed because of missing periods at various modes.
c       In addition, observed data can have multiple observations
c       at each period for a given mode
c
c       The mitigating factor is that both are arranged in order
c       of increasing period, as
c           LOVE
c               PHASE
c                   MODES
c               GROUP
c                   MODES
c           RAYLEIGH
c               PHASE
c                   MODES
c               GROUP
c                   MODES
c
c-----
c
c       initialize the search
c-----
        tp1 = 0.0
        md1 = 0
        k1 = 0
        ihave = -1
   10   continue
c-----
c       get observed data
c-----
        read(1,end=40)ifn,k,md,tp,vobs,sd
C       WRITE(6,*)'***', ifn,k,md,tp,vobs,sd
        if(vobs.eq.0.0)go to 10
c-----
c       the mod(k,2) is used to for phase/group
c-----
        iwant = ifn*1000 +mod(k,2)*100 + md
   11   continue
        rewind 2
        call gttheo(ifn,k,md,tp,itst,k1,md1,tp1,vpred,iret,
     1      m,wref,invcsl,invdep,iwant,ihave)
CRBH        WRITE(6,*)'ihave,iwant,iret:',ihave,iwant,iret
        if(iret.lt.0)then
            rewind 2
            tp1 = 0.0
            md1 = 0
            k1 = 0
            ihave = -1
            go to 10
        endif
C       if(ihave.gt.iwant)then
C           rewind 2
C           go to 10
C       endif

c-----
c       if model does not support a higher mode drop from listing
c-----
        if(itst.eq.0)go to 10
c-----
c       There is a match , print out the information
c-----
c   k = 3 is gamma
CRBHN        WRITE(6,*)'---',ifn,k,md,tp,vobs,vpred,m
        call outputpar(ifn,k,md,tp,vobs,vpred,m)
    
        go to 10
      write(LOT,*)'mismatch between observed and computed values'
      write(LOT,*)'in main program - job aborted'
      stop
  40  continue
      close(1,status='keep')
      close(2,status='keep')
        return
        end

        subroutine outputpar(ilvry,icug,mode,period,vobs,vpred,m)
c-----
c       This is a general routine for listing the output
c       Note this only gives the values when there is a valid
c       comparison. For example if the data has the 1st higher mode, e.g.,
c               mode = 2 and this does not exist for the current model and
c               period, then nothing will be output
c
c       ilvry   - I  : 1 = Love
c                      2 = Rayleigh
c       icug    - I  : type of value
c                      1 = phase velocity
c                      2 = group velocity
c                      3 = gamma in exp ( -gamma r)
c       mode    - I  : mode number, e.g.,  
c                      1 = fundamental
c                      2 = first
c     period    - R  : period of observation in seconds
c       vobs    - R  : observed value
c      vpred    - R  : model predicted value
c          M    - I  : number of layers in the model
c-----
      implicit none
        integer  LOT
        parameter(LOT=6)
c-----
c     subroutine variables
c-----
      integer ilvry, icug, mode, m
      real  period, vobs,vpred
c-----
c       common block for the values of the partials for this data value
c-----
        integer NL
        parameter(NL=200)
        common/vpartial/dvda(NL),dvdb(NL),dvdh(NL),
     1             dvdqainv(NL),dvdqbinv(NL)
        real dvda, dvdb, dvdh, dvdqainv, dvdqbinv

      integer i
      real residual

CRBHN      WRITE(6,*)'+++',ilvry,icug,mode,period,vobs,vpred,m
      residual = vobs - vpred
c-----
c     tailor this output to your needs. This is a simple ASCII output
c     Note that some of the output will be zero, e.g. dcda or duda for
c     Love waves since Love waves do not depend on the P velocity
c-----
c     Each observed value that is matched by the model will
c     be output on 5 lines
c-----
c     ilvry    1 = Love, 2 = Rayleigh
c     icug     1 = phase, 2 = group,3 = gamms
c     mode     0 = fundamental, 1 = dirst etc
c     period   period in sec
c     vobs     observed value
c     vpred    model predicted value
c     residual = vobs - vpred
c     m        number of layers in the model
c-----
      WRITE(LOT,*)ilvry,icug,mode,period,vobs,vpred,residual,m
c-----
c     output the partials with 
c        dvda     - partial of value with respect to P velocity
c        dvdb     - partial of value with respect to S velocity
c        dvdh     - partial of value with respect to layer thickness
c        dvdqainc - partial of value with respect to inverse Qp 
c        dvdqbinv - partial of value with respect to inverse Qs
c-----
      WRITE(LOT,*)(dvda(i),i=1,m)
      WRITE(LOT,*)(dvdb(i),i=1,m)
      WRITE(LOT,*)(dvdh(i),i=1,m)
      WRITE(LOT,*)(dvdqainv(i),i=1,m)
      WRITE(LOT,*)(dvdqbinv(i),i=1,m)
      
        return
        end

        subroutine gttheo(ifn,k,md,tp,itst,k1,md1,tp1,c1,iret,
     1  m,wref,invcsl,invdep,iwant,ihave)
        parameter (NL=200,NL2=NL+NL)
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),ct(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc
        real*4 dcdb(NL2),dcda(NL2),dgdq(NL2),dudb(NL2),duda(NL2),
     1      dudq(NL2),dcdq(NL2),dgdv(NL2),dcdh(NL2), dudh(NL2),
     2      dgdh(NL2), dcdqb(NL2), dcdqa(NL2),dudqb(NL2),dudqa(NL2),
     3      dgdqa(NL2),dgdqb(NL2),dgdb(NL2),dgda(NL2)
        save dcdb,dcda,dgdq,dudb,duda,dudq,
     1      dcdq,dgdv,dcdh,dudh,
     2      dgdh, dcdqb, dcdqa,dudqb,dudqa,
     3      dgdqa,dgdqb,dgdb,dgda
        save cvel,gvel,gam
c-----
c       common block for the values of the partials for this data value
c-----
        common/vpartial/dvda(NL),dvdb(NL),dvdh(NL),
     1             dvdqainv(NL),dvdqbinv(NL)
        real dvda, dvdb, dvdh, dvdqainv, dvdqbinv
c-----
c       search through theoretical file to find a match
c       to a particular observation
c
c       OBSERVED
c           ifn : 1 Love 2 Rayl
c           k   : 1 Phase 2 Group
c           md  : Mode 
c           tp  : Period
c       THEORETICAL
c           itst    : Whether dispersion point found
c               : 0 mode does not exist, 1 = Love, 2 = Rayleigh
c           k1  : 1 Phase 2 Group
c           md1 : Mode
c           tp1 : Period
c           c1  : Velocity
c           dd  : Partial Derivatives
c-----
c
c       check for possibility of a double observation
c
c-----
c-----
c       initialize the dd array
c-----
        call zero(dd,1,m+m)
c-----
c       test to see if request is for information already in hand
c       NOTE: the test mod(k,2) == mod(k1,2) arises since
c       we only need dcdb, phase velocity and perehaps dcda to get
c       the anelastic attenutation coefficient gamma
c-----
        iret = 0
        dif = abs( (tp-tp1)/tp )
        jret = -1
c-----
c       test to see if information is in hand for a mode
c       not computed
c       E.g., 1 st mode not found, but we are asked for
c       phase and gamma, so read in another observed data
c       value rather than reading another theoretical value
c       and getting the two lists out of synchronization
c-----
        if(itst.eq.0 .and. md.eq.md1 .and. dif.lt.1.0e-5
     1      .and. mod(k,2).eq.mod(k1,2) )then
            jret = 1
        endif
        if(jret.eq.1)return
c-----
c       since we will not reuse previous value
c       read in new theoretical observation
c-----
   10   continue
        iret = 0
        read(2,end=30) itst,k1,md1,tp1
        ihave = 1000*itst + 100*mod(k1,2) + md1
        dif = abs ( (tp-tp1)/tp)
        if(k1.eq.2)then
            read(2,end=30)gvel,(dudb(j),j=1,m)
            read(2,end=30) (dudh(j),j=1,m)
            read(2,end=30) cvel,(dcdb(j),j=1,m)
            read(2,end=30) (dcdh(j),j=1,m)
            if(itst.eq.2)then
                 read(2,end=30)(dcda(j),j=1,m)
                 read(2,end=30)(duda(j),j=1,m)
            else
                 call zero(dcda,1,m)
                 call zero(duda,1,m)
            endif
        else
            read(2,end=30)cvel,(dcdb(j),j=1,m)
            read(2,end=30) (dcdh(j),j=1,m)
            if(itst.eq.2)then
                 read(2,end=30)(dcda(j),j=1,m)
            else
                 call zero(dcda,1,m)
            endif
            call zero(duda,1,m)
            call zero(dudb,1,m)
            call zero(dudh,1,m)
        endif
        if(mod(k,2).ne.mod(k1,2).or.md.ne.md1.or.dif.gt.1.0e-5)then
            go to 10
        endif
c-----
c       define correct parameter
c-----

        call getgam(dcdb,b,qbinv,dcda,a,qainv,m,itst,cvel,tp1,
     1      gam)
        call zero(dcdq,1,m)
        call zero(dudq,1,m)
        call zero(dgdv,1,m)
CRBH        if(invcsl.ge.1)then
        call cslmod(gvel,cvel,gam,wref,tp1,m,
     1      qbinv,qainv,dcdb,dcda,dudb,duda,itst,a,b,2,
     1      dcdh,dudh,dgdh,
     1      dcdqb,dcdqa,dudqb,dudqa,dgdb,dgda,dgdqb,dgdqa)
CRBH        endif
c-----
c       setup causal, uncoupled if invcsl = 1
c-----
C       if(invcsl.eq.1 .or. invdep.eq.0)then
C           call zero(dudq,1,m)
C           call zero(dcdq,1,m)
C           call zero(dgdv,1,m)
C       endif
c-----
c       copy the values to the generic dcdv, dudv etc
c-----
        if(k.eq.1)then
c----
c       phase velocity
c-----
             c1 = cvel
             call cpy(dcda ,dvda    ,0,m)
             call cpy(dcdb ,dvdb    ,0,m)
             call cpy(dcdh ,dvdh    ,0,m)
             call cpy(dcdqa,dvdqainv,0,m)
             call cpy(dcdqb,dvdqbinv,0,m)
        else if(k.eq.2)then
c-----
c       group velocity
c-----
             c1 = gvel
             call cpy(duda ,dvda    ,0,m)
             call cpy(dudb ,dvdb    ,0,m)
             call cpy(dudh ,dvdh    ,0,m)
             call cpy(dudqa,dvdqainv,0,m)
             call cpy(dudqb,dvdqbinv,0,m)
        else if(k.eq.3)then
c-----
c       gamma 
c-----
             c1 = gam
             call cpy(dgda ,dvda    ,0,m)
             call cpy(dgdb ,dvdb    ,0,m)
             call cpy(dgdh ,dvdh    ,0,m)
             call cpy(dgdqa,dvdqainv,0,m)
             call cpy(dgdqb,dvdqbinv,0,m)
        endif
        if(ihave.ne.iwant)go to 10
        return
   30   continue
            iret = -1
        return
        end


        subroutine zero(x,m,n)
        real*4 x(*)
c-----
c       set x(m) = 0, x(m+1)=0, ..., x(n)=0
c-----
        do 100 i=m,n
            x(i) = 0.0
  100   continue
        return
        end

        subroutine cpy(x,y,m,n)
        real*4 x(*),y(*)
c-----
c       copy contents of array x into array y with offset m
c-----
        do 100 i=1,n
            y(i+m) = x(i)
  100   continue
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

        subroutine gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,invdep,pval,
     3      sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c     07 DEC 2012  - place the refdef in
c          the rdum5 position of tmpsrfi.00 so that
c          refdef is preserved in all model files
c-----
      common/depref/refdep

c-----
c       read control file
c-----
        integer NL
        parameter (NL=200)
        integer nf10(NL)
        open(1,file='tmpsrfi.00',form='unformatted',access='sequential')
        rewind 1
        read(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,dlam
     1      ,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,
     2      nurftn,invdep,pval,
     3      sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
        refdep = rdum5
        close(1,status='keep')
        return
        end

        subroutine cslmod(u,c,gam,wref,tp,m,
     1      qbinv,qainv,dcdb,dcda,dudb,duda,ilvry,a,b,k1,
     1      dcdh,dudh,dgdh,
     1     dcdqb,dcdqa,dudqb,dudqa,dgdb,dgda,dgdqb,dgdqa)
        integer NL,NLAY, NL2
        parameter (NL=200,NLAY=200,NL2=NL+NL)
        real u,c,gam,wref,tp
        integer m,k1,ilvry
        real*4 a(NL),b(NL)
        real*4 qbinv(NL),qainv(NL)
        real*4 dcdh(NL2), dudh(NL2), dgdh(NL2)
        real*4 dcdb(NL2),dcda(NL2),dudb(NL2),duda(NL2)
        real*4 dcdqb(NL2),dcdqa(NL2),dudqb(NL2),dudqa(NL2),dgdb(NL2),
     1         dgda(NL2),dgdqb(NL2),dgdqa(NL2)


c-----
c       modify phase, group velocities and partials for
c       constant Q causality
c
c       ilvry 1 = L, 2 = R
c           k1  : 1 Phase 2 Group
c           tp : Period
c-----
        c0 = c
        u0 = u
        f  = 1./tp
        pi = 3.1415927
        omega = 6.2831853*f
        faclog = alog(f/wref)/pi
        facg = 2.*gam*u0/(pi*omega)
c-----
c       get correct phase velocity
c-----
        c = c0 + 2.*gam*c0*c0*faclog/omega
c-----
c       get correct group velocity
c-----
        cmc0c0 = (c-c0)/c0
        if(k1.eq.2)then
            u0c0 = u0/c0
            u = u0 *( 1. + (2. - u0c0)*cmc0c0 + facg ) 
            uu0 = u/u0
        endif
        faca = 0.0
        fars = 0.0
        
c-----
c       get correct partials
c-----
        do 100 i=1,m
c-----
c          save values that will be overwritten
c          dcdh is not changed
c          dudh is not changed
c-----
            dcdp = dcda(i)
            dcds = dcdb(i)
            dcdb(i) = dcds*(1. + faclog*qbinv(i) )
            dcda(i) = dcdp*(1. + faclog*qainv(i) )
            dcdqa(i) = faclog*dcdp*a(i)
            dcdqb(i) = faclog*dcds*b(i)
            dgdqa(i) = 0.5*omega*dcdp*a(i)/(c0*c0)
            dgdqb(i) = 0.5*omega*dcds*b(i)/(c0*c0)
            if(k1.eq.2)then
                 duds = dudb(i)
                 dudp = duda(i)
                 dudb(i) = duds*(uu0 -u0c0*cmc0c0 + facg )
     1            + dcds*u0c0*(-2.*facg +u0c0*qbinv(i)/pi +
     2               (2.-u0c0)*(faclog*qbinv(i) -cmc0c0) +
     3               u0c0 *cmc0c0 )
                 duda(i) = dudp*(uu0 -u0c0*cmc0c0 + facg )
     1            + dcdp*u0c0*(-2.*facg +u0c0*qainv(i)/pi +
     2               (2.-u0c0)*(faclog*qainv(i) -cmc0c0) +
     3               u0c0 *cmc0c0 )
                 dudqa(i) = 
     1                   u0c0*(2.-u0c0)*dcdqa(i) 
     1                   +
     1                u0c0*u0c0*dcdp*a(i)/pi
                 dudqb(i) =               
     1                   u0c0*(2.-u0c0)*dcdqb(i) 
     1                  +
     1                u0c0*u0c0*dcds*b(i)/pi
            endif
  100   continue
        return
        end

        subroutine getgam(dcdb,b,qbinv,dcda,a,qainv,m,ilvry,cc,tp,
     1          gam)
        real*4 dcdb(*),b(*),dcda(*),a(*),qbinv(*), qainv(*)
c-----
c       determine spatial attenuation at constant frequency
c       and partial derivatives
c-----
        omega = 6.2831853/tp
        factr = 0.5*omega/(cc*cc)
        sumg = 0.0
        do 100 i=1,m
            sumg = sumg + dcdb(i)*b(i)*qbinv(i)
            if(ilvry.eq.2)then
                sumg = sumg + dcda(i)*a(i)*qainv(i)
            endif
  100   continue
        gam = factr * sumg
        return
        end
