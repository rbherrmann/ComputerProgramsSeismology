        program saccompress
c-----
c       This program undisperses a signal given the pahse velocity
c       and outputs as a sacfile with the suffix .cmp
c
c       
c-----
        
c-----
c       waveform
c-----
        integer NMAX
        parameter (NMAX=16384)
        real data(NMAX)
        complex zdata(NMAX), zzdat(NMAX)
        character fname*80
        integer nerr, npts
        real dist, dt
        real btime, otime

c-----
c       dispersion
c-----
c-----
c       read a SURF96 dispersion file
c-----
c       fildsp  Ch* name of dispersion file
c       lun I*4 logical unit for read
c       ierr    I*4 Error code: -1 file does not exist
c       ndsp    I*4 number of dispersion values
c
c       ftype   A   should be SURF96
c       jlorr() I*4 1 = Love, 2 = Rayleigh
c       jobs()  I*4 1 = phase velocity, 2=group velocity, 3=gamma
c       jobsyn() I*4    1 = observed, 2 synthetic
c       jmode() I*4 mode: 0=Fund, 1=1 st, etc
c       fper()  R*4 period
c       fobs()  R*4 Array of ground velocities
c       fobserr()   R*4 Error in u()
c-----
        integer NOBS
        parameter(NOBS=1000)
        integer lun, ierr, ndsp
        integer*4 jlorr(NOBS), jobs(NOBS), jobsyn(NOBS), jmode(NOBS)
        real*4  fper(NOBS)   , fobs(NOBS)   , fobserr(NOBS)
        character dispfile*80


        real fval(NOBS), cval(NOBS)
        real df, t0, fac
        integer ls, n2pts, iu
        complex zfac
        character instr*80
        character*80 dogeom
c-----
c       function
c-----
        integer n2, n21
     
c-----
c       initially hardwire
c-----
C        fname = 'IUANTO_HHR__.SAC.cut'
C       read(5,'(a)')fname
C       dispfile='C.dsp'
        call mgtarg(1,fname)
        call mgtarg(2,instr)
        read(instr,'(bn,f15.0)')dist
        call mgtarg(3,dispfile)
        call mgtarg(4,dogeom)
        fmin = 0.01
        fmax = 0.05
        call rddisp(dispfile,2,nerr,ndsp,NOBS,
     1      jlorr,jobs,jobsyn,jmode,
     2      fper,fobs,fobserr)


        ls = lgstr(fname)
        call brsac (1,NMAX,fname,data,nerr)
        call getfhv('DELTA   ',dt,nerr)
        call getfhv('B       ',btime,nerr)
        call getfhv('O       ',otime,nerr)
        call getnhv('NPTS    ',npts,nerr)
        call bwsac(2,NMAX,fname(1:ls)//'.cmp',data)
c-----
c       travel time of the first sample is t0
c-----
        t0 = btime - otime
c-----
c       fill the time series and perhaps taper
c-----
        n2 = npts
        call npow2(n2)
        do i =1,n2
           if(i.le.npts)then
              zdata(i) = cmplx(data(i),0.0)
           else
              zdata(i) = cmplx(0.0,0.0)
           endif
        enddo
        call zfour(zdata,n2,-1,dt,df) 

c-----
c       the model of the observation in the frequency domain is that
c       H(f) = A(f) exp(- i k r )     where k is the frequency dependent wavenumebr
c       To get a synthetic to start at a t0, the model is
c       H(f) = A(f) exp ( i [omega t0 - k r] )
c       Thus if the time series starts at t0, and a fourier transform is taken,then
c       the resultang Fourier transform will be  H(f) = A(f) exp ( i [omega t0 - k r] )
c
c       To compress the signal, we form
c       exp( - i [omega t0 +  k r ]) H(f)
c
c       finally for stacking and analysis we should center the resuling compressed signal. This is
c       done in the transfrom domain by creating a new spectra which is H(00, -H(1), H(2),-H93) ...
c
c-----
c       compress here
c-----
        n21 = n2/2 + 1

        flip = 1.0
        do i=1,n21
             freq = (i-1)*df
             call cinterp(fper,fobs,ndsp,freq,c)
             if(c.lt.0.0)c= 1000000.
             fac = 6.2831853*freq*(-t0 + dist/c)
             zfac=cmplx(cos(fac),sin(fac))
             zzdat(i) = zdata(i)*zfac
             call taper(freq,fmin,fmax,ftaper)
             zzdat(i) = zzdat(i)*ftaper
             if(i.gt.1)then
                 zzdat(n+2-i) = conjg(zzdat(i))
             endif
             zzdat(i) = flip*zzdat(i)
             flip = - flip
         enddo
c-----
c       inverse fft
c-----
        call zfour(zzdat,n2,+1,dt,df) 
c-----
c       output the result
c-----
        call setnhv('NPTS',n2,nerr)
        call setfhv('B       ',-(n2/2 -1)*dt,nerr)
        call setfhv('O       ',0.0          ,nerr)
        write(6,*)dogeom
        do i =1,n2
           data(i) = real(zzdat(i))
          if(dogeom.eq.'TRUE')then
               data(i) = data(i)*dist
          endif
        enddo
        call bwsac(2,NMAX,fname(1:ls)//'.cmp',data)
        end

        subroutine taper(freq,flow,fhigh,ftaper)
c-----
c       10% cubic taper
c-----
        df = (fhigh-flow)
        f1 = flow
        f2 = flow  + 0.1*df
        f3 = fhigh - 0.1*df
        f4 = fhigh
        if(freq.lt.f1)then
           ftaper = 0.0
        else if(freq.ge.f1.and.freq.lt.f2)then
            p = (freq-f1)/(f2-f1)
            p = 2*p -1
            ftaper = 0.5+0.75*p*(1.0-p*p/3.0)
        else if(freq.ge.f2.and.freq.lt.f3)then
           ftaper = 1.0
        else if(freq.ge.f3.and.freq.lt.f4)then
            p = (freq-f4)/(f2-f4)
            p = 2*p -1
            ftaper = 0.5+0.75*p*(1.0-p*p/3.0)
        else if(freq.gt.f4)then
           ftaper = 0.0
        endif
        return 
        end

        subroutine cinterp(fper,fobs,nper,freq,c)
        real fper(*), fobs(*)
        integer nper
        real freq, c

        integer i
        real per, pval
        logical reverse
c-----
c       get the order for tghe table
c-----
        if(fper(2).gt.fper(1))then
            reverse = .false.
        else
            reverse = .true.
        endif
       
c-----
c       interpolate - assume that the periods are ordered
c-----
        per = 1.0/freq
        c = -1
        do i=1,nper -1
           if(reverse)then
c               12 11 10
              if(per.le.fper(i) .and. per.ge.fper(i+1))then
                 pval = (per-fper(i))/(fper(i+1)-fper(i))
                 c = pval*fobs(i) + ( 1.0-pval)*fobs(i+1)
              endif            
           else
c               10 11 12
              if(per.ge.fper(i) .and. per.le.fper(i+1))then
                 pval = (per-fper(i))/(fper(i+1)-fper(i))
                 c = pval*fobs(i) + ( 1.0-pval)*fobs(i+1)
              endif            
   
           endif
        enddo
        return
        end

        subroutine rddisp(fildsp,lun,ierr,ndsp,NOBS,
     1      jlorr,jobs,jobsyn,jmode,
     2      fper,fobs,fobserr)
c-----
c       CHANGES
c       19 JAN 2002 - Make obsyn=2 for synthetic, 1 for observed
c                     to be compatible with f96subs.f 2=syn/1=obs
c       05 MAR 2002 - synthetic denoted by 'T' and observed by 'X'
c-----
c       read a SURF96 dispersion file
c-----
c       fildsp  Ch* name of dispersion file
c       lun I*4 logical unit for read
c       ierr    I*4 Error code: -1 file does not exist
c       ndsp    I*4 number of dispersion values
c
c       ftype   A   should be SURF96
c       jlorr() I*4 1 = Love, 2 = Rayleigh
c       jobs()  I*4 1 = phase velocity, 2=group velocity, 3=gamma
c       jobsyn() I*4    1 = observed, 2 synthetic
c       jmode() I*4 mode: 0=Fund, 1=1 st, etc
c       fper()  R*4 period
c       fobs()  R*4 Array of ground velocities
c       fobserr()   R*4 Error in u()
c       
c
c
CSURF96 R U T 0  350.000     2.694     0.635  
c-----
c-----
c       command line variables
c-----
        character fildsp*(*)
        integer lun, ierr, ndsp,  NOBS
        integer*4 jlorr(NOBS), jobs(NOBS), jobsyn(NOBS), jmode(NOBS)  
        real*4  fper(NOBS)   , fobs(NOBS)   , fobserr(NOBS)
c-----
c       internal variables
c-----
        logical ext
        character instr*132
        character ic*1
        integer ilorr, imode, iobs
        ierr = -1
        ext = .false.
        inquire(file=fildsp,exist=ext)
        if(.not. ext)then
            ierr = -1
        else
            open(lun,file=fildsp,access='sequential',
     1      form='formatted',status='old')
            rewind lun
c-----
c           read values and parse them correctly
c-----
            ndsp = 0
 1000       continue
            read(lun,'(a)',end=9999)instr
            ls = lgstr(instr)
c-----
c           do the parsing
c-----
            if(instr(1:6).eq.'surf96' .or.
     1          instr(1:6).eq.'SURF96')then
c-----
c               now get wave type
c-----      
                lsep = 6
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'R'. or. ic(1:1).eq.'r')then
                    ilorr = 2
                else if(ic(1:1).eq.'L'. or. ic(1:1).eq.'l')then
                    ilorr = 1
                else if(ic(1:1).eq.'A'. or. ic(1:1).eq.'a')then
                    ilorr = 0
                else
                    go to 1000
                endif
c-----
c               now get observation type
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'C'. or. ic(1:1).eq.'c')then
                    iobs = 1
                else if(ic(1:1).eq.'U'. or. ic(1:1).eq.'u')then
                    iobs = 2
                else if(ic(1:1).eq.'G'. or. ic(1:1).eq.'g')then
                    iobs = 3
                else
                    go to 1000
                endif
c-----
c               now get whether observed or synthetic
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'T'. or. ic(1:1).eq.'t')then
                    iobsyn = 2
c-----
c       the F designation is deprecated
c-----
                else if(ic(1:1).eq.'F'. or. ic(1:1).eq.'f')then
                    iobsyn = 1
                else if(ic(1:1).eq.'X'. or. ic(1:1).eq.'x')then
                    iobsyn = 1
                else
                    go to 1000
                endif
c-----
c-----
c               now get the values using list directed IO
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                read(instr(lnobl:ls),*)
     1          imode,per,obs,obserr
c-----
c       test for array dimensions
c-----
                if(ndsp.eq.NOBS)go to 9999
                ndsp = ndsp + 1
                jlorr(ndsp) = ilorr
                jobs(ndsp)  = iobs
                jobsyn(ndsp)  = iobsyn
                jmode(ndsp) = imode
                fper(ndsp)  = per
                fobs(ndsp)  = obs
                fobserr(ndsp)  = obserr
            endif
            
            go to 1000
        endif
 9999       continue
            close (lun)
            ierr = 0
        return
        end


        subroutine getblnk(instr,lsep,ls,lnobl)
c-----
c       determine first non-blank character
c
c       instr   Ch* Character string to be parsed
c       lsep    I*4 index of last non blank character
c       ls  I*4 length of input string
c       lnobl   I*4 index of first non blank character
c-----
        character instr*(*)
        integer lsep,ls,lnobl
        character tab*1
        tab=char(9)
        lnobl = lsep+1
        igotit = 0
        do 1000 i=lsep+1,ls
            if(igotit.eq.0)then
            if(instr(i:i).ne.' ' .and. instr(i:i).ne.tab)then
                lnobl = i
                igotit = 1
            endif
            endif
 1000   continue
        return
        end
c-----
c       CHANGES
c
c       10 APR 2003 - argument to etoh corrected
c       01 APR 2004 - forced wsac1 to set depmax, depmin depmin
c       07 JAN 2005 - forces bwsac and bwsac2 to 
c           delete the file before writing. This is
c               required if the file already exists and 
c           is larger than necessary
c               Check the close(UNIT,status='delete') across compilers
c       07 JUL 2007 - change nrec to nbytes in brsac to emphasize the
c           number of bytes to read. Note that that open(...,recl=N) is
c           compiler dependent, sometimes N=bytes and other times words,
c           whatever a word is
c       24 JAN 2009 - corrected rsac1 so prevent overflow - it returns the
c           MIN(number of points in the file, buffer size)
c       10 MAR 2009 - added  (Gusst Nolet) correctly set the 
c              setihv - set enumerated header value for IZTYPE as
c               IUNKN (Unknown)     - integer 5
c               IB (Begin time)     - integer 9
c               IDAY (Midnight of refernece GMT day) - integer 10
c               IO (Event origin time) - integer 11
c               IA (First arrival time) - integer 12
c               ITn (User defined time pick n, n=0,9) - integer 13-22
c       11 JUN 2009 - hermann.zeyen@u-psud.fr found the following problem
c               in subroutine ARSAC: if the number of data is multiple 
c               of five, the program tries to read one line too much.
c               This is not a problem for sacsubc.c
c       26 APR 2016 - fixed error in brsach in which the character header was not read
c-----
        subroutine brsac (IRU,LN,name,data,nerr)
c-----
c       IRU I*4 logical unit for IO
c       LN  I*4 length of data array
c       name    C*  Name of file to be opened
c       rhdr    R*4 Real header
c       ihdr    I*4 Integer Header
c       chdr    C*  Character Header
c       data    R*4 Array of trace values
c       nerr    I*4 -1 file does not exist
c               -2 data points in file exceed dimension
c
c       NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c       LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
        real*4 data(LN)

        logical ext

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*(*) name
c-----
c  Read real and integer header blocks to find actual number
c  of waveform data points which is stored in ihdr(10).
c
c-----
        inquire(file=name,exist=ext)
        if(.not. ext)then
            ihdr(10) = 0
            nerr = -1
            return
        endif
            nerr = 0
            open (IRU,file=name,form='unformatted',
     &          access='direct',recl=440,status='old')
                read (IRU,rec=1) (rhdr(i), i=1,70),
     &                   (ihdr(i), i=1,40)
            close (IRU)
c-----
c
c  Read header and waveform data blocks using recored 
c           length of 158*4=632.
c
c-----
            if(ihdr(10).gt.LN)then
                maxpts = LN
                ihdr(10) = LN
                nerr = -2
            else 
                maxpts = ihdr(10)
                nerr = 0
            endif
c-----
            nbytes=632+4*maxpts
            nread = 0
c-----
c       because of SUNOS Fortran problems with IO transfers 
c       more than 2048 bytes, read these  chunks in 
c----- 
            ndat = maxpts
            if(nbytes.gt.2048)then
                open (IRU,file=name,form='unformatted',
     &              access='direct',recl=2048)
                ndat1 = (2048 - 632) / 4
                irec = 1
                read (IRU,rec=irec,err=1001) (rhdr(i), i=1,70),
     &                   (ihdr(i), i=1,40),
     &                   (chdr(i), i=1,24),
     &                   (data(i), i=1,ndat1)
                nread = nread + ndat1
 1000           continue
                nl = nread + 1
                nh = nl + 512 - 1
                if(nh.gt.ndat)then
                    nh = ndat
                endif
                if(nl.gt.ndat)go to 1001
                irec = irec + 1
                read (IRU,rec=irec,err=1001) (data(i), i=nl,nh)
                nread = nread + (nh-nl+1)

                go to 1000
 1001           continue
            close (IRU)
            else
                open (IRU,file=name,form='unformatted',
     &              access='direct',recl=nbytes)
                read (IRU,rec=1) (rhdr(i), i=1,70),
     &                   (ihdr(i), i=1,40),
     &                   (chdr(i), i=1,24),
     &                   (data(i), i=1,ndat)
            close (IRU)
            endif
            if(ihdr(10).gt.LN)then
                maxpts = LN
                ihdr(10) = LN
            else 
                maxpts = ihdr(10)
            endif
            ihdr(10) = maxpts
        return
        end

        subroutine brsach(IRU,name,nerr)
c-----
c       IRU I*4 logical unit for IO
c       name    C*  Name of file to be opened
c       rhdr    R*4 Real header
c       ihdr    I*4 Integer Header
c       chdr    C*  Character Header
c       nerr    I*4 -1 file does not exist
c               -2 data points in file exceed dimension
c
c       NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c       LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----

        logical ext

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*(*) name
c-----
c  Read real and integer header blocks to find actual number
c  of waveform data points which is stored in ihdr(10).
c
c-----
        inquire(file=name,exist=ext)
        if(.not. ext)then
            ihdr(10) = 0
            nerr = -1
            return
        endif
            nerr = 0
            open (IRU,file=name,form='unformatted',
     &          access='direct',recl=632,status='old')
                read (IRU,rec=1) (rhdr(i), i=1,70),
     &                   (ihdr(i), i=1,40),
     &                   (chdr(i), i=1,24)
            close (IRU)
        return
        end

        subroutine arsac (IRU,LN,name,data,nerr)
c-----
c       IRU I*4 logical unit for IO
c       LN  I*4 length of data array
c       name    C*  Name of file to be opened
c       rhdr    R*4 Real header
c       ihdr    I*4 Integer Header
c       chdr    C*  Character Header
c       data    R*4 Array of trace values
c       nerr    I*4 -1 file does not exist
c               -2 data points in file exceed dimension
c
c       NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c       LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  This routine reads files written in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
        logical ext
        real*4 data(LN)
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*(*) name
c-----
        inquire(file=name,exist=ext)
        if(.not. ext)then
            ihdr(10) = 0
            nerr = -1
            return
        endif
        nerr = 0
            open (IRU,file=name,status='old',access='sequential')
            rewind IRU
c----- 
c  Read real header block.
c-----
            j1=1
            j2=5
            do 1110  i=1,14
                read (IRU,'(5g15.0)') (rhdr(j), j=j1,j2)
                j1=j1+5
                j2=j2+5
 1110       continue
c-----
c  Read integer header block.
c-----
            j1=1
            j2=5
            do 1120 i=1,8
                read (IRU,'(5i10)') (ihdr(j), j=j1,j2)
                j1=j1+5
                j2=j2+5
 1120       continue
c-----
c  Read character header block.
c-----
            j1=1
            j2=3
            do 1130 i=1,8
                read (IRU,'(3a8)') (chdr(j), j=j1,j2)
                j1=j1+3
                j2=j2+3
 1130       continue
            if(ihdr(10).gt.LN)then
                maxpts = LN
                ihdr(10) = LN
                nerr = -2
            else 
                maxpts = ihdr(10)
                nerr = 0
            endif
c-----
c  Read waveform data organized in a five columns block.
c-----
c      FIX by Hermann Zeyen 11 JUN 2009
c           nrow=(maxpts/5)+1
c-----
            nrow=maxpts/5
            if(nrow*5 .lt. maxpts) nrow=nrow+1
            do 1140 i=1,nrow
                l=i*5
                k=l-4
                read (IRU,'(5g15.0)') (data(j), j=k,l)
 1140       continue
            close (IRU)
            if(ihdr(10).gt.LN)then
                maxpts = LN
                ihdr(10) = LN
            else 
                maxpts = ihdr(10)
            endif
        return
        end

        subroutine arsach(IRU,name,nerr)
c-----
c       IRU I*4 logical unit for IO
c       name    C*  Name of file to be opened
c       rhdr    R*4 Real header
c       ihdr    I*4 Integer Header
c       chdr    C*  Character Header
c       nerr    I*4 -1 file does not exist
c               -2 data points in file exceed dimension
c
c       NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c       LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  This routine reads files written in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----
        logical ext
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*(*) name
c-----
        inquire(file=name,exist=ext)
        if(.not. ext)then
            ihdr(10) = 0
            nerr = -1
            return
        endif
        nerr = 0
            open (IRU,file=name,status='old',access='sequential')
            rewind IRU
c-----
c  Read real header block.
c-----
            j1=1
            j2=5
            do 1110  i=1,14
                read (IRU,'(5g15.0)') (rhdr(j), j=j1,j2)
                j1=j1+5
                j2=j2+5
 1110       continue
c-----
c  Read integer header block.
c-----
            j1=1
            j2=5
            do 1120 i=1,8
                read (IRU,'(5i10)') (ihdr(j), j=j1,j2)
                j1=j1+5
                j2=j2+5
 1120       continue
c-----
c  Read character header block.
c-----
            j1=1
            j2=3
            do 1130 i=1,8
                read (IRU,'(3a8)') (chdr(j), j=j1,j2)
                j1=j1+3
                j2=j2+3
 1130       continue
            close (IRU)
        return
        end

        subroutine getfhv(strcmd,fval,nerr)
c-----
c       Get float header value
c
c       strcmd  C*8 String to key on
c       val R*4 Real value returned
c       nerr    I*4 Error condition
c                   0 no error
c                   1336 Value not defined
c                   1337 Header variable does not exist
c-----
        character strcmd*(*)
        logical streql
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output real header
c-----
        fval = -12345.
        nerr = -1
        do 1000 i=1,70
            if(streql(strcmd,rstr(i))) then
                nerr = 0
                fval = rhdr(i)
            endif
 1000   continue
        return
        end

        subroutine getnhv(strcmd,ival,nerr)
c-----
c       Get integer header value
c
c       str C*8 String to key on
c       ival    R*4 integer value returned
c       nerr    I*4 Error condition
c                   0 no error
c                   1336 Value not defined
c                   1337 Header variable does not exist
c-----
        character strcmd*(*)
        logical streql
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output integer header
c-----
        ival = -12345
        nerr = -1
        do 2000 i=1,40
            if(streql(strcmd,istr(i))) then
                nerr = 0
                ival = ihdr(i)
            endif
 2000   continue
        return
        end

        subroutine getkhv(strcmd,cval,nerr)
c-----
c       Get character header value
c
c       strcmd  C*8 String to key on
c       cval    C*8 character value returned
c       nerr    I*4 Error condition
c                   0  no error
c                   1336 Value not defined
c                   1337 Header variable does not exist
c-----
        character strcmd*(*)
        logical streql

        character cval*8
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output character header
c-----
        nerr = -1
        do 3000 i=1,24
            if(streql(strcmd,cstr(i))) then
                nerr = 0
                cval = chdr(i)
c-----
c               safety 14 AUG 2007 - get rid of non printing
c-----
                do j=1,8
                     if(ichar(cval(j:j)).lt.31)then
                         cval(j:j) = ' '
                     endif
                enddo

            endif
 3000   continue
        return
        end

        subroutine getlhv(strcmd,lval,nerr)
        character strcmd*(*)
        logical lval
        
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       input logical header
c-----
        call getnhv(strcmd,ival,nerr)
        if(ival.eq.0)then
            lval = .false.
        else
            lval = .true.
        endif
        return
        end

c---------------------------------------------------------
        subroutine bwsac (IWU,LN,name,data)
c---------------------------------------------------------

c
c  This routine writes out a waveform data in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character name*(*)
        real*4 data(LN)
c-----
c       remove the original file so that the output length is
c       never greater than desired. Else the dregs of the 
c           first will remain
c-----
            open (IWU,file=name,form='unformatted',
     &          access='sequential',status='unknown')
            rewind IWU
            close (IWU,status='delete')
c
c  The actual number of waveform data points is stored in integer
c  header 10. The file recored length is 158*4=632.
c
            nrec=632+4*ihdr(10)
            open (IWU,file=name,form='unformatted',
     &          access='direct',recl=nrec,status='unknown')
            write (IWU,rec=1) (rhdr(i),i=1,70),
     &                (ihdr(k),k=1,40),
     &                (chdr(j),j=1,24),
     &                (data(l), l=1,ihdr(10))
            close (IWU)
        return
        end
c---------------------------------------------------------
        subroutine awsac (IWU,LN,name,data)
c---------------------------------------------------------
c
c  This routine writes out files in SAC alpha (ASCII) format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
        real*4 data(LN)
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character name*(*)
c
            open (IWU,file=name,status='unknown',
     &               access='sequential')
            rewind IWU
c
c  Write real header block.
c
            j1=1
            j2=5
            do 1100 i=1,14
                write (IWU,'(5g15.7)') (rhdr(j), j=j1,j2)
                j1=j1+5
                j2=j2+5
 1100       continue
c
c  Write integer header block.
c
            j1=1
            j2=5
            do 1110 i=1,8
                write (IWU,'(5i10)') (ihdr(j), j=j1,j2)
                j1=j1+5
                j2=j2+5
 1110       continue
c
c  Write character header block.
c
            j1=1
            j2=3
            do 1120 i=1,8
                write (IWU,'(3a8)') (chdr(j), j=j1,j2)
                j1=j1+3
                j2=j2+3
 1120       continue
c
c  Ensure the last row is padded with zeros, if actual number of
c  waveform points is less than the product of number of rows
c  and number of columns which constitutes the data block.
c
            nrow=(ihdr(10)/5)+1
            nrc5=nrow*5
            if (nrc5 .gt. ihdr(10)) then
                nrcx=ihdr(10)+1
                do 1140 i=nrcx,nrc5
                    data(i)=0.0
 1140           continue
            end if
c
c  Write waveform data in five columns format.
c
            do 1150 i=1,nrow
                k=i*5
                j=k-4
                write (IWU,'(5g15.7)') (data(l), l=j,k)
 1150       continue
            close (IWU)
        return
        end

        subroutine setfhv(strcmd,fval,nerr)
c-----
c       Set float header value
c
c       strcmd  C*8 String to key on
c       fval    C*8 real value set
c       nerr    I*4 Error condition
c                   0  no error
c                   1337 Header variable does not exist
c-----
        character strcmd*(*)
        logical streql
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output real header
c-----
        do 1000 i=1,70
            if(streql(strcmd,rstr(i))) rhdr(i) = fval
 1000   continue
        return
        end

        subroutine setnhv(strcmd,ival,nerr)
c-----
c       Set integer header value
c
c       strcmd  C*8 String to key on
c       ival    C*8 integer value set
c       nerr    I*4 Error condition
c                   0  no error
c                   1337 Header variable does not exist
c-----
        character strcmd*(*)
        integer ival, nerr

        logical streql
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output integer header
c-----
        do 2000 i=1,40
            if(streql(strcmd,istr(i))) ihdr(i) = ival
 2000   continue
        return
        end

        subroutine setkhv(strcmd,cval,nerr)
c-----
c       Set character header value
c
c       strcmd  C*8 String to key on
c       cval    C*8 character value set
c       nerr    I*4 Error condition
c                   0  no error
c                   1337 Header variable does not exist
c-----
        character strcmd*(*)
        character cval*8
        logical streql
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output character header
c-----
        do 3000 i=1,24
            if(streql(strcmd,cstr(i))) chdr(i) = cval
 3000   continue
        return
        end


        subroutine setlhv(strcmd,lval,nerr)
        character strcmd*(*)
        logical lval
        
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        character*8 rstr(70), istr(40), cstr(24)
        data (rstr(i),i=1,45)/
     1  'DELTA   ', 'DEPMIN  ', 'DEPMAX  ', 'SCALE   ', 'ODELTA  ', 
     1  'B       ', 'E       ', 'O       ', 'A       ', 'FMT     ', 
     1  'T0      ', 'T1      ', 'T2      ', 'T3      ', 'T4      ', 
     1  'T5      ', 'T6      ', 'T7      ', 'T8      ', 'T9      ', 
     1  'F       ', 'RESP0   ', 'RESP1   ', 'RESP2   ', 'RESP3   ', 
     1  'RESP4   ', 'RESP5   ', 'RESP6   ', 'RESP7   ', 'RESP8   ', 
     1  'RESP9   ', 'STLA    ', 'STLO    ', 'STEL    ', 'STDP    ', 
     1  'EVLA    ', 'EVLO    ', 'EVEL    ', 'EVDP    ', 'MAG     ', 
     1  'USER0   ', 'USER1   ', 'USER2   ', 'USER3   ', 'USER4   ' /
        data (rstr(i),i=46,70)/
     1  'USER5   ', 'USER6   ', 'USER7   ', 'USER8   ', 'USER9   ', 
     1  'DIST    ', 'AZ      ', 'BAZ     ', 'GCARC   ', 'SB      ', 
     1  'SDELTA  ', 'DEPMEN  ', 'CMPAZ   ', 'CMPINC  ', 'XMINIMUM', 
     1  'XMAXIMUM', 'YMINIMUM', 'YMAXIMUM', 'ADJTM   ', 'TIMMAX  ', 
     1  'TIMMIN  ', 'FHDR67  ', 'FHDR68  ', 'FHDR69  ', 'FHDR70  ' /
        data istr/
     1  'NZYEAR  ', 'NZJDAY  ', 'NZHOUR  ', 'NZMIN   ', 'NZSEC   ', 
     1  'NZMSEC  ', 'NVHDR   ', 'NINF    ', 'NHST    ', 'NPTS    ', 
     1  'NSNPTS  ', 'NSN     ', 'NXSIZE  ', 'NYSIZE  ', 'NHDR15  ', 
     1  'IFTYPE  ', 'IDEP    ', 'IZTYPE  ', 'IHDR4   ', 'IINST   ', 
     1  'ISTREG  ', 'IEVREG  ', 'IEVTYP  ', 'IQUAL   ', 'ISYNTH  ', 
     1  'IHDR11  ', 'IHDR12  ', 'IHDR13  ', 'IHDR14  ', 'IHDR15  ', 
     1  'IHDR16  ', 'IHDR17  ', 'IHDR18  ', 'IHDR19  ', 'IHDR20  ', 
     1  'LEVEN   ', 'LPSPOL  ', 'LOVROK  ', 'LCALDA  ', 'LHDR5   ' /
        data cstr/
     1  'KSTNM   ', 'KEVNM   ', 'KEVNMC  ', 'KHOLE   ', 
     1  'KO      ', 'KA      ', 'KT0     ', 'KT1     ', 
     1  'KT2     ', 'KT3     ', 'KT4     ', 'KT5     ', 
     1  'KT6     ', 'KT7     ', 'KT8     ', 'KT9     ', 
     1  'KF      ', 'KUSER0  ', 'KUSER1  ', 'KUSER2  ', 
     1  'KCMPNM  ', 'KNETWK  ', 'KDATRD  ', 'KINST   '
     1  /
c-----
c       output logical header
c-----
        if(lval)then
            call setnhv(strcmd,1,nerr)
        else
            call setnhv(strcmd,0,nerr)
        endif
        return
        end


        subroutine newhdr()
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        call inihdr()
c       ITIME
        ihdr(16) = 1
c       LEVEN = TRUE
        ihdr(36) = 1
c       LOVROK = TRUE
        ihdr(38) = 1
c       LCALDA = TRUE
        ihdr(39) = 1
c       IZTYPE = IUNKN
        ihdr(18) = 5

        return
        end

        subroutine inihdr()
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
        do 1100 i=1,70
            rhdr(i)= -12345.0
 1100   continue
        do 1110 i=1,35
            ihdr(i)= -12345
 1110   continue
        ihdr(7)=6
        ihdr(8)=0
        ihdr(9)=0
        do 1120 i=36,40
            ihdr(i)=0
 1120   continue
        do 1130 i=1,24
            chdr(i)='-12345  '
 1130   continue
        return
        end


        logical function streql(str1,str2)
        character str1*(*), str2*(*)
        character nstr1*8, nstr2*8
c-----
c       determine if two strings are equal
c-----
        nstr1 = ' '
        l1 = lgstr(str1)
        nstr1(1:l1) = str1(1:l1)
        nstr2 = ' '
        l2 = lgstr(str2)
        nstr2(1:l2) = str2(1:l2)
c-----
        if(nstr1 .eq. nstr2)then
            streql = .true.
        else
            streql = .false.
        endif
        return 
        end

        subroutine getihv(strcmd,strval,nerr)
c-----
c       Get enumerated header value
c
c       strcmd  C*8 String to key on
c       strval  C*8 real value set
c       nerr    I*4 Error condition
c               0  no error
c               1336 Header variable undefined
c               1337 Header variable does not exist
c-----
        character strcmd*(*), strval*8
        parameter (NEVAL=50)
        character*8 eval(NEVAL)
c-----
c       header integer equivalents of enumerated values
c       e.g., IDISP == 2
c-----
      data eval/'ITIME   ','IRLIM   ','IAMPH   ','IXY     ','IUNKN   ', 
     1  'IDISP   ', 'IVEL    ', 'IACC    ', 'IB      ', 'IDAY    ', 
     2  'IO      ', 'IA      ', 'IT0     ', 'IT1     ', 'IT2     ', 
     3  'IT3     ', 'IT4     ', 'IT5     ', 'IT6     ', 'IT7     ', 
     4  'IT8     ', 'IT9     ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5  'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6  'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7  'ISRO    ', 'INUCL   ', 'IPREN   ', 'IPOSTN  ', 'IQUAKE  ', 
     8  'IPREQ   ', 'IPOSTQ  ', 'ICHEM   ', 'IOTHER  ', 'IGOOD   ', 
     9  'IGLCH   ', 'IDROP   ', 'ILOWSN  ', 'IRLDTA  ', 'IVOLTS  '/
c-----
            call getnhv(strcmd,nval,nerr)
            strval = '        '
            if(nerr.eq.0)then
                if(nval.ge.1 .and. nval.le.NEVAL)then
                    strval = eval(nval)
                endif
            endif
        return
        end

        subroutine setihv(strcmd,strval,nerr)
c-----
c       Set enumerated header value
c
c       strcmd  C*8 String to key on
c       strval  C*8 real value set
c       nerr    I*4 Error condition
c               0  no error
c               1336 Header variable undefined
c               1337 Header variable does not exist
c-----
        character strcmd*(*), strval*8
        parameter (NEVAL=50)
        character*8 eval(NEVAL)
        character*8 ival(NEVAL)
        logical streql
c-----
c       header integer equivalents of enumerated values
c       e.g., IDISP == 2
c-----
      data eval/'ITIME   ','IRLIM   ','IAMPH   ','IXY     ','IUNKN   ', 
     1  'IDISP   ', 'IVEL    ', 'IACC    ', 'IB      ', 'IDAY    ', 
     2  'IO      ', 'IA      ', 'IT0     ', 'IT1     ', 'IT2     ', 
     3  'IT3     ', 'IT4     ', 'IT5     ', 'IT6     ', 'IT7     ', 
     4  'IT8     ', 'IT9     ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5  'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6  'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7  'ISRO    ', 'INUCL   ', 'IPREN   ', 'IPOSTN  ', 'IQUAKE  ', 
     8  'IPREQ   ', 'IPOSTQ  ', 'ICHEM   ', 'IOTHER  ', 'IGOOD   ', 
     9  'IGLCH   ', 'IDROP   ', 'ILOWSN  ', 'IRLDTA  ', 'IVOLTS  '/
c-----
c       equivalence of header field position and enumerated value
c       e.g., IDEP can be IUNKN, IDISP, IVEL, IVOLTS or IACC
c-----
      data ival/'IFTYPE  ','IFTYPE  ','IFTYPE  ','IFTYPE  ','IDEP    ', 
     1  'IDEP    ', 'IDEP    ', 'IDEP    ', 'IZTYPE  ', 'IZTYPE  ', 
     2  'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 
     3  'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 'IZTYPE  ', 
     4  'IZTYPE  ', 'IZTYPE  ', 'IRADNV  ', 'ITANNV  ', 'IRADEV  ', 
     5  'ITANEV  ', 'INORTH  ', 'IEAST   ', 'IHORZA  ', 'IDOWN   ', 
     6  'IUP     ', 'ILLLBB  ', 'IWWSN1  ', 'IWWSN2  ', 'IHGLP   ', 
     7  'ISRO    ', 'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 
     8  'IEVTYP  ', 'IEVTYP  ', 'IEVTYP  ', 'IQUAL   ', 'IQUAL   ', 
     9  'IQUAL   ', 'IQUAL   ', 'IQUAL   ', 'ISYNTH  ', 'IDEP    '/

c-----
c       now do the work, parse the table for the match
c           strcmd = ival and strval = eval, then 
c           using the table index, I, 
c               do a call setnhv(strcmd,I,nerr)
c
c       However, the IUNKN is used in both IDEP and IZTYPE 
c       and IOTHER is used in both IEVTYP and IQUAL
c       
c-----
        nerr = 0
        if(streql(strcmd,'IDEP    ')
     1          .and. streql(strval,'IUNKN   '))then
            call setnhv('IDEP    ',5,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IUNKN   '))then
            call setnhv('IZTYPE  ',5,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IB      '))then
            call setnhv('IZTYPE  ',9,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IDAT    '))then
            call setnhv('IZTYPE  ',10,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IO      '))then
            call setnhv('IZTYPE  ',11,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IA      '))then
            call setnhv('IZTYPE  ',12,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT0     '))then
            call setnhv('IZTYPE  ',13,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT1     '))then
            call setnhv('IZTYPE  ',14,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT2     '))then
            call setnhv('IZTYPE  ',15,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT3     '))then
            call setnhv('IZTYPE  ',16,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT4     '))then
            call setnhv('IZTYPE  ',17,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT5     '))then
            call setnhv('IZTYPE  ',18,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT6     '))then
            call setnhv('IZTYPE  ',19,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT7     '))then
            call setnhv('IZTYPE  ',20,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT8     '))then
            call setnhv('IZTYPE  ',21,nerr)
        else if(streql(strcmd,'IZTYPE  ')
     1          .and. streql(strval,'IT9     '))then
            call setnhv('IZTYPE  ',22,nerr)
        else if(streql(strcmd,'IEVTYP  ')
     1          .and. streql(strval,'IOTHER  '))then
            call setnhv('IEVTYP  ',44,nerr)
        else if(streql(strcmd,'IQUAL   ')
     1          .and. streql(strval,'IOTHER  '))then
            call setnhv('IQUAL   ',44,nerr)
        else
            nerr = 1336
c-----
c       IFTYPE
c-----
            do 100 i=1,NEVAL
                if(
     1              streql(strval,eval(i)))then
                    call setnhv(strcmd,i,nerr)
                endif
  100   continue
        endif
        return
        end
            
        subroutine rsac1(infile,y,npts,btime,dt,maxpts,nerr)
        character infile*(*)
        real y(maxpts)
        integer npts
        real btime
        real dt
        integer maxpts
        integer nerr
c-----
c       PURPOSE READ AN EVENLY SPACED SAC FILE
c
c       read a binary sac file with evenly sampled data
c       infile  Char    name of file
c       y   R   array of values
c       npts    I   number of points in data
c       btime   R   start time
c       dt  R   sample interval
c       maxpts  I   maximum number of points to read
c       nerr    I   error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

        npts = maxpts
c-----
c       rad up to maxpts, indicate the number actually read
c-----
        call brsac(1,maxpts,infile,y,nerr)
        call getnhv('NPTS    ',npts,ierr)
        call getfhv('DELTA   ',dt  ,ierr)
        call getfhv('B       ',btime  ,ierr)
        return
        end

        subroutine wsac0(ofile,x,y,nerr)
        character ofile*(*)
        real x(*)
        real y(*)
        integer nerr

        logical leven
c-----
c       WRITE  SAC FILE USING CURRENT HEADER VALUES
c       however we should look at the header value LEVEN
c       to decide if we should write out at x,y pairs or as
c       y values only
c
c       RBH 2000 08 31 look at the header value,
c       then do a wsac1 or a wsac2
c
c       write a binary sac file with evenly sampled data
c
c       ofile   Char    name of file
c       ydummy  R   array of values
c       y   R   array of values
c       nerr    I   error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

        npts = ihdr(10)
c-----
c       this is wrong
c       'To write a SAC file to disk using current header values'
C       well I am forcing even spaced data
c           had to do this to kludge owens SAC
c-----
        call getlhv('LEVEN   ',leven,nerr)
c-----
c       the value of npts is set in the brsac if the actual number
c       of points exceeds the dimension limit
c-----
        call getnhv('NPTS    ',npts,nerr)
        call getfhv('DELTA   ',dt,nerr)
        call getfhv('B       ',b ,nerr)
        if(leven)then
            call wsac1(ofile,y,npts,b,dt,nerr)
        else
            call wsac2(ofile,x,npts,y,nerr)
        endif
C       call bwsac(1,npts,ofile,y)
        nerr = 0
        return
        end

        subroutine wsac1(ofile,y,npts,btime,dt,nerr)
        character ofile*(*)
        real y(npts)
        integer npts
        real btime
        real dt
        integer nerr 
c----- 
c       PURPOSE: WRITE AN EVENLY SPACED SAC FILE
c
c       write a binary sac file with evenly sampled data
c       ofile   Char    name of file
c       y   R   array of values
c       npts    I   number of points in data
c       btime   R   start time
c       dt  R   sample interval
c       maxpts  I   maximum number of points to read
c       nerr    I   error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

        integer indmax, indmin
c-----
c       create a new header
c-----
C       call newhdr()
c-----
c       place essential values into the new header
c-----
        call scmxmn(y,npts,depmax,depmin,depmen,indmax,indmin)
        call setfhv('DEPMAX', depmax, ierr)
        call setfhv('DEPMIN', depmin, ierr)
        call setfhv('DEPMEN', depmen, ierr)
        call setnhv('NPTS    ',npts,nerr)
        call setfhv('DELTA   ',dt  ,nerr)
        call setfhv('B       ',btime  ,nerr)
        call setfhv('TIMMAX  ',btime + indmax*dt  ,nerr)
        call setfhv('TIMMIN  ',btime + indmin*dt  ,nerr)
        call setihv('IFTYPE  ','ITIME   ',nerr)
        e = btime + (npts -1 )*dt
        call setfhv('E       ',e     ,nerr)
        call setlhv('LEVEN   ',.true.,nerr)
        call setlhv('LOVROK  ',.true.,nerr)
        call setlhv('LCALDA  ',.true.,nerr)
        call bwsac(1,npts,ofile,y)
        nerr = 0
        return
        end

        subroutine wsac2(ofile,x,npts,y,nerr)
        character ofile*(*)
        real y(npts), x(npts)
        integer npts
        integer nerr 
c----- 
c       PURPOSE: WRITE AN UNEVENLY SPACED SAC FILE
c
c       write a binary sac file with evenly sampled data
c       ofile   Char    name of file
c       y   R   array of values
c       npts    I   number of points in data
c       x   R   array of independent variable
c       nerr    I   error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c       create a new header
c-----
        call newhdr()
c-----
c       place essential values into the new header
c-----
        call setnhv('NPTS    ',npts,nerr)
        call setihv('IFTYPE  ','ITIME   ',nerr)
        call setfhv('B       ',y(1)  ,nerr)
        call setfhv('E       ',y(npts)  ,nerr)
        call bwsac2(1,npts,ofile,x,y,npts)
        nerr = 0
        return
        end

        subroutine rsac2(ofile,x,npts,y,maxpts,nerr)
        character ofile*(*)
        real y(npts), x(npts)
        integer npts
        integer nerr 
        integer maxpts
c----- 
c       PURPOSE: WRITE AN UNEVENLY SPACED SAC FILE
c
c       write a binary sac file with evenly sampled data
c       ofile   Char    name of file
c       y   R   array of values
c       npts    I   number of points in data
c       x   R   array of independent variable
c       nerr    I   error return
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c       create a new header
c-----
        
c-----
c       place essential values into the new header
c-----
        call brsac2(1,maxpts,ofile,x,y,npts)
        nerr = 0
        return
        end

c---------------------------------------------------------
        subroutine bwsac2 (IWU,LN,name,x,y,npts)
c---------------------------------------------------------

c
c  This routine writes out a waveform data in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

        character name*(*)
        real*4 x(LN)
        real*4 y(LN)
c-----
c       remove the original file so that the output length is
c       never greater than desired. Else the dregs of 
c           the first will remain
c-----
            open (IWU,file=name,form='unformatted',
     &          access='sequential',status='unknown')
            rewind IWU
            close (IWU,status='delete')
c
c  The actual number of waveform data points is stored in integer
c  header 10. The file recored length is 158*4=632.
c
            nrec=632+8*npts
            open (IWU,file=name,form='unformatted',
     &          access='direct',recl=nrec,status='unknown')
            write (IWU,rec=1) (rhdr(i),i=1,70),
     &                (ihdr(k),k=1,40),
     &                (chdr(j),j=1,24),
     &                (x(l), l=1,npts),
     &                (y(l), l=1,npts)
            close (IWU)
        return
        end
c---------------------------------------------------------
        subroutine brsac2 (IRU,LN,name,x,y,npts)
c---------------------------------------------------------
c-----
c       IRU I*4 logical unit for IO
c       LN  I*4 length of data array
c       name    C*  Name of file to be opened
c       rhdr    R*4 Real header
c       ihdr    I*4 Integer Header
c       chdr    C*  Character Header
c       x   R*4 Array of x values
c       y   R*4 Array of y values
c       nerr    I*4 -1 file does not exist
c               -2 data points in file exceed dimension
c
c       NOTE IF THE BINARY FILE HAS MORE THAN LN POINTS, THEN ONLY
c       LN POINTS ARE USED
c-----
c
c  This routine reads waveform data written in SAC binary format.
c
c  Written by Hafidh A. A. Ghalib, 1988.
c
c-----

        logical ext

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)

        character*(*) name
        integer LN
        real*4 x(LN), y(LN)
        integer npts

        
c-----
c  Read real and integer header blocks to find actual number
c  of waveform data points which is stored in ihdr(10).
c
c-----
        inquire(file=name,exist=ext)
        if(.not. ext)then
            ihdr(10) = 0
            nerr = -1
            return
        endif
            nerr = 0
            open (IRU,file=name,form='unformatted',
     &          access='direct',recl=440,status='old')
                read (IRU,rec=1) (rhdr(i), i=1,70),
     &                   (ihdr(i), i=1,40)
            close (IRU)
c-----
c
c  Read header and waveform data blocks using recored 
c           length of 158*4=632.
c
c-----
            if(ihdr(10).gt.LN)then
                maxpts = LN
                ihdr(10) = LN
                nerr = -2
            else 
                maxpts = ihdr(10)
                nerr = 0
            endif
            npts = maxpts
            nrec=632+8*npts
            open (IRU,file=name,form='unformatted',
     &          access='direct',recl=nrec,status='unknown')
            rewind IRU
            read (IRU,rec=1) (rhdr(i),i=1,70),
     &                (ihdr(k),k=1,40),
     &                (chdr(j),j=1,24),
     &                (x(l), l=1,npts),
     &                (y(l), l=1,npts)
        return
        end

C       character kkdate*20
C       character kktime*20
C       ncdate = 20
C       call kadate(1998,133,ncdate,kkdate,nerr)
C       write(6,*)kkdate
C       call kadate(1998, 23,ncdate,kkdate,nerr)
C       write(6,*)kkdate
C       call kadate(1998,  3,ncdate,kkdate,nerr)
C       write(6,*)kkdate
C      call katime(23,33,43,890,12,kktime,nerr)
C       write(6,*)kktime
C      call katime( 3, 3, 3, 90,12,kktime,nerr)
C       write(6,*)kktime
C      call katime( 3, 3, 3,  0,12,kktime,nerr)
C       write(6,*)kktime
C      end

      subroutine katime(ihour,imin,isec,imsec,nctime,kktime,nerr)
*=====================================================================
* PURPOSE: To convert four integer fields representing hour, minute,
*          second, and millisecond to an alphanumeric equivalent
*          of the form:  HH:MM:SS.SSS
*=====================================================================
* INPUT ARGUMENTS:
*    IHOUR:   Integer hour field.
*    IMIN:    Integer minute field.
*    ISEC:    Integer second field.
*    IMSEC:   Integer millisecond field.
*    NCTIME:  Maximum length of output alphanumeric field.
*             NCTIME must be at least 12.
*=====================================================================
* OUTPUT ARGUMENTS:
*    KKTIME:  Equivalent alphanumeric field.
*             Set to all asterisks if an error occurs.
*    NERR:    Error flag. Set to 0 if no error occurred.
*             Potential error numbers: 0905, 0907.
*===================================================================
        integer ihour, imin, isec, imsec,nctime,nerr
        character kktime*(*)
        kktime(1:nctime) = ' '
        nerr = 0
        sec = isec + 0.001*imsec
        write(kktime,'(i2,1x,i2,1x,f6.3)')ihour, imin, sec
        kktime(3:3) = ':'
        kktime(6:6) = ':'
        do 1000 i=1,12
            if(kktime(i:i).eq.' ')kktime(i:i) = '0'
 1000   continue
        return
        end

      subroutine kadate(iyear,ijday,ncdate,kkdate,nerr)
*=====================================================================
* PURPOSE: To convert two integer fields representing year and
*          julian day to an alphanumeric equivalent
*          of the form:  MMM DD (JJJ), YYYY
*=====================================================================
* INPUT ARGUMENTS:
*    IYEAR:   Integer year field.
*    IJDAY:   Integer day field.
*    NCDATE:  Maximum length of output alphanumeric field.
*             NCDATE must be at least 18.
*=====================================================================
* OUTPUT ARGUMENTS:
*    KKDATE:  Equivalent alphanumeric date field.
*             Set to all asterisks if an error occurred.
*    NERR:    Error flag. Set to 0 if no error occurred.
*             Potential error numbers:
*=====================================================================
        integer iyear, ijday, nerr
        character kkdate*(*)
        character kmonth*4
       dimension kmonth(12)
       data kmonth/'JAN ','FEB ','MAR ','APR ','MAY ','JUN ',
     1            'JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/


        nerr = 0
        call mnthdy(iyear,ijday,imonth,iday)
        kkdate(1:ncdate) = ' '
        kkdate(1:4)=kmonth(imonth)
        if(iday.gt.10)then
            write(kkdate(5:6),'(i2)')iday
        else
            kkdate(5:5) = '0'
            write(kkdate(6:6),'(i1)')iday
        endif
        kkdate(7:14)=' (   ), '
        if(ijday.gt.99)then
            write(kkdate(9:11),'(i3)')ijday
        else if(ijday.le.99 .and .ijday.ge.10)then
            kkdate(9:9) = '0'
            write(kkdate(10:11),'(i2)')ijday
        else
            kkdate(9:10)='00'
            write(kkdate(11:11),'(i1)')ijday
        endif
        
        write(kkdate(15:18),'(i4)')iyear

        return
        end

        subroutine etoh(epoch,date,str,doy,
     1      year,month,day,hour,minute,second)
        implicit none
c-----
c       convert from epoch time to human time
c
c       epoch   - R*8 time in seconds relative to 0000 1 Jan 1970
c       date    - I*4 Julian date
c       str - C*  printable string C*32
c-----
        real*8 epoch
        integer*4 date
        character str*(*)
        integer*4 diy, doy, hour, minute, year, month
        integer*4 day
        real*4 second
        real*8 seclft
        integer i
        integer leapdy
        

        str=' '
        doy = (epoch/86400.0d+00)
        seclft = dmod(epoch, 86400.0d+00)
        hour = 0
        minute = 0
        second = 0.00
c-----
c       compute hours minutes seconds
c-----
        if(seclft .ne. 0.00d+00)then
c-----
c                   before 1970 subtract and add a day
c-----
            if(seclft .lt. 0.0d+00)then
                doy = doy - 1
                seclft = seclft + 86400.00d+00
            endif
            hour = (seclft/3600.00d+00)
            seclft = dmod(seclft,3600.0d+00)
            minute = seclft/60.0d+00
            second = dmod(seclft,60.0d+00)
        endif

        if(doy .ge. 0)then
            year = 1970
 1000       continue
                diy =  leapdy(year)
                if(doy .lt. diy)go to 2000
                doy = doy - diy
                year = year + 1
            go to 1000
        else
            year = 1969
 1100       continue
                diy =  leapdy(year)
                doy = doy + diy
                if( doy .ge. 0 ) go to 2000
                year = year - 1
            go to 1100
        endif
 2000   continue
        doy = doy + 1
        date = year*1000 + doy
        call mnthdy(year,doy,month,day)
c-----
c       the minimum length of the string is 32
c-----
        write(str,110) year,doy,year,month,day,hour,minute,second
  110   format(i4,'/',i3,' ',i4,'/',i2,'/',i2,' ',i2,':',i2,':',f6.3)
c-----
c       guarantee that there are no blanks in the string str
c-----
        do 2100 i=1,32
            if(str(i:i).eq.' ')str(i:i)='0'
 2100   continue
c-----
c       except between date and time
c-----
        str( 9: 9) = ' '
        str(20:20) = ' '
        return
        end

        function leapdy(yr)
        integer*4 yr
        logical t1, t2, t3
        t1 = mod(yr,4).ne.0
        t2 = mod(yr,100).ne.0
        t3 = mod(yr,400).ne.0
        if( .not.t1 .and. t2)then
            isleap = 1
            leapdy = 366
        elseif( .not.t3)then
            isleap = 1
            leapdy = 366
        else
            isleap = 0
            leapdy = 365
        endif
        return
        end

        subroutine mnthdy(year,doy,month,day)
c-----
c       given YEAR and DOY, return MONTH and DAY
c-----
        integer*4 year, doy, month, day
        integer*4 i, dim, leap
        integer*4 dmnth(12)
        data dmnth/31,28,31,30,31,30,31,31,30,31,30,31/
        if(leapdy(year).eq.366)then
            leap = 1
        else
            leap = 0
        endif
        day = doy
        do 100 i=1,12
            month = i
            dim = dmnth(i)
            if(leap.eq.1 .and. i.eq.2)dim = dim + 1
            if(day .le.dim)goto 1000
            day = day - dim 
  100   continue
 1000   continue
        return
        end


        subroutine htoe(year,month,day,hour,minute,second,epoch)
c-----
c       convert calendar date to epoch time since January 1, 1970
c-----
c       year    - I*4   year
c       month   - I*4   month
c       day - I*4   day
c       hour    - I*4   hour
c       minute  - I*4   minute c    second  - I*4   second
c       second  - R*4   seconds
c       epoch   - R*8   time in seconds relative to 00:00 01 Jan 1970
c-----
        integer*4 year, month, day, hour, minute, date, diy
        real*4 second
        real*8 epoch, dtoepo
        integer*4 daymon(12)
        data daymon/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     1      304, 334/
        diy = daymon(month) + day
        if(leapdy(year).eq.366 .and. month.gt.2)diy=diy+1
        date = 1000*year + diy
c       write(6,*)'date=',date
        epoch = dtoepo(date) + hour * 3600.0d+00 + 
     1      minute * 60.0d+00 +dble(second)
        return
        end

c-----
c       convert julian date to epoch time
c-----
        function dtoepo(date)
        real*8 dtoepo
        integer*4 date, diy, cnt, days

        cnt = date / 1000
        days = 0
        if (cnt .gt. 1970)then
 1000       continue
            cnt = cnt -1
            if(cnt.lt.1970)go to 2000
                days = days + leapdy(cnt)
            go to 1000
        else if (cnt .lt. 1970)then
 1100       continue
            if(cnt.ge.1970)goto 2000
                days = days - leapdy(cnt)
                cnt = cnt + 1
            go to 1100
        endif
 2000   continue
        diy = (date -1) / 1000
        diy = (date -1 ) -  1000*diy
c       write(6,*)'days=',days,' diy=',diy
        dtoepo = (days + diy) * 86400.0d+00
        return
        end
            
        subroutine scmxmn(x,npts,depmax,depmin,depmen,indmax,indmin)
c-----
c       get extremal values of the time series
c-----
        real*4 x(*)
        real*4 depmax,depmin,depmen
        integer*4 npts
        depmax = -1.0e+38
        depmin =  1.0e+38
        sum = 0.0
        do 1000 i=1, npts
            if( x(i) .gt. depmax) then
                depmax = x(i)
                indmax = i-1
            endif
            if( x(i) .lt. depmin) then
                depmin = x(i)
                indmin = i-1
            endif
            sum = sum + x(i)
 1000   continue
        if(npts.gt.0)then
            depmen = sum / npts
        else
            depmen = -12345.
        endif
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
        function mnmarg()
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MNMARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        implicit none
        integer iargc
        integer mnmarg
            mnmarg = iargc()
        return
        end
        subroutine mgtarg(i,name)
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINE: MGTARG                                             c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c       return the i'th command line argument
c
c       This version works on SUN, IBM RS6000
c-----
        implicit none
        integer i
        character name*(*)
            call getarg(i,name)
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
        subroutine npow2(npts)
c-----  
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
        end
