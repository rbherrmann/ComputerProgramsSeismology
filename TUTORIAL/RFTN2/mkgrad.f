        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/depref/refdep
        real refdep

        integer iunit, iiso, iflsph, idimen,icnvel
        logical lverby
c-----
c     initialize
c-----
        refdep = 0.0
        s1 = 1./6.0
        s2 = 1./8.0
        dH = 1.0
        H1 = 40.0
        H2 = 40.0
        mmax = 0.0
c-----
c       get the number of layers
c-----
        call gcmdln(N)
c-----
c       get a finer resolution
c------
        N = N * 2
        dH = dH / 2.0
c-----
c      first layer
c-----
        mmax = mmax + 1
        d(mmax) = H1 - N*dH
        a(mmax) = 1./s1
        do i = 1,2*N
           mmax = mmax  + 1
           d(mmax) = dH
           s = (2*N+1 -i)*s1/(2*N + 1 ) +  i*s2/(2*N + 1)
           a(mmax) = 1./s
        enddo
        mmax = mmax + 1
        d(mmax) = H2 - N*dH
        a(mmax) = 1./s2
c-----
c       use the Harkrider routines to get Vs and Density
c----
        ALP = 0.0
        BET = 0.0
        DENS = 0.0
        SALP = 1.0
        SBET = 1.0
        SDENS = 1.0

        do i = 1, mmax
          ALP = a(i)
          BET = 0.0
          DENS = 0.0
          CALL NAFTAB(ALP,BET,DENS,SALP,SBET,SDENS)
          IF(BET.EQ.0.) THEN
              DENS=1.0
          ENDIF
          b(i) = BET
          rho(i) = DENS
          qa(i) = 0.0
          qb(i) = 0.0
          frefp(i) = 1.0
          frefs(i) = 1.0
          etap(i) = 0.0
          etas(i) = 0.0
        enddo

        iunit = 0
        iiso = 0
        iflsph = 0
        idimen = 1
        icnvel = 0
        lverby = .false.
        call putmod(1,'model.out',mmax,'Simple Gradient',
     1      iunit,iiso,iflsph,
     1      idimen,icnvel,lverby)

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
c-----
c       The following code is from David G. Harkrider in his program rocks.f
C       If only the P-velocity is non-zero, the program
c       returns the S-velocity and density
c-----
      SUBROUTINE NAFTAB(ALP,BET,RHO,SALP,SBET,SRHO)
      DIMENSION DVP(2,34),DVS(2,24)
      DATA DVP/9.9500E-01,1.5175E 00,1.0824E 00,1.4813E 00,1.1363E 00,
     *1.4706E 00,1.2151E 00,1.4570E 00,1.4376E 00,1.4615E 00,
     *1.5930E 00,1.5092E 00,1.7291E 00,1.6214E 00,1.8628E 00,
     *1.8561E 00,1.9979E 00,2.2034E 00,2.1268E 00,2.6226E 00,
     *2.1791E 00,2.8630E 00,2.2230E 00,3.0768E 00,2.2656E 00,
     *3.3098E 00,2.3457E 00,3.8236E 00,2.4708E 00,4.6909E 00,
     *2.5804E 00,5.4304E 00,2.6398E 00,5.7769E 00,2.6704E 00,
     *5.8992E 00,2.6975E 00,6.0599E 00,2.7223E 00,6.1683E 00,
     *2.8213E 00,6.5794E 00,2.8809E 00,6.8258E 00,2.9842E 00,
     *7.1986E 00,3.0551E 00,7.4412E 00,3.1732E 00,7.7815E 00,
     *3.2630E 00,8.0486E 00,3.3909E 00,8.3979E 00,3.5421E 00,
     *8.8049E 00,3.7065E 00,9.2254E 00,3.9935E 00,9.9638E 00,
     *4.2246E 00,1.0560E 01,4.5599E 00,1.1420E 01,4.6988E 00,
     *1.1801E 01,4.7789E 00,1.2011E 01/
      DATA DVS/1.4987E 00,4.5433E-04,1.5724E 00,1.6394E-02,1.6612E 00,
     *8.0173E-02,1.7580E 00,2.5805E-01,1.9398E 00,6.3669E-01,
     *2.0583E 00,9.9000E-01,2.1924E 00,1.4252E 00,2.2895E 00,
     *1.8114E 00,2.4065E 00,2.3509E 00,2.4763E 00,2.6876E 00,
     *2.5779E 00,3.1022E 00,2.6738E 00,3.3931E 00,2.7785E 00,
     *3.6246E 00,2.8973E 00,3.8761E 00,3.0015E 00,4.0663E 00,
     *3.1566E 00,4.3207E 00,3.3236E 00,4.5602E 00,3.5516E 00,
     *4.8919E 00,3.8125E 00,5.2744E 00,4.0642E 00,5.6351E 00,
     *4.2642E 00,5.9253E 00,4.4870E 00,6.2366E 00,4.6469E 00,
     *6.4843E 00,4.7818E 00,6.6793E 00/
      NPS=24
      NPP=34
      IF(ALP.NE.0.) GO TO 100
      IF(BET.EQ.0.) GO TO 200
C*****  ALP=0. AND BET.NE.0.
C**** INTERP FOR ALP AND ITS SLOPE WRST. BET
      CALL WSPCCP(RHOX,SRHOX,BET,DVS,2,NPS)
      CALL WSPCCP(ALP,SALP,RHOX,DVP,1,NPP)
      IF(RHO.NE.0) RETURN
C*** INCLUDE RHO AND ITS SLOPE WRST. BET
      RHO=RHOX
      SRHO=SRHOX
      RETURN
C****** ALP.NE.0.
C***  INTERP FOR BET GIVEN ALP, THEN OBTAIN ALP SLOPE WRST. BET
  100 CALL WSPCCP(RHOX,SRHOX,ALP,DVP,2,NPP)
      IF(BET.EQ.0.) CALL WSPCCP(BET,SBETX,RHOX,DVS,1,NPS)
      CALL WSPCCP(ALPSX,SALP,RHOX,DVP,1,NPP)
      IF(RHO.NE.0.) RETURN
C**** INCLUDE RHO AND ITS SLOPE WRST. BET
      CALL WSPCCP(RHO,SRHO,BET,DVS,2,NPS)
      RETURN
C***** INTERP ALP,BET AND THEIR SLOPES WRST. BET
  200 IF(RHO.EQ.0.) RETURN
      IF(RHO.LT.1.) THEN
        RHO=0.92
        ALP=3.78
        BET=1.84
      ELSE
        CALL WSPCCP(ALP,SALP,RHO,DVP,1,NPP)
        CALL WSPCCP(BET,SBET,RHO,DVS,1,NPS)
      ENDIF
      RETURN
      END
      SUBROUTINE WSPCCP(Y,SY,X,DV,IX,NP)
C******************************************************************
C*
C*     X :  INPUT VALUE. EXPECTS DV(IX,I) TO BE MOMATONIC.
C*     Y :  Harmonic mean WEIGHTED SLOPE PIECEWISE CONTINUOUS CUBIC 
C*          POLYNOMIAL INTERPOLATION OF DATA PAIRS DV(2,NP) FOR A GIVEN X.
C*    SY :  INTERPOLATED SLOPE OF Y WRSPT X.
C*    IX :  SPECIFIES X DATA, IE., DV(IX,NP).
C*    NP :  NUMBER OF DATA PAIRS.
C*
C******************************************************************
      DIMENSION DV(2,*)
      IY=2
      IF(IX.EQ.2)IY=1
      ODIF=X-DV(IX,1)
      I1=1
      IF(ODIF.EQ.0.) GO TO 11
      DO 10 I=2,NP
      I1=I
      IF(X.EQ.DV(IX,I)) GO TO 11
      DIF=X-DV(IX,I)
      CHS=DIF*ODIF
      IF(CHS.LT.0.) GO TO 12
      ODIF=DIF
   10 CONTINUE
      RETURN
   11 Y=DV(IY,I1)
      I0=I1-1
      I2=I1+1
      IF(I1.NE.NP)PS2=(DV(IY,I2)-Y)/(DV(IX,I2)-X)
      IF(I1.NE.1) GO TO 20
      SY=PS2
      RETURN
   20 PS1=(Y-DV(IY,I0))/(X-DV(IX,I0))
      IF(I1.NE.NP) GO TO 21
      SY=PS1
      RETURN
   21 SY=0.
      S1I=0.
      IF(PS1*PS2.GT.0.)S1I=.5/PS1+.5/PS2
      IF(S1I.NE.0.)SY=1./S1I
      RETURN
   12 IM1=I1-2
      I0=I1-1
      I2=I1+1
      Y1=DV(IY,I1)
      Y0=DV(IY,I0)
      X1=DV(IX,I1)
      X0=DV(IX,I0)
      PS1=(Y1-Y0)/(X1-X0)
      PS2=PS1
      IF(I1.NE.NP)PS2=(DV(IY,I2)-Y1)/(DV(IX,I2)-X1)
      S1=0.
      S1I=0.
      IF(PS1*PS2.GT.0.) S1I=.5/PS1+.5/PS2
      IF(S1I.NE.0.) S1=1./S1I
      PS0=PS1
      IF(I0.NE.1)PS0=(Y0-DV(IY,IM1))/(X0-DV(IX,IM1))
      S0=0.
      S0I=0.
      IF(PS1*PS0.GT.0.) S0I=.5/PS1+.5/PS0
      IF(S0I.NE.0.) S0=1./S0I
      CALL CUBSLP(Y,SY,X,Y1,Y0,X1,X0,S1,S0)
      RETURN
      END
      SUBROUTINE CUBSLP(Y,S,X,Y1,Y0,X1,X0,S1,S0)
      DTX=X1-X0
      DTXI=1./DTX
      DTX2I=DTXI*DTXI
      DTX3I=DTXI*DTX2I
      DX=X-X0
      C1=Y1-Y0-DTX*S0
      C2=S1-S0
      DTXC2=DTX*C2
      A0=DTX2I*(3.*C1-DTXC2)
      B0=-DTX3I*(2.*C1-DTXC2)
      Y=Y0+DX*(S0+DX*(A0+DX*B0))
      S=S0+DX*(2.*A0+3.*DX*B0)
      RETURN
      END
      FUNCTION LSTCHR(ID)
      CHARACTER*(*) ID
      MAXCHR=LEN(ID)
      J=MAXCHR
      DO 10 I=1,MAXCHR
      IF(ID(J:J).NE.' ') GO TO 20
      J=J-1
   10 CONTINUE
   20 LSTCHR=J
      RETURN
      END
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
      subroutine gcmdln(N)
      integer N
c-----
c     N   I   - number of layers to add above and below Moho.
c               Total number of layers inserted is 2N
c-----
c     parse command line arguments and return control
c     parameters
c
c     requires subroutine mgtarg(i,name) to return
c           the i'th argument in the string name
c
c     and the function numarg() to return the number
c           of arguments excluding program name
c           The first argument is i = 1
c
c-----
      character*20 name
      integer*4 mnmarg
      N = 0

      nmarg = mnmarg()
      i = 0
   11 i = i + 1
      if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-N')then
                  i=i+1
                  call mgtarg(i,name)
                  read(name,'(bn,i10)')N
            endif
            go to 11
   13 continue
      return
      end
