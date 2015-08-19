        PROGRAM eSPec
      IMPLICIT NONE
c     **
*     ..
*     Purpose
*     =======
*
*     ..
*     Arguments
*     =========
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes - email:freddy@ufmg.br
*     Amary Cesar Ferreira - email:yrra@zeus.qui.ufmg.br
*     Vinicius Vaz da Cruz - email: viniciusvcruz@gmail.com
*     Viviane Costa Felicissimo - 
*     ..
*     Historic
*     ========
*     (13/03/2003) First version of espec written by Freddy
*     In this first version the program is able to canstruct the spectrum by 
*     time dependent formalism.
*     (20/06/2003) Inclusion of the time independent spectrum calculation 
*     by Freddy.
*     (30/06/2003) Inclusion of pulse control by Freddy.
*     (-/09/2014) Some bugs fixed in version 07 by Vinicius and Freddy
*     and inclusion of cross term hamiltonian
* 
c     **
c     ** Scalar arguments
      LOGICAL       ABSORB, CHANGE, PRTEGVC, PRTPOT, PRTVEFF, PRTEIGVC2
      LOGICAL       PRTREGION,PRTONLYREIM
      LOGICAL       PRTDIPL, PRTPULS, LHWHM, WARN, FSTOP, SHIFTP,PRTABS
      CHARACTER*1   DR
      CHARACTER*2   CHAUX
      CHARACTER*3   TDI
      CHARACTER*5   DIM
      CHARACTER*12  OUTFAT, PRTCRL, TPDIAG, TPPROPG
      CHARACTER*12  TPTRANS, TPWIND, EFC, INIEIGVC, TPDIPL
      CHARACTER*14  INPUT, TPCLC
      CHARACTER*16  TPABSOR
      CHARACTER*25  POTFILE, POTFILEF, FILEAUX, DMFILE, POTCH, POTCHF
      CHARACTER*172 TITLE
      INTEGER       I, IIL, IIU, IFL, IFU, INFO, J, K, MC, NDAUX,NREG
      INTEGER       MF, MP, MPR, MPI, NC, NFS, NSEED, NIS, ND, NDT
      INTEGER       NT, NREORT, NSHOT, NPR, MAXINT, KP, KL, KLX, ISC
      REAL*8        ABSTOL, CS, DE, DT, DTF, T, TI, TF, TOL, SC, WP
      REAL*8        FATE, FATL, FATT, E0, T0, TD, TP, OMG, SNI, E1I
      REAL*8        CHWHM, ALPHA, ADE, RABI, XM, XMIN, XMINB
      REAL*8        GAMMA, OMEGA, DMX, T0X, TPX, AA, QC, DELQ
c     ** 
c     ** Parameters 
      INCLUDE "param.f" !This file contain the parameters used by the main program, 
c for change in the dimensions have a look on. There you can change everything you 
c want and recompile espec using 'make veryclean; make' command. 
      REAL*8        ZERO, ONE, TWO, SHBAR, SHCGS, CCGS, A0CGS, A0A, TAU
      REAL*8        TWOPI
      PARAMETER     (
     &     ZERO  = +0.0D+0, 
     &     ONE   = +1.0D+0, 
     &     TWO   = +2.0D+0,
     &     SHBAR = +1.0D+0, 
     &     SHCGS = +6.62606876D-27,
     &     CCGS  = +2.99792458D+10, 
     &     A0A   = +5.291772083D-1,  
     &     A0CGS = +5.291772083D-9,
     &     TAU   = +2.4188843243D-17,
     &     TWOPI = +6.2831853071795864769252867663D+0
     &     )
c     **
c     ** Array arguments
      INTEGER       IWK(MXDCT),  NP(MXDM)
      REAL*8        CSTI(MXCST), CSTF(MXCST), CSTFB(MXCST), CSTDM(MXCST)  
      REAL*8        EIGVL(LMTREORT), PL(MXDM), SH(MXDM), SHM(MXDM)
      REAL*8        SM(MXDM), SMF(MXDM), SHF(MXDM),U1(MXDCT), V1(MXDCT)
      REAL*8        VPOT(MXDCT), VABC(MXDCT), DM(MXDCT), XI(MXDM)
      REAL*8        XF(MXDM), XP(MXDM), WK1(MXDCT), WK2(MXDCT)
      REAL*8        WK3(MXDCT), WK4(MXDCT), AL(MXDM), AR(MXDM)
      REAL*8        VOI(MXDM), X0(MXDM), rK(MXDM), rA(MXDM), VAR(MXDCT)
      REAL*8        EIGVC(MXDCT,LMTREORT), LANCZ(MXDCT,LMTREORT)
      REAL*8        RANGE(5,7)
c      REAL*8        WMTX(MXDCT,10)
c     **
c     ** External functions 
cdel      CHARACTER
      INTEGER       ICHLENGTH, NISPG, NFSPG
      REAL*8        DLAMCH
c     **
c     ** External subroutines 
      EXTERNAL      RDINPUT, MAKECST, RDPT, GVPOT, MTRXDIAG, LANCZS
      EXTERNAL      DFDXI, RIF, PSOD, PLNZ, SPECTRUMTD, SPECTRUMTI
      EXTERNAL      GETCORR, PRTPT, PRTEIGVC, PPSOD, PSPOFFT, PPLNZ
      EXTERNAL      DIPLMMT, HWHM, ABSORBINGBC, S2PPSOD, S2PPABM2
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT
c     ..
c     .. Starting program
      WRITE(*,1001)
      WRITE(*,1002)
      WRITE(*,1001)
      WRITE(*,1003)'Initializing...'
c     .. 
c     .. Starting default logical values
      CHANGE    = .FALSE.
      PRTEGVC   = .FALSE.
      PRTPOT    = .FALSE.
      PRTVEFF   = .FALSE.
      PRTEIGVC2 = .FALSE.
      PRTPULS   = .FALSE.
      PRTDIPL   = .FALSE.
      ABSORB    = .FALSE.
      LHWHM     = .FALSE.
      WARN      = .FALSE.
      FSTOP     = .TRUE.
      SHIFTP    = .TRUE.
      PRTREGION = .FALSE.
      PRTONLYREIM = .FALSE.
c     .. 
c     .. Starting character default values
      INPUT    = 'input.spc   ' 
      EFC      = '.NONE   '
      OUTFAT   = '.AU   '
      PRTCRL   = '.NO   '
      TPCLC    = '.SPECTRUM   '
      TDI      = '.NC   '
      TPDIAG   = '.MTRXDIAG   '
      POTCH    = '.NULL   '
      POTCHF   = '.NULL   '
      POTFILE  = '.NULL   '
      POTFILEF = '.NULL   '
      TPPROPG  = '.NULL   '
      TPTRANS  = '.ONE   '
      TPWIND   = '.NONE   '
      TPDIPL   = '.NULL   '
      INIEIGVC = '.CALC   '
c     .. 
c     .. Starting scalar default values
      IIL = ZERO
      IIU = ZERO
      IFL = ZERO
      IFU = ZERO
      MAXINT = 300
      MF     = +1.3D+1
      MP     = +1.0D+1
      NIS    = ONE
      NFS    = ZERO
      NSEED  = ONE
      NREORT = LMTREORT
      NISPG  = ZERO
      NFSPG  = ZERO
      KP     = ONE
      KL     = ZERO
      KLX    = ZERO
      ABSTOL = 1D-6!TWO*DLAMCH('S')
      CS     = ONE
      DT     = 1.0D-4
      FATE   = ONE
      FATL   = ONE
      FATT   = ONE
      TI     = ZERO
      TF     = DT
      TOL    = +3.0D+0
      CHWHM  = ZERO
      NSHOT  = +1.0D+1
      NPR    = ONE
      DE     = ZERO
      E1I    = ZERO
      ADE    = ZERO
      E0     = ZERO
      OMG    = ZERO 
      TP     = ZERO
      SNI    = ZERO
      T0     = ZERO
      TD     = ZERO
      GAMMA  = ZERO
      OMEGA  = ZERO
      DMX    = ZERO
      TPX    = ONE
      T0X    = ZERO
c     ..
c     .. Starting integer arrays values
      DATA IWK / MXDCT* ZERO /
      DATA NP  / MXDM* ONE / 
c     .. 
c     .. Starting real arrays values
c      DATA ((EIGVC(I,J),I=1,MXDCT,1),J=1,LMTREORT,1) / MXAUX * 0.0E+0 /
      DATA AL / MXDM* ZERO /
      DATA AR / MXDM* ZERO /
      DATA rK / MXDM* ZERO /    
      DATA rA / MXDM* ZERO /      
      DATA X0 / MXDM* ZERO /
      DATA CSTI  / MXCST* ZERO /
      DATA CSTF  / MXCST* ZERO /
      DATA CSTFB / MXCST* ZERO /
      DATA DM    / MXDCT* ZERO /
      DATA VAR   / MXDCT* ZERO /
      DATA EIGVL / LMTREORT* ZERO /
      DATA PL  / MXDM* ONE /
      DATA SH  / MXDM* ZERO /
      DATA SHM / MXDM* ZERO /
      DATA SM  / MXDM* ZERO /
      DATA SMF / MXDM* ONE /
      DATA SHF / MXDM* ONE /
      DATA U1  / MXDCT* ZERO /
      DATA V1  / MXDCT* ZERO /
      DATA VABC / MXDCT* ZERO /
      DATA VPOT / MXDCT* ZERO /
      DATA VOI  / MXDM* ZERO /
      DATA XI   / MXDM* ZERO /
      DATA XF   / MXDM* ZERO /
      DATA XP   / MXDM* ZERO /
      DATA WK1  / MXDCT* ZERO /
      DATA WK2  / MXDCT* ZERO /
      DATA WK3  / MXDCT* ZERO /
      WRITE(*,*)'Initialized.'
c     ..
c     .. Getting input parameters
      WRITE(*,1003)'Getting input parameters...'
      CALL RDINPUT(ABSORB, CHANGE, PRTEGVC, PRTPOT, PRTVEFF, 
     &     PRTEIGVC2, PRTPULS, PRTDIPL, LHWHM, FSTOP, SHIFTP, INPUT,  
     &     TITLE, DIM, FILEAUX, POTCH, POTCHF, POTFILE, POTFILEF, 
     &     PRTCRL, TDI, TPABSOR, TPCLC, 
     &     TPDIAG, TPPROPG, TPTRANS, TPWIND, TPDIPL, EFC, INIEIGVC, 
     &     DMFILE, IIL, IIU, NIS, IFL, IFU, MF, MP, NFS, NSEED, NSHOT, 
     &     NPR, MXCST, MXDM, NISPG, NFSPG, NREORT, MAXINT, KP, KL, KLX, 
     &     NP, ABSTOL, DT, DTF, TI, TF, TOL, WP, E0, TP, TD, T0, SNI, 
     &     OMG, CHWHM, DE, E1I, ADE, GAMMA, OMEGA, DMX, TPX, T0X, AA, 
     &     QC, DELQ, CSTI, CSTF, CSTFB, CSTDM, SM, SMF, XI, XF, AL, AR,
     &     rK, rA, X0, VOI, VAR,SHF,PRTREGION,NREG,RANGE,PRTONLYREIM)
C,PRTABS)
      WRITE(*,*)'Input parameters got.'
c
      IIL = IIL + ONE
      IIU = IIU + ONE
      IFL = IFL + ONE
      IFU = IFU + ONE
      NISPG = NISPG + ONE
      NFSPG = NFSPG + ONE
c     ..
c     .. Calculating the total dimension
      NT = ONE
      ND = ZERO
      DO I=1,MXDM,1
         IF(NP(I).EQ.1)GOTO 10
         ND = ND + ONE
         NT = NT*NP(I)
      ENDDO
 10   CONTINUE
      IF(DIM(1:4).EQ.'.2DC' .OR. DIM(1:4).EQ.'.2D' .AND.
     &     (POTCH(1:8).EQ.'.LEPS_CC' .OR. POTCHF(1:8).EQ.'.LEPS_CC')
     &     )THEN
         NDAUX = ND + ONE
      ELSE
         NDAUX = ND
      ENDIF
c
      IF(TPPROPG(1:3).EQ.'.PF' .OR.
     &     TPPROPG(1:4).EQ.'.PPF' .OR.
     &     TPPROPG(1:5).EQ.'.P2PF')THEN
         ISC = NINT(LOG(FLOAT(NT))/LOG(TWO))
         NC = TWO**ISC
         IF((NT - NC).NE.ZERO)THEN
            IF(CHANGE .AND. DIM(1:3).EQ.'.1D' .AND. 
     &           POTCH(1:5).EQ.'.FILE')THEN
               NP(1) = NC
               NT = NC
               WRITE(*,*)'-> Warning! Total number of points changed to'
     &              , NC
            ELSE
               WRITE(*,*)'This progation requires FFT please use a grid 
     &of points proportional to 2^M'
            ENDIF
         ENDIF
      ENDIF
c
      WRITE(*,*)
      WRITE(*,*)'Input description:'
      WRITE(*,*)'""""""""""""""""""'
      WRITE(*,*)' General:'
      WRITE(*,*)' ````````'
      NC = ICHLENGTH(TITLE, 5)
      IF(NC.NE.ZERO)WRITE(*,*)'    Title: ',TITLE(1:NC)
      NC = ICHLENGTH(TPCLC, 0)
      IF(NC.GT.172)NC=172
      WRITE(*,*)'    This is a ', TPCLC(2:NC),' calculation;' 
c
      IF(TPCLC(1:9).NE.'.ONLYSPEC')THEN
         WRITE(*,*)'    Using ',TDI(2:3),' formalism;' 
c
         NC = ICHLENGTH(DIM, 0)
         WRITE(*,'(A16,A4,A3)')'    Dimension = ',DIM(2:NC),'  ;'
c
         IF(DIM(1:5).EQ.'.2DCT')THEN
            WRITE(*,'(A82)') "This calculation includes a cross term 
     &derivative in the kinetic energy operator"
            WRITE(*,*) " Cross term value ",SHF(3)
            WRITE(*,*) " Cross term mass  ",SM(3)
            WRITE(*,*)
         ENDIF
c
         NC = ICHLENGTH(TPCLC, 0)
         WRITE(*,1004)('R',I,NP(I), 
     &        I=ND,1,-1)
c
         WRITE(*,*)'    Total number of points =',NT,';'
         IF(CHANGE)THEN
            WRITE(*,*)'    Changing environment turn on;'
         ELSE
            WRITE(*,*)'    Changing environment turn off;'
         ENDIF
         IF(SHIFTP)THEN
            WRITE(*,*)'    Shift potential environment turn on;'
         ELSE
            WRITE(*,*)'    Shift potential environment turn off;'
         ENDIF
         IF(FSTOP)THEN
            CONTINUE
         ELSE
            WRITE(*,*)'    DONOTSTOP environment turn on;'
            WRITE(*,*)'       The program will try to go to the end even 
     &if mistakes are found.'
         ENDIF
c
c         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*)' System:'
         WRITE(*,*)' ```````'
c     
 5001    FORMAT(4X,3(1X,A13,I1,1X,A1,1X,F9.6,A1))
 5002    FORMAT(4X,3(1X,A20,I1,1X,A1,1X,F9.6,A1))
         WRITE(*,5001)('Reduced mass_',I,'=', SM(I),';',I=NDAUX,1,-1)
         IF(SMF(1).NE.ONE .OR. SMF(2).NE.ONE .OR. 
     &        SMF(3).NE.ONE)THEN
            WRITE(*,5002)('Rescaling mass factor_',I,'=', SMF(I),';',
     &           I=ND,1,-1)
         ENDIF
c     
c     .. Calculate step double precision array sh(i)
         IF(POTCH(1:5).NE.'.FILE' .OR. POTCHF(1:5).NE.'.FILE' .OR.
     &        TPCLC(1:9).EQ.'.ONLYSPEC')THEN
            DO I=1,ND,1
               SH(I) = (XF(I) - XI(I))/(NP(I) - ONE)
            ENDDO
            WRITE(*,*)'    Grid ranges:'
            WRITE(*,1007)
            WRITE(*,*)'    | Dimension |   Initial/A    |    Final/A',
     &           '     |     Step/A     |'
            WRITE(*,*)'    =========================================',
     &           '======================='
            DO I=1,ND,1
               WRITE(*,1008)I, XI(I),XF(I),SH(I)
               WRITE(*,1007)
            ENDDO
         ENDIF
      ENDIF
c
      IF(TPCLC(1:9).NE.'.ONLYSPEC')THEN
         NC = ICHLENGTH(POTCH, 0)
         IF(POTCH(1:NC).NE.'.NULL')
     &        WRITE(*,*)'    Initial chosen potential: ',POTCH(2:NC),';'
         NC = ICHLENGTH(POTCHF, 0)
         IF(POTCHF(1:NC).NE.'.NULL')
     &        WRITE(*,*)'    Final chosen potential: ',POTCHF(2:NC),';'
      ENDIF
c     
      IF(TPCLC(1:9).EQ.'.ONLYSPEC')THEN
         NC = ICHLENGTH(FILEAUX, 0)
         WRITE(*,*)'       Name of correlation function ',
     &        'input file: ',FILEAUX(1:NC),';'
      ELSE
         NC = ICHLENGTH(POTFILE, 0)
         IF(POTFILE(1:NC).NE.'.NULL')
     &        WRITE(*,*)'       Name of initial potential input file: ',
     &        POTFILE(1:NC),';'
         NC = ICHLENGTH(POTFILEF, 0)
         IF(POTFILEF(1:NC).NE.'.NULL')
     &        WRITE(*,*)'       Name of final potential input file: ',
     &        POTFILEF(1:NC),';'
      ENDIF
c
      IF(POTCH(1:5).NE.'.FILE' .AND. POTCH(1:5).NE.'.NULL')THEN
         WRITE(*,*)'       Parameters initial potential:'         
         WRITE(*,1005)
         WRITE(*,*)'       |  Constant  |      Value     |  Constant  |
     &    Value      |'
         WRITE(*,*)'       ==============================================
     &================'
         DO I=1,MXCST/2,1
            WRITE(*,1006)I,CSTI(I),I+MXCST/2,CSTI(I+MXCST/2)
            WRITE(*,1005)
         ENDDO
      ENDIF
      IF(POTCHF(1:5).NE.'.FILE' .AND. POTCHF(1:5).NE.'.NULL')THEN
         WRITE(*,*)'       Parameters final potential:'
         WRITE(*,1005)
         WRITE(*,*)'       |  Constant  |      Value     |  Constant  |
     &    Value      |'
         WRITE(*,*)'       ==============================================
     &================'
         DO I=1,MXCST/2,1
            WRITE(*,1006)I,CSTF(I),I+MXCST/2,CSTF(I+MXCST/2)
            WRITE(*,1005)
         ENDDO
      ENDIF
c     ..
c     .. Diagonalization description
      IF(POTCH(1:5).NE.'.FILE' .AND. INIEIGVC(1:5).NE.'.CALC')THEN
         WRITE(*,*)'-> Warning! Replacing *INIEIGVC to .CALC'
      ENDIF
c     ..
c     .. Diagonalization description
      IF(INIEIGVC(1:5).EQ.'.CALC' .AND. TPCLC(1:9).NE.'.ONLYSPEC'
     &     .AND. TPCLC(1:10).NE.'.COLLISION')THEN
c         WRITE(*,*)
         NC = ICHLENGTH(TPDIAG, 0)
         WRITE(*,*)
         WRITE(*,*)' Initial state(s):'
         WRITE(*,*)' `````````````````'
         WRITE(*,*)'    First initial state =',
     &        IIL-1,'; last initial state =',IIU-1,';'
         WRITE(*,*)'       total number of initial states =',NIS,';'
         WRITE(*,*)'    Diagonalization of the hamiltonian ',
     &        'will be done by: ',TPDIAG(2:NC),';'
         IF(TPDIAG(1:7).EQ.'.LANCZS')THEN
            WRITE(*,*)'        with tolerance =',ABSTOL,'; ', 
     &           'and ', NREORT,' partial reorthogonalizations;'
         ELSE
            WRITE(*,*)'        with tolerance =',ABSTOL,';'
         ENDIF
         IF(TDI(1:3).EQ.'.TI')THEN
            WRITE(*,*)
            WRITE(*,*)' Final state(s):'
            WRITE(*,*)' ```````````````'
            WRITE(*,*)'    First final state =',
     &        IFL-ONE,'; last final state =',IFU-ONE,';'
            WRITE(*,*)'       total number of final states =',NIS,';'
         ENDIF
      ENDIF
c     ..
c     .. Description of the collison parameters

      IF(TPCLC(1:10).EQ.'.COLLISION')THEN
        WRITE(*,*)
        WRITE(*,*)' Collision parameters:'
        WRITE(*,*)' `````````````````````'
        NC = ZERO
        DO I=1,MXDM,1
           IF(NP(I).GT.ONE)NC = NC + 1
        ENDDO
        WRITE(*,*)'    Collison energy =',(rK(I),' a.u., ',I=1,NC,1)
        WRITE(*,*)'    Gaussian size =',(rA(I),' a.u., ',I=1,NC,1)
        WRITE(*,*)'    Initial possition of the WP =',
     &       (X0(I),' a.u., ',I=1,NC,1)
      ENDIF

c     ..
c     .. Propagation description
      IF(TDI(1:3).EQ.'.TD')THEN
         WRITE(*,*)
         WRITE(*,*)' Propagation:'
         WRITE(*,*)' ````````````'
         NC = ICHLENGTH(TPPROPG, 0)
         WRITE(*,*)'    The propagation will be done by: ',
     &        TPPROPG(2:NC),'; energy and norma tolerance =',TOL,'%;'
         IF(TPPROPG(1:5).EQ.'.PLNZ') WRITE(*,*)'    Number of ',
     &        'vectors in lanczos propagation =',MP,';'
         WRITE(*,*)'    It will be propagate state:',NISPG-ONE,';' 
         WRITE(*,*)'    Propagation initial time =',TI,
     &        '; final time =',TF,';'
         WRITE(*,*)'     and delta time =',DT,';'
         NC = ICHLENGTH(EFC, 0)
         IF(TPPROPG(1:3).EQ.'.PP' .AND. EFC(1:NC).NE.'.NONE' 
     &        .AND. TPDIPL(2:ICHLENGTH(TPDIPL, 0)).NE.'.NULL' 
     &        .OR. TPPROPG(1:3).EQ.'.P2' )THEN
            NC = ICHLENGTH(EFC, 0)
            WRITE(*,*)'    Using laser pulse: ',EFC(2:NC),
     &           '; with intensity (W/cm^2) = ', E0,';' 
            WRITE(*,*)'       frequency/eV = ', OMG,'; ',
     &           'pulse time/fs (HWHM) = ', TP,';' 
            WRITE(*,*)'       phase = ', SNI,';',
     &           ' max. pulse time = ', T0
            NC = ICHLENGTH(TPDIPL, 0)
            WRITE(*,*)'    Dipole moment: ',TPDIPL(2:NC),';'
            IF(TPDIPL(1:NC).EQ.'.READ')THEN
               NC = ICHLENGTH(DMFILE, 0)
               WRITE(*,*)'       Name of dipole moment input file: ',
     &              DMFILE(1:NC),';' 
            ENDIF
         ENDIF
         NC = ICHLENGTH(TPTRANS, 0)
         WRITE(*,*)'    It will be used transition operator: ',
     &        TPTRANS(2:NC),';'
         IF(ABSORB)THEN
            NC = ICHLENGTH(TPPROPG, 0)
            WRITE(*,*)'    Absorbing bondary conditions ', 
     &           TPABSOR(2:NC+1),' turn on;'
         ENDIF
      ENDIF
c     ..
c     .. FFT description
      IF(TPCLC(1:9).EQ.'.SPECTRUM' .AND. TDI(1:3).NE.'.TI' 
     &     .OR. TPCLC(1:9).EQ.'.ONLYSPEC')THEN
         WRITE(*,*)
         WRITE(*,*)' Spectra:'
         WRITE(*,*)' ````````'
         NC = ICHLENGTH(TPWIND, 0)
         IF(TPWIND(1:NC).NE.'.NONE')THEN
            WRITE(*,*)'    Windowing by function: ',TPWIND(2:NC),';'
            WRITE(*,*)'    Windowing parameter =',WP,';'
         ENDIF
         NC = TWO**MF
         IF(NC.GT.MXDCT)THEN          
            WRITE(*,*)'<<<>>> Number of points in the correlation ',
     &           'function, ',NC,', bigger than maximum number of ',
     &           'points in discretization,',MXDCT,'<<<>>>'
            IF(FSTOP) STOP
         ELSE
            
            WRITE(*,*)'    Number of points in the Fourier transform =',
     &           NC,';'
         ENDIF
         IF(LHWHM)THEN
            WRITE(*,*)'    Half width at half maximum of the system =',
     &           CHWHM,' eV;'
            WRITE(*,*)'    Shift in the spectrum energies =',DE,' eV;'
            WRITE(*,*)'    Energy of the initial state =', E1I,' eV;'  
            WRITE(*,*)'    Difference of energy minimum-minimum',
     &           ' between the two electronic states =',ADE,' eV;'
            

         ENDIF
      ENDIF
c
      WRITE(*,*)
      WRITE(*,*)' Prints:'
      WRITE(*,*)' ```````'
      IF(PRTCRL(1:4).EQ.'.YES')THEN
         WRITE(*,*)'    Print correlation function switch on;'
      ELSEIF(PRTCRL(1:8).EQ.'.PARTIAL')THEN
         WRITE(*,*)'    Print partial correlation function switch on;'
      ELSE
         WRITE(*,*)'    Print partial correlation function switch ',
     &        'off;'
      ENDIF
c     
      IF(PRTEGVC)THEN
         WRITE(*,*)'    Print eigenvectors switch on;'
      ELSE
         WRITE(*,*)'    Print eigenvectors switch off;'
      ENDIF
c
      IF(PRTREGION)THEN
         WRITE(*,*)'    Print only regions of the wavepacket'
         DO I=1, NREG,1
            WRITE(*,*) '    region ',I
            WRITE(*,'(6ES15.7)') (RANGE(I,J),J=1,6,1)
            WRITE(*,*)
         ENDDO
      ENDIF
c     
      IF(PRTPOT)THEN
         WRITE(*,*)'    Print potential switch on;'
      ELSE
         WRITE(*,*)'    Print potential switch off;'
      ENDIF
c
      IF(TDI(1:3).EQ.'.TD')THEN
         IF(TPPROPG(1:3).EQ.'.PP')THEN
            IF(PRTVEFF)THEN
               WRITE(*,*)'    Print effetive potential time by time',
     &              ' switch on;'
            ELSE
               WRITE(*,*)'    Print effetive potential time by time',
     &              ' switch off;'
            ENDIF
            IF(PRTEIGVC2)THEN
               WRITE(*,*)'    Print norm of eigenvector time by time', 
     &              ' switch on;'
            ELSE
               WRITE(*,*)'    Print norm of eigenvector time by time', 
     &              ' switch off;'
            ENDIF 
            IF(PRTPULS)THEN
               WRITE(*,*)'    Print pulse time by time switch',
     &              ' on;'
            ELSE
               WRITE(*,*)'    Print pulse time by time switch',
     &              ' off;'
            ENDIF
            IF(PRTDIPL)THEN
               WRITE(*,*)'    Print dipole momentum switch on;'
            ELSE
               WRITE(*,*)'    Print dipole momentum switch off;'
            ENDIF
         ENDIF
      ENDIF
c     ..
c     .. redireciona out
      IF(TPCLC(1:9).EQ.'.ONLYSPEC') GOTO 400
c     ..
c     .. Converting constants to a.u. units
      WRITE(*,1003)'Converting constants to atomic units...'
c      WRITE(*,*)
      WRITE(*,*)'Constants in a.u.'
      WRITE(*,*)'"""""""""""""""""'
      WRITE(*,*)'System:'
c     ..
c     .. Calculate mass in a.u.
c     .. m_e/m_p = 1836.1526675(39) +- 2.1E-9 (1998)
      DO I=1,NDAUX,1
         SM(I) = 1.836152667539D+3*SM(I)
         WRITE(*,1009)I,SM(I)
      ENDDO
c     ..
c     .. Calculate positions in a.u.
      DO I=1,ND,1
         XI(I) = XI(I)/A0A
         XF(I) = XF(I)/A0A
         SH(I) = SH(I)/A0A
         write(*,*)'XI(I),XF(I)',XI(I),XF(I)
      ENDDO
      IF(POTCH(1:5).NE.'.FILE' .OR. POTCHF(1:5).NE.'.FILE' )THEN
         DO I=1,ND,1
c            XI(I) = XI(I)/A0A
c            XF(I) = XF(I)/A0A
c            SH(I) = SH(I)/A0A
            SHM(I) = SHF(I)/(2.4D+1*SM(I)*SMF(I)*SH(I)*SH(I))
            WRITE(*,1010)I,XI(I),I,XF(I),I,SH(I)
            IF(SH(I).GT.4.0D-2)THEN
               WRITE(*,*)'<<<>>> Step bigger than', 4.0D-2,
     &              ' a.u. <<<>>>'
               IF(FSTOP) STOP
            ELSEIF(SH(I).GT.1.0D-2)THEN
               WRITE(*,*)'-> Warning step bigger than', 1.0D-2,
     &              ' a.u. !'
            ENDIF    
         ENDDO 
         IF(DIM(1:4).EQ.'.2DC')THEN
            SHM(3) = SHF(3)/(2.4D+1*SM(3)*SMF(3)*SH(1)*SH(2))
         ENDIF
         IF(DIM(1:5).EQ.'.2DCT')THEN
            SHM(3) = SHF(3)/(9.6D+1*SM(3)*SMF(3)*SH(1)*SH(2))
         ENDIF
      ENDIF
c     ..
      IF(ABSORB)THEN
         WRITE(*,*)'Absorbing bondary conditions:'
         IF(POTCH(1:5).EQ.'.FILE' .OR. POTCHF(1:5).EQ.'.FILE'
     &        .OR. TPABSOR(1:11).EQ.'.SMOOTHW_au' .OR.
     &        TPABSOR(1:11).EQ.'.SMRS_au' .OR. 
     &        TPABSOR(1:7).EQ.'.VOPTIC')THEN
            DO I=1,ND,1
               WRITE(*,1011)I,AL(I),I,AR(I)
            ENDDO 
         ELSE
            DO I=1,ND,1
               AL(I) = AL(I)/A0A
               AR(I) = AR(I)/A0A
               WRITE(*,1011)I,AL(I),I,AR(I)
            ENDDO
         ENDIF
      ENDIF
c     ..
c     .. Getting output conversion factors
      IF(OUTFAT(1:3).EQ.'.AU')THEN
         FATE = ONE
         FATL = ONE
         FATT = ONE
c      ELSEIF(OUTFAT(1:).EQ.'.CGS')THEN
c         FATE = 
c         FATL = 
c         FATT =  
c      ELSEIF(OUTFAT(1:).EQ.'.SI')THEN
c         FATE = 
c         FATL = 
c         FATT = 
      ENDIF
c     ..
c     .. Convert potential parameters to a.u. unit
      IF(POTCH(1:5).NE.'.FILE' .AND. POTCH(1:5).NE.'.NULL')THEN
         WRITE(*,*)'Initial potential:'
         CALL MAKECST(DIM, POTCH, SM, CSTI)
      ENDIF
      IF(POTCHF(1:5).NE.'.FILE' .AND. POTCHF(1:5).NE.'.NULL')THEN 
         WRITE(*,*)'Final potential:'
         CALL MAKECST(DIM, POTCHF, SM, CSTF)
      ENDIF
c     ..
c     .. End of units conversion.
      WRITE(*,*)
      WRITE(*,*)'Units converted.'
c     ..
c     .. Getting initial potential;
      IF(POTCH(1:5).EQ.'.FILE')THEN 
         IF(INIEIGVC(1:5).EQ.'.CALC')THEN 
            IF(TPCLC(1:10).EQ.'.COLLISION')THEN
               WRITE(*,1003)'Computing the initial wave packet...'
               DO I=1,ND,1
                  IF(RK(I).LT.ZERO)THEN
CNEEDS MODIFICATION, EXPECTED VALUE FOR GAUSSIAN ENERGY WRONG <<<<<<<<<<<<<<<<<<< vinicius
C ACTUALLY, THE PROGRAM READS THE ENERGY ASSOCIATED WITH THE PLANE WAVE WITH K0 in the gaussian wavepacket distribution << vinicius
                     RK(I) = -SQRT(-TWO*SM(I)*SMF(I)*RK(I))
                  ELSE
                     RK(I) =  SQRT(TWO*SM(I)*SMF(I)*RK(I))
                  ENDIF
                  SH(I) = (XF(I) - XI(I))/(NP(I) - ONE)
               ENDDO
d               write(*,*)rk(1),rk(2)
               CALL INIT_COND(DIM, INIEIGVC, NP, rK, rA, SH, XI, 
     &              X0, U1, V1, WK1, WK2, VAR)
               WRITE(*,*)'Initial wave packet computed.'
               GOTO 300
            ELSE
               WRITE(*,1003)'Getting initial potential...' 
               CALL RDPT(SHIFTP, POTFILE, NT, ND, XI, XF, XM, VPOT)
               WRITE(*,*)'Initial potential got.'
            ENDIF
         ELSEIF(INIEIGVC(1:5).EQ.'.GETC' .OR. 
     &           INIEIGVC(1:6).EQ.'.READC')THEN
            IF(TPCLC(1:10).EQ.'.COLLISION' .AND. 
     &           TPCLC(1:14).NE.'.COLLISION_RST')THEN
               WRITE(*,1003)'Getting initial complex eigenvector ...' 
               IF(INIEIGVC(1:6).EQ.'.GETC2')THEN
                  CALL RDPT2(POTFILE, NP(2), 1, XI(2), XF(2), U1, V1)  
 1             ELSEIF(INIEIGVC(1:5).EQ.'.GETC' .OR. 
     &              INIEIGVC(1:6).EQ.'.GETC1')THEN
                  CALL RDPT2(POTFILE, NP(1), 1, XI(1), XF(1), U1, V1)
               ELSE
                  WRITE(*,*)'INIEIGVC not found, ',INIEIGVC
                  STOP
               ENDIF
               WRITE(*,*)'Initial complex eigenvector got.'
c     
               WRITE(*,1003)'Computing the initial wave packet...'
c               write(*,*)'rk(1),rk(2)',rk(1),rk(2)
               DO I=1,ND,1
                  IF(RK(I).LT.ZERO)THEN
                     RK(I) = -SQRT(-TWO*SM(I)*SMF(I)*RK(I))
                  ELSE
                     RK(I) = SQRT(TWO*SM(I)*SMF(I)*RK(I))
                  ENDIF
                  SH(I) = (XF(I) - XI(I))/(NP(I) - ONE)
               ENDDO
c               write(*,*)'2 rk(1),rk(2)',rk(1),rk(2)
               CALL INIT_COND(DIM, INIEIGVC, NP, rK, rA, SH, XI, 
     &              X0, U1, V1, WK1, WK2, VAR)
               WRITE(*,*)'Initial wave packet computed.'
            ELSEIF(TPPROPG(1:4).EQ.'.P2P')THEN
               WRITE(*,1003)'Getting initial complex eigenvector A...' 
               CALL RDPT2(POTFILE, NT, ND, XI, XF, U1(1), V1(1))      
               WRITE(*,*)'Initial complex eigenvector A got.'
               WRITE(*,*)
               WRITE(*,1003)'Getting initial complex eigenvector B...' 
               CALL RDPT2(POTFILE, NT, ND, XI, XF, U1(NT+1), V1(NT+1))      
               WRITE(*,*)'Initial complex eigenvector B got.'
            ELSE
               WRITE(*,1003)'Getting initial complex eigenvector...' 
               CALL RDPT2(POTFILE, NT, ND, XI, XF, U1, V1)      
               WRITE(*,*)'Initial complex eigenvector got.'
            ENDIF
         ELSEIF(INIEIGVC(1:5).EQ.'.GETR' .OR. 
     &           INIEIGVC(1:6).EQ.'.READR')THEN
            IF(TPCLC(1:10).EQ.'.COLLISION' .AND.
     &           TPCLC(1:14).NE.'.COLLISION_RST')THEN
               WRITE(*,1003)'Getting initial real eigenvector ...'
               IF(INIEIGVC(1:6).EQ.'.GETR2')THEN
                  CALL RDPTE(POTFILE, NP(2), 1, XI(2), XF(2), XMIN, U1) 
               ELSEIF(INIEIGVC(1:5).EQ.'.GETR' .OR. 
     &              INIEIGVC(1:6).EQ.'.GETR1') THEN
                  CALL RDPTE(POTFILE, NP(1), 1, XI(1), XF(1), XMIN, U1)
               ELSE
                  WRITE(*,*)'INIEIGVC not found, ',INIEIGVC
                  STOP
               ENDIF
               WRITE(*,*)'Initial real eigenvector got.'
c     
               WRITE(*,1003)'Computing the initial wave packet...'
               DO I=1,ND,1
                  write(*,*)'xi(i),xf(i)',xi(i),xf(i)
                  IF(RK(I).LT.ZERO)THEN
                     RK(I) = -SQRT(-TWO*SM(I)*SMF(I)*RK(I))
                  ELSE
                     RK(I) = SQRT(TWO*SM(I)*SMF(I)*RK(I))
                  ENDIF
                  SH(I) = (XF(I) - XI(I))/(NP(I) - ONE)
               ENDDO
               CALL INIT_COND(DIM, INIEIGVC, NP, rK, rA, SH, XI, X0, 
     &              U1, V1, WK1, WK2, VAR)
               WRITE(*,*)'Initial wave packet computed.'
            ELSEIF(TPPROPG(1:4).EQ.'.P2P')THEN
               WRITE(*,1003)'Getting initial real eigenvector A...' 
               CALL RDPTE(POTFILE, NT, ND, XI, XF, XMIN, U1(1))
               WRITE(*,*)'Initial real eigenvector A got.'
               WRITE(*,*)
               WRITE(*,1003)'Getting initial real eigenvector B...' 
               CALL RDPTE(POTFILE, NT, ND, XI, XF, XMIN, U1(NT+1))
               WRITE(*,*)'Initial real eigenvector B got.'
            ELSE
               WRITE(*,1003)'Getting initial real eigenvector...' 
               CALL RDPTE(POTFILE, NT, ND, XI, XF, XMIN, U1)
               WRITE(*,*)'Initial real eigenvector got.'
            ENDIF
         ELSE
            WRITE(*,*)'The program does not how to initiate the wave pack
     &et',INIEIGVC
            STOP
         ENDIF
c     ..
c     .. Calculate step double precision array sh(i);
c     .. and print information concerning to the potential of input file
         WRITE(*,*)
         WRITE(*,*)'Constants of the external input file in a.u.'
         WRITE(*,*)'""""""""""""""""""""""""""""""""""""""""""""'
         DO I=1,ND,1
            SH(I) = (XF(I) - XI(I))/(NP(I) - ONE)
            SHM(I) = SHF(I)/(2.4D+1*SM(I)*SMF(I)*SH(I)*SH(I))
            WRITE(*,1010)I,XI(I),I,XF(I),I,SH(I)

            IF(SH(I).GT.4.0D-2)THEN

               WRITE(*,*)'<<<>>> Step bigger than', 4.0D-2,
     &              ' a.u. program will stop <<<>>>'
               IF(FSTOP) STOP
            ELSEIF(SH(I).GT.1.0D-2)THEN
               WRITE(*,*)'-> Warning step bigger than', 1.0D-2,
     &              ' a.u. !'
            ENDIF
         ENDDO 
c
         IF(DIM(1:4).EQ.'.2DC')THEN
            SHM(3) = SHF(3)/(2.4D+1*SM(3)*SMF(3)*SH(1)*SH(2))
         ENDIF   
         IF(DIM(1:5).EQ.'.2DCT')THEN
            SHM(3) = SHF(3)/(9.6D+1*SM(3)*SMF(3)*SH(1)*SH(2))
         ENDIF
c
         IF(INIEIGVC(1:5).EQ.'.READ' .OR.
     &        INIEIGVC(1:4).EQ.'.GET')THEN
            NISPG = ONE
            IIL = ZERO
            IIU = ZERO
            INFO = 0
            GOTO 300
         ENDIF
      ELSE
         IF(TPCLC(1:10).EQ.'.COLLISION')THEN
            WRITE(*,1003)'Computing the initial wave packet...'
c            print *, '1 rk(1), rk(2)', rk(1), rk(2)
            DO I=1,ND,1
c               write(*,*)'RK(I)',RK(I)
               IF(RK(I).LT.ZERO)THEN
                  RK(I) = -SQRT(-TWO*SM(I)*SMF(I)*RK(I))
               ELSE
                  RK(I) = SQRT(TWO*SM(I)*SMF(I)*RK(I))
               ENDIF
            ENDDO
c            print *, '2 rk(1), rk(2)', rk(1), rk(2)
            CALL INIT_COND(DIM, INIEIGVC, NP, rK, rA, SH, XI, 
     &           X0, U1, V1, WK1, WK2, VAR)
            WRITE(*,*)'Initial wave packet computed.'
            GOTO 300
         ELSE
            WRITE(*,1003)'Generating initial potential...' 
            CALL GVPOT(SHIFTP, POTCH, DIM, NP, XI, XF, XMIN, SM, SH, 
     &           CSTI, VPOT)
            WRITE(*,*)'Initial potential generated.'
         ENDIF
      ENDIF  
c
      IF(PRTPOT)THEN
         WRITE(*,*)
         WRITE(*,*)'Initial potential points:'
         WRITE(*,*)'"""""""""""""""""""""""""'
         CALL PRTPT(POTCH, 6, ND, ZERO, NP, XP, XI, SH, VPOT)
      ENDIF
c
      IF(TPCLC(1:4).EQ.'.PES' .OR. 
     &     TPCLC(1:10).EQ.'.POTENTIAL')GOTO 9999 
c     ..
c     .. Calculate eignvalues and eignvectors of initial state

      WRITE(*,1003)'Calculating eigenvalue(s) and eigenvector(s)'
      WRITE(*,1013)'for the initial state...'
c
      IF(CHANGE)THEN         
         IF(NT.GT.3.0D+4)THEN
            TPDIAG = '.LANCZSG'
            WRITE(*,*)'-> Type of diagonalization changed to: ',
     &           TPDIAG(2:8)
c            WRITE(*,*)'-> .LANCZSG allows only the ground state ',
d     &           'eigenvectors calculation'
         ENDIF
         IF(TPDIAG(1:8).EQ.'.LANCZSG')THEN
            WRITE(*,*)'-> .LANCZSG allows only the ground state ',
     &           'eigenvectors calculation'
            IIL = ONE 
            IIU = ONE
            NIS = ONE
            WRITE(*,*)'-> Total number of eigenvectors changed to ',IIL
         ENDIF
      ENDIF
c
      IF(TPDIAG(1:9).EQ.'.MTRXDIAG')THEN
         WRITE(*,*)'Using complete matrix diagonalization procedure,'
         CALL MTRXDIAG(DIM, IIL, IIU, INFO, MXDCT, NT, ABSTOL, IWK, 
     &        NP, EIGVL, SHM, VPOT, WK1, WK2, EIGVC)
      ELSEIF(TPDIAG(1:8).EQ.'.LANCZSG')THEN
         WRITE(*,*)'Using interative lanczos tri-diagonalization ',
     &        'procedure,'
         CALL LANCZSG(CHANGE, DIM, IIL, IIU, INFO, NREORT, LMTREORT, NT, 
     &     NSEED, MXDCT, ABSTOL, IWK, MAXINT, NP, U1, V1, EIGVL, SHM, 
     &     WK1, WK2, VPOT, WK3, WK4, VAR, EIGVC, LANCZ)
      ELSEIF(TPDIAG(1:7).EQ.'.LANCZS')THEN
         WRITE(*,*)'Using lanczos tri-diagonalization procedure,'
         CALL LANCZS(CHANGE, DIM, IIL, IIU, INFO, NREORT, LMTREORT, NT, 
     &     NSEED, MXDCT, ABSTOL, IWK, NP, U1, V1, EIGVL, SHM, WK1, 
     &     WK2, VPOT, WK3, WK4, VAR, EIGVC, LANCZ)
      ELSE
         WRITE(*,*)'<<<>>> Initial eigenvalue(s) and eigenvector(s)', 
     &        'error <<<>>>'
         WRITE(*,*)TPDIAG
         STOP
      ENDIF
      E1I = EIGVL(NISPG-IIL+1)*27.2113834D0
      IF(INFO.EQ.0)THEN
         WRITE(*,*)
         WRITE(*,*)'Eigenvalues:'
         WRITE(*,*)'""""""""""""'
         WRITE(*,1014)
         WRITE(*,*)'|  Eigenstate  |  Eigenvalues   |'
         WRITE(*,*)'================================='
         DO I = IIL,IIU,1
            WRITE(*,1015)I-1,EIGVL(I-IIL+1)*FATE
            WRITE(*,1014)
         ENDDO
         IF(PRTEGVC .AND. INFO.EQ.0)THEN
            WRITE(*,*)
            WRITE(*,*)'Eigenvectors:'
            WRITE(*,*)'"""""""""""""'
            CALL PRTEIGVC(TPDIAG, 6, MXDCT, ND, NIS, NP, XP, XI, SH, 
     &           EIGVC)
         ENDIF
      ELSE
         WRITE(*,*)'<<<>>> It wasn`t possible obtain eignvalues and ',
     &        'eignvectors correctly <<<>>>'
         WRITE(*,*)'Exiting.'
         IF(FSTOP) STOP
      ENDIF
c     ..
c     .. If desired only energy program stop
      IF(TPCLC(1:7).EQ.'.ENERGY')GOTO 9999 
 300  CONTINUE
c     ..
c     .. Getting final potential
      WRITE(*,*)

      IF(POTCHF(1:6).EQ.'.FILE')THEN
         IF(TPPROPG(1:4).EQ.'.P2P')THEN
            WRITE(*,1003)'   Getting final potential A...'
            CALL RDPTE(POTFILEF, NT, ND, XI, XF, XMIN, VPOT(1))
            WRITE(*,*)'Final potential A got.'
c
            WRITE(*,1003)'   Getting final potential B...' 
            CALL RDPTE(POTFILEF, NT, ND, XI, XF, XMINB, VPOT(NT+1))
            WRITE(*,*)'Final potential B got.'
         ELSE
            WRITE(*,1003)'   Getting final potential...'
            CALL RDPT(SHIFTP, POTFILEF, NT, ND, XI, XF, XMIN, VPOT(1))
            WRITE(*,*)'Final potential got.'
            CONTINUE
         ENDIF    
      ELSE
         IF(TPPROPG(1:4).EQ.'.P2P')THEN
            WRITE(*,1003)'   Generating final potential A...'
            CALL GVPOT(SHIFTP, POTCHF, DIM, NP, XI, XF, XMIN, SM, SH,  
     &           CSTI, VPOT(1)) 
            WRITE(*,*)'Final potential A generated.'
            WRITE(*,1003)'   Generating final potential B...' 
            CALL GVPOT(SHIFTP, POTCHF, DIM, NP, XI, XF, XMINB, SM, SH, 
     &           CSTF, VPOT(NT+1))
            WRITE(*,*)'Final potential B generated.'
         ELSE
            WRITE(*,1003)'   Generating final potential...'
            CALL GVPOT(SHIFTP, POTCHF, DIM, NP, XI, XF, XMINB, SM, SH, 
     &           CSTF, VPOT)
            IF(ABS(XMINB-XMIN).LT.1d-6)WRITE(*,*)'-> Warning: The initia 
     &l and final potential have the same minimum possition.'
            WRITE(*,*)'Final potential generated.' 
         ENDIF      
      ENDIF    
c
      IF(PRTPOT)THEN
         WRITE(*,*)
         IF(TPPROPG(1:4).EQ.'.P2P')THEN
            WRITE(*,*)'Final potential A values:'
            WRITE(*,*)'"""""""""""""""""""""""""'
            CALL PRTPT(POTCHF, 6, ND, ZERO, NP, XP, XI, SH, VPOT)
            WRITE(*,*)
            WRITE(*,*)'Final potential B values:'
            WRITE(*,*)'"""""""""""""""""""""""""'
            CALL PRTPT(POTCHF, 6, ND, ZERO, NP, XP, XI, SH, 
     &           VPOT(NT+1))
         ELSE
            WRITE(*,*)'Final potential values:'
            WRITE(*,*)'"""""""""""""""""""""""'
            CALL PRTPT(POTCHF, 6, ND, ZERO, NP, XP, XI, SH, VPOT)
         ENDIF
      ENDIF  
c
 400  CONTINUE
c     ..
c     .. Change between time-independent  time-dependent formalism 
      IF(TDI(1:3).EQ.'.TI')THEN
c     .............................
c     .. Time-independent formalism
c     .............................
c     ..
c     .. Kepping eigenvalues and eignvector of initial potential
         IF(TPCLC(1:9).EQ.'.SPECTRUM')THEN
            OPEN(1, STATUS='UNKNOWN', FILE='initial_spc.aux')
            WRITE(1,'(A)')'Initial state eigenvalues and eigenvectors:'
            WRITE(1,'(A)')'"""""""""""""""""""""""""""""""""""""""""""'
            WRITE(1,*)
            DO J=1,IIU-IIL+1,1
               WRITE(1,1045)IIL + J - 1
               WRITE(1,1046)EIGVL(J)
               WRITE(1,1047)IIL + J - 1
               WRITE(1,1048)(EIGVC(I, J), I=1,NT,1)
            ENDDO
            CLOSE(1)
         ENDIF
c     ..
c     .. Calculate eignvalues and eignvectors of final state
         WRITE(*,1003)'   Calculating eigenvalue(s) and eigenvector(s)'
         WRITE(*,1013)'      from final state...'
         IF(TPDIAG(1:9).EQ.'.MTRXDIAG')THEN
            WRITE(*,*)'Using complete matrix diagonalization procedure.'
            CALL MTRXDIAG(DIM, IFL, IFU, INFO, MXDCT, NT, ABSTOL, IWK, 
     &           NP, EIGVL, SHM, VPOT, WK1, WK2, EIGVC)
         ELSEIF(TPDIAG(1:7).EQ.'.LANCZS')THEN
            WRITE(*,*)'Using lanczos tri-diagonalization procedure.'
            CALL LANCZS(CHANGE, DIM, IFL, IFU, INFO, NREORT, LMTREORT, 
     &           NT, NSEED, MXDCT, ABSTOL, IWK, NP, U1, V1, EIGVL, 
     &           SHM, WK1, WK2, VPOT, WK3, WK4, VAR, EIGVC, LANCZ) 
         ENDIF
c     
         IF(INFO.EQ.0)THEN
            WRITE(*,*)
            WRITE(*,*)'Eigenvalues:'
            WRITE(*,*)'""""""""""""'
            WRITE(*,1014)
            WRITE(*,*)'|  Eigenstate  |   Eigenvalues   |'
            WRITE(*,*)'=================================='
            DO I = IFL,IFU,1
               WRITE(*,1015)I-1,EIGVL(I-IFL+1)/FATE
               WRITE(*,1014)
            ENDDO
            IF(PRTEGVC)THEN
               WRITE(*,*)
               WRITE(*,*)'Eigenvectors:'
               WRITE(*,*)'"""""""""""""'
               CALL PRTEIGVC(TPDIAG, 6, MXDCT, ND, NFS, NP, XP, XI, SH, 
     &              EIGVC)
            ENDIF
         ELSE
            WRITE(*,*)'<<<>>> It wasn`t able obtain eignvalues and ',
     &           'eignvectors correctly <<<>>>'
            WRITE(*,*)'EXIT'
            IF(FSTOP) STOP
         ENDIF
c
c     .. Making spectrum 
         CALL SPECTRUMTI(IIL, IIU, IFL, IFU, MXDCT, NT, WK1, EIGVL, 
     &        EIGVC)

      ELSEIF(TDI(1:3).EQ.'.TD')THEN
c     ........................... 
c     .. Time-dependent formalism
c     ...........................
         K = NISPG - IIL + ONE
         NC = K - ONE
         IF(NC.EQ.1)THEN
            CHAUX = 'st'
         ELSEIF(NC.EQ.2)THEN
            CHAUX = 'nd'
         ELSEIF(NC.EQ.3)THEN
            CHAUX = 'rd'
         ELSE
            CHAUX = 'th'
         ENDIF
         WRITE(*,1031)NC,CHAUX
c     
         DO I=1,NT,1
            WK1(I) = ZERO
            WK2(I) = ZERO
            WK3(I) = ZERO
            WK4(I) = ZERO
         ENDDO
c     
         IF(TPTRANS(1:4).EQ.'.ONE')THEN
            NDT = ONE
         ELSE
            NDT = TWO*ND
         ENDIF
         WRITE(*,*)NC,CHAUX(1:2),' state started.'
c     
         IF(PRTEGVC .AND. INIEIGVC(1:4).EQ.'.GET' .OR. 
     &       PRTEGVC .AND. INIEIGVC(1:5).EQ.'.READ' .OR.
     &        PRTEGVC .AND. TPCLC(1:10).EQ.'.COLLISION')THEN
            WRITE(*,*)
            IF(TPPROPG(1:4).EQ.'.P2P')THEN
               WRITE(*,*)'Eigenvector A:'
               WRITE(*,*)'""""""""""""""'
            ELSE
               WRITE(*,*)'Eigenvector:'
               WRITE(*,*)'""""""""""""'
            ENDIF
            IF(INIEIGVC(1:5).EQ.'.GETR' .AND. 
     &           TPCLC(1:10).NE.'.COLLISION' .OR. 
     &           INIEIGVC(1:6).EQ.'.READR' .AND. 
     &           TPCLC(1:10).NE.'.COLLISION')THEN
               CALL PRTPT(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, U1(1))
            ELSE
               CALL PRPT2(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, U1(1), 
     &              V1(1))
            ENDIF    
            WRITE(*,*)
C     ??????????????????? Eigenvector B
C     vinicius > i comented everything from 1 to 2
C 1
C            WRITE(*,*)'Eigenvector B:'
C            WRITE(*,*)'""""""""""""""'
C            IF(INIEIGVC(1:5).EQ.'.GETR' .OR. 
C     &           INIEIGVC(1:6).EQ.'.READR')THEN
C               CALL PRTPT(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, U1(NT+1))
C            ELSE
C               CALL PRPT2(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, U1(NT+1), 
C     &              V1(NT+1))
C            ENDIF
         ENDIF
c     
         IF(ABSORB)THEN
            CALL ABSORBINGBC(DIM, TPABSOR, NT, NP, AL, AR, SH, XI, 
     &           XF, VOI, VABC, ALPHA)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
cccccccccccccc Vinicius 18/10/2011
            IF(PRTABS) THEN
            write(*,*)'***********Absorbing Boundary Condition'
            CALL PRTPT(POTCHF, 6, ND, ZERO, NP, XP, XI, SH, VABC)
            write(*,*)'***********'
            ENDIF
cccccccccccccc Vinicius 18/10
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc

c            read(*,*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c            IF(TPPROPG(1:7).EQ.'.P2PSOD')THEN
c                CALL ABSORBINGBC(DIM, TPABSOR, NT, NP, AL, AR, SH, XI, 
c     &           XF, VOI, VABC(NT+1), ALPHA)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               write(*,*)'***********cond. absorvente'
c               CALL PRTPT(POTCH, 6, ND, ZERO, NP, XP, XI, SH, 
c     &              VABC(NT+1))
c               write(*,*)'***********cond. absorvente final'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            ENDIF
cXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         ENDIF
c     
         IF(TPPROPG(1:3).EQ.'.PP' .AND. TPDIPL(1:5).NE.'.NULL'
     &        .OR. TPPROPG(1:3).EQ.'.P2' .AND. TPDIPL(1:5).NE.'.NULL' 
     &        )THEN
            WRITE(*,*)
            IF(TPPROPG(1:4).EQ.'.P2P')THEN
               WRITE(*,1003)'   Getting dipole momentum A...'
            ELSE
               WRITE(*,1003)'   Getting dipole momentum...'
            ENDIF
            CALL DIPLMMT(TPDIPL, DMFILE, NT, ND, CSTI, CSTF, CSTDM, 
     &           XI, XF, SH, DM)
CCC
CC          VINICIUS: calculates the transition dipole moments
            CALL TRANSDIPOLE(IIL, IIU, IFL, IFU, MXDCT, NT, WK1, EIGVL, 
     &           EIGVC,DM)
CC
            WRITE(*,*)
            NC = NINT((XMIN - XI(1))/SH(1) + ONE)
c     write(*,*)'nc',nc,e0,DM(NC),xmin
            RABI = ABS(DM(NC))*SQRT(E0)*2.1132D-9*27.2114D0
            WRITE(*,*)'   Dipole moment =',ABS(DM(NC))
            WRITE(*,*)'   Rabi frequency =',RABI,' eV'
            IF(RABI.GE.ZERO)WRITE(*,*)'   Rabi period =',
     &           TWOPI/RABI*0.6582,' fs'
            WRITE(*,*)'   Period of oscillations of the explicit ', 
     &           'field =', TWOPI/OMG*0.6582,' fs' 
            WRITE(*,*)
            IF(TPPROPG(1:4).EQ.'.P2P')THEN
               WRITE(*,*)'Dipole momentum A got.'
            ELSE
               WRITE(*,*)'Dipole momentum got.'
            ENDIF
            
            WRITE(*,*)
            IF(TPPROPG(1:4).EQ.'.P2P')THEN
               WRITE(*,1003)'   Getting dipole momentum B...'
               CALL DIPLMMT(TPDIPL, DMFILE, NT, ND, CSTI, CSTF, CSTDM, 
     &              XI, XF, SH, DM(NT+1))
               WRITE(*,*)
               NC = NINT((XMINB - XI(1))/SH(1) + ONE) + NT
c               write(*,*)'nc',nc,e0,DM(NC),xminb
               RABI = ABS(DM(NC))*SQRT(E0)*2.1132D-9*27.2114D0
               WRITE(*,*)'Rabi frequency =',RABI,' eV'
               IF(RABI.GE.ZERO)WRITE(*,*)'Rabi period =',
     &              TWOPI/RABI*0.6582,' fs'
               WRITE(*,*)'Period of oscillations of the explicit ', 
     &              'field =', TWOPI/OMG*0.6582,' fs'
               RABI = DMX
               WRITE(*,*)'Rabi frequency X-ray =',RABI,' eV'
               IF(RABI.GE.ZERO)WRITE(*,*)'Rabi period X-ray =',
     &              TWOPI/RABI*0.6582,' fs'
c               WRITE(*,*)'Period of oscillations of the  ', 
c     &              'field =', TWOPI/OMG*0.6582,' fs'
               WRITE(*,*)
               WRITE(*,*)'Dipole momentum B got.'
               WRITE(*,*)
            ENDIF

            IF(PRTDIPL)THEN
               WRITE(*,*)
               IF(TPPROPG(1:4).EQ.'.P2P')THEN
                  WRITE(*,*)'Dipole momentum potential A:'
                  WRITE(*,*)'""""""""""""""""""""""""""""'
               ELSE
                  WRITE(*,*)'Dipole momentum potential:'
                  WRITE(*,*)'""""""""""""""""""""""""""'
               ENDIF
               CALL PRTPT(TPDIPL, 6, ND, ZERO, NP, XP, XI, SH, DM)
               IF(TPPROPG(1:4).EQ.'.P2P')THEN
                  WRITE(*,*)
                  WRITE(*,*)'Dipole momentum potential B:'
                  WRITE(*,*)'""""""""""""""""""""""""""""'
                  CALL PRTPT(TPDIPL, 6, ND, ZERO, NP, XP, XI, SH, 
     &                 DM(NT+1))
               ENDIF
            ENDIF
         ENDIF    
c     
         DO J=1,NDT,2
c     
            IF(INIEIGVC(1:4).NE.'.GET' .AND.
     &           INIEIGVC(1:5).NE.'.READ' .AND. 
     &           TPCLC(1:10).NE.'.COLLISION')THEN
c               write(*,*)'cannot enter here'
c               stop
               DO I=1,NT,1
                  U1(I) = EIGVC(I,K)
                  V1(I) = ZERO
               ENDDO
            ENDIF
c     
            IF(J.EQ.1)THEN
               DR = 'i'
            ELSEIF(J.EQ.3)THEN
               DR = 'j'
            ELSEIF(J.EQ.5)THEN
               DR = 'k'
            ENDIF
c     ..
c     .. Acting transition operator
            IF(TPTRANS(1:9).EQ.'.VELOCITY')THEN
               WRITE(*,1003)'Applying velocity transition operator...'
               CALL DFDXI(DR, NT, ND, NP, CS, SH, PL, WK4, U1) 
            ELSEIF(TPTRANS(1:7).EQ.'.DIPOLE')THEN
               WRITE(*,1003)'Applying dipole transition operator...'
               CALL RIF(DR, NT, ND, NP, CS, XI, XF, SH, PL, WK4, U1)
            ELSEIF(TPTRANS(1:4).EQ.'.ONE')THEN
               CONTINUE
            ELSE
               WRITE(*,*)'<<<>>> Transition operator error <<<>>>'
               IF(FSTOP) STOP
            ENDIF
            IF(PRTEGVC .AND. TPTRANS(1:9).EQ.'.VELOCITY' .OR.
     &           TPTRANS(1:7).EQ.'.DIPOLE')THEN
               WRITE(*,*)'Transition operator applied.'
               WRITE(*,*)
               IF(TPPROPG(1:4).EQ.'.P2P')THEN
                  WRITE(*,*)'Eigenvector A:'
                  WRITE(*,*)'""""""""""""""'
               ELSE
                  WRITE(*,*)'Eigenvector:'
                  WRITE(*,*)'""""""""""""'
               ENDIF
               IF(INIEIGVC(1:5).EQ.'.GETC' .OR.
     &              INIEIGVC(1:6).EQ.'.READC' .OR.
     &              TPCLC(1:10).NE.'.COLLISION')THEN
                  CALL PRPT2(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, U1(1), 
     &                 V1(1))
               ELSE
                  CALL PRTPT(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, U1(1))
               ENDIF
               IF(TPPROPG(1:4).EQ.'.P2P')THEN
                  WRITE(*,*)'Eigenvector B:'
                  WRITE(*,*)'""""""""""""""'
                  IF(INIEIGVC(1:5).EQ.'.GETC' 
     &                 .OR. INIEIGVC(1:6).EQ.'.READC')THEN
                     CALL PRPT2(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, 
     &                    U1(NT+1), V1(NT+1))
                  ELSE
                     CALL PRTPT(TPDIAG, 6, ND, ZERO, NP, XP, XI, SH, 
     &                    U1(NT+1))
                  ENDIF
               ENDIF
            ENDIF
c     ..
c     .. Calculating DTF parameter for total correlation
            WARN = .FALSE.
            IF(TPCLC(1:9).EQ.'.SPECTRUM' .OR.
     &           TPCLC(1:7).EQ.'.PROPAG' .OR. 
     &           TPCLC(1:10).EQ.'.COLLISION' )THEN
 500           CONTINUE
c
               SC  = TWO*(TF - TI)/DT
               NC  = INT(SC/TWO**MF)
c
               IF(NC.LT.ZERO)THEN
                  NC = 2147483645
                  SC = 2147483645
                  WARN = .TRUE.
                  WRITE(*,*)'-> Warning! The program will print the prop
     &agation with an interval',NC*DT,'fs. This correspond to an interva
     &l of',NC,' points in time'
               ENDIF
               IF(ABS(SC - NC*TWO**MF).GE.1.0D-10 .AND. CHANGE)THEN
                  TF = NC*TWO**(MF - ONE)*DT                  
                  WRITE(*,*)'-> Warning! Final time changed to', TF
                  GOTO 500
               ELSEIF(WARN)THEN
                  CONTINUE
               ELSEIF(ABS(SC - NC*TWO**MF).GE.1.0D-10 .AND. 
     &                 TPCLC(1:9).EQ.'.SPECTRUM')THEN
                  WRITE(*,*)'<<<>>> Final time is not adequate for FFT <
     &<<>>>'
                  IF(FSTOP) STOP
               ENDIF
            ENDIF

c     ..
c     .. Propagation of wave function
            WRITE(*,1003)'Propagating eigenvector...'
            MPR = MP + J
            MPI = MP + J + ONE
            IF(TPPROPG(1:5).EQ.'.PSOD')THEN
               WRITE(*,*)'Using SOD procedure.'
               CALL PSOD(ABSORB, DIM, PRTCRL, TPABSOR, INFO, NC, MP, 
     &              MPI, MPR, MXDCT, NT, DT, TI, TF, TOL, NP, SHM, U1, 
     &              WK1, V1, WK2, XI, XF, VPOT, VABC, WK3, VAR, EIGVC)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
            ELSEIF(TPPROPG(1:6).EQ.'.PFSOD')THEN
               WRITE(*,*)'Using SOD procedure with FFT.'
               CALL PSODFFT(ABSORB, DIM, PRTCRL, TPABSOR, INFO, NC, MP, 
     &              MPI, MPR, MXDCT, NT, DT, TI, TF, TOL, NP, SHM, SH, 
     &              SM, U1, WK1, V1, WK2, XI, XF, DM, VPOT, VABC, WK3, 
     &              WK4, VAR, EIGVC) 
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
c            ELSEIF(TPPROPG(1:6).EQ.'.PWSOD')THEN
c               WRITE(*,*)'Using SOD procedure with Wavelet.'
c               CALL PSODWLET(DIM, PRTCRL, INFO, NC, MP, MPI, MPR, MXDCT, 
c     &              NT, DT, TI, TF, TOL, NP, SHM, SH, SM, U1, WK1, V1, 
c     &              WK2, XI, XF, DM, VPOT, WK3, WK4, EIGVC) 
c               IF(INFO.NE.ZERO)THEN
c                  WRITE(*,*)'<<<>>> Propagation operator error ',
c     &                 '<<<>>>'
c                  WRITE(*,*)'INFO =',INFO
c                  STOP
c               ENDIF
c               IF(TPCLC(1:12).EQ.'.PROPAGATION')GOTO 9999
            ELSEIF(TPPROPG(1:6).EQ.'.PPSOD')THEN
               WRITE(*,*)'Using time dependent SOD procedure with pulse 
     &and printing facilities.'
               CALL PPSOD(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, EFC, 
     &              PRTCRL, TPABSOR, INFO, NC, MP, MPI, MPR, MXDCT, NT, 
     &              NSHOT, NPR, KP, NISPG, E0, T0,  TD, TP, OMG, SNI, 
     &              KL, DT, TI, TF, TOL, NP, ND, SH, SHM, U1, WK1, V1, 
     &              WK2, XI, XP, XF, DM, VPOT, VABC, WK3, VAR, EIGVC,
     &              PRTREGION,NREG,RANGE,PRTONLYREIM) 
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ELSEIF(TPPROPG(1:7).EQ.'.P2PSOD')THEN
               WRITE(*,*)'Using time dependent SOD for coupled equations
     &','   with pulse and printing facilities.'
c               IF(GAMMA.GT.1.0D-3)WRITE(*,*)'<<<>>> Alert! this method i 
c     &s, in general, unstablew for GAMMA values bigger than 0.001 <<<>>>
c     &'
               CALL S2PPSOD(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, 
     &              EFC, PRTCRL, TPABSOR, INFO, NC, MP, MPI, MPR, MXDCT, 
     &              NT, NSHOT, NPR, KP, NISPG, KL, E0, T0, TD, TP, 
     &              OMG, SNI, DT, TI, TF, TOL, DELQ, QC, AA, 
     &              NP, ND, SH, SHM, U1, WK1, V1, WK2, XI, XP, 
     &              XF, SM, DM, VPOT, VABC, WK3, VAR, EIGVC)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ELSEIF(TPPROPG(1:8).EQ.'.P2PABM2')THEN
cxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xxxxxxxxx
c234567890123456789012345678901234567890123456789012345678901234567890123456789
               WRITE(*,*)
               WRITE(*,*)'Using time dependent ABM procedure with pulse 
     &and printing ',
     &              '   facilities for two PES coupled by the vibronic o
     &perator.'
               CALL S2PPABM2(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, 
     &              EFC, PRTCRL, TPABSOR, INFO, NC, MP, MPI, MPR, MXDCT, 
     &              NT, NSHOT, NPR, KP, NISPG, KL, KLX, E0, T0, TD, TP, 
     &              OMG, SNI, DT, TI, TF, TOL, GAMMA, OMEGA, DMX, T0X, 
     &              TPX, AA, QC, DELQ, NP, ND, SH, SHM, U1, WK1, V1, 
     &              WK2, XI, XP, XF, DM, VPOT, VABC, WK3, SM, VAR, 
     &              EIGVC)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ELSEIF(TPPROPG(1:7).EQ.'.P2PABM')THEN
cxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xxxxxxxxx
c234567890123456789012345678901234567890123456789012345678901234567890123456789
               WRITE(*,*)'Using time dependent ABM for coupled equations
     & procedure with pulse and', 
     &              '   printing facilities for coupled equation.'
               CALL S2PPABM(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, 
     &              EFC, PRTCRL, TPABSOR, INFO, NC, MP, MPI, MPR, MXDCT, 
     &              NT, NSHOT, NPR, KP, NISPG, KL, KLX, E0, T0, TD, TP, 
     &              OMG, SNI, DT, TI, TF, TOL, GAMMA, OMEGA, DMX, T0X, 
     &              TPX, NP, ND, SH, SHM, U1, WK1, V1, WK2, XI, XP, 
     &              XF, DM, VPOT, VABC, WK3, VAR, EIGVC)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ELSEIF(TPPROPG(1:6).EQ.'.PFSPO')THEN
               WRITE(*,*)'Using SPO procedure with FFT.'               
               CALL PSPOFFT(ABSORB, DIM, PRTCRL, TPABSOR, INFO, NC, MP, 
     &              MPI, MPR, MXDCT, NT, DT, TI, TF, TOL, NP, SHM, SH, 
     &              SM, U1, V1, XI, XF, WK1, VPOT, VABC, EIGVC, WK2, 
     &              WK3, VAR)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ',
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF     
c            ELSEIF(TPPROPG(1:6).EQ.'.PWSPO')THEN
c               WRITE(*,*)'Using SPO procedure with Wavelet.'               
c               CALL PSPOWLET(DIM, PRTCRL, INFO, NC, MP, MPI, MPR, MXDCT, 
c     &              NT, DT, TI, TF, TOL, NP, SHM, SH, SM, U1, V1, XI, 
c     &              XF, WK1, VPOT, EIGVC, WK2)
c               IF(INFO.NE.ZERO)THEN
c                  WRITE(*,*)'<<<>>> Propagation operator error ',
c     &                 '<<<>>>'
c                  WRITE(*,*)'INFO =',INFO
c                  STOP
c               ENDIF
c               IF(TPCLC(1:12).EQ.'.PROPAGATION')GOTO 9999
            ELSEIF(TPPROPG(1:5).EQ.'.PSIL' .OR.
     &              TPPROPG(1:5).EQ.'.PLNZ')THEN
               WRITE(*,*)'Using SIL procedure.'
               CALL PLNZ(ABSORB, CHANGE, DIM, PRTCRL, TPABSOR, INFO, NC, 
     &              LMTREORT, MP, MPI, MPR, MXDCT, NT, ABSTOL, DT, TI, 
     &              TF, TOL, NP, IWK, SHM, VPOT, VABC, U1, V1, XI, XF, 
     &              WK1, WK2, WK3, WK4, VAR, EIGVC)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ', 
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
            ELSEIF(TPPROPG(1:6).EQ.'.PPSIL' .OR.
     &         TPPROPG(1:6).EQ.'.PPLNZ')THEN
               WRITE(*,*)'Using time dependent SIL procedure with pulse 
     &and printing facilities.'
               CALL PPLNZ(ABSORB, CHANGE, PRTVEFF, PRTEIGVC2, PRTPULS, 
     &              DIM, EFC, PRTCRL, TPABSOR, INFO, NC, LMTREORT, MP, 
     &              MPI, MPR, MXDCT, NT, NSHOT, KL, NPR, KP, ABSTOL, DT, 
     &              TI, TF, TOL, E0, T0, TD, TP, OMG, SNI, NP, ND, IWK, 
     &              SH, SHM, VPOT, VABC, U1, V1, XI, XF, XP, DM, WK1, 
     &              WK2, WK3, WK4, VAR, EIGVC)
c               CALL PPLNZ(ABSORB, CHANGE, PRTVEFF, PRTEIGVC2, PRTPULS, 
c     &              DIM, EFC, PRTCRL, TPABSOR, INFO, LC, LMTREORT, MP, 
c     &              MPI, MPR, MXDCT, N, NSHOT, KL, NPR, KP, ABSTOL, DT, 
c     &              TI, TF, TOL, E0, T0, TD, TP, OMG, SNI, NP, ND, IWORK, 
c     &              SH, SHM, VPOT,       U,   V, XI, XF, XP, DM,  WU, 
c     &               WV,  WP, WORK, LNZVC)
               IF(INFO.NE.ZERO)THEN
                  WRITE(*,*)'<<<>>> Propagation operator error ', 
     &                 '<<<>>>'
                  WRITE(*,*)'INFO =',INFO
                  IF(FSTOP) STOP
               ENDIF
            ELSE
               NC = ICHLENGTH(TPPROPG,0)
               WRITE(*,*)'<<<>>> Propagator',TPPROPG(1:NC),
     &              ' does not exist! <<<>>>' 
               IF(FSTOP) STOP
            ENDIF
            WRITE(*,*)'Eigenvector propagated.'
         ENDDO 
c     ..
c     .. Redirect certain types of calculations to end up the program
         IF(TPCLC(1:12).EQ.'.PROPAGATION' .OR.
     &        TPCLC(1:10).EQ.'.COLLISION')GOTO 9999

c     ..
c     .. Calculating total correlation
         WRITE(*,1003)'Calculating total correlation function...'
         MC = TWO**MF/TWO
         DTF = TF/MC
         DO I=1,2*MC,1
            U1(I) = ZERO
            V1(I) = ZERO
         ENDDO
         DO J=1,2*ND,2 
            MPR = MP + J
            MPI = MP + J + ONE
            DO I=1,MC,1
c     ND so deve ser usado para calculos com diferentes momentos de dipolo
               U1(MC+I) = U1(MC+I) + EIGVC(I,MPR)!/ND ND so deve ser usado 
               U1(MC-I+2) = U1(MC+I)
               V1(MC+I) = V1(MC+I) + EIGVC(I,MPI)!/ND
               V1(MC-I+2) = -V1(MC+I)
            ENDDO
            U1(1) = U1(1) + EIGVC(MC,MPR)!/ND
            V1(1) = V1(1) + EIGVC(MC,MPI)!/ND
         ENDDO
         V1(1) = - V1(1)
c
         WRITE(*,*)'Total correlation function calculated.'
         WRITE(*,*)
         WRITE(*,*)'Correlation Function:'
         WRITE(*,*)'"""""""""""""""""""""'
         WRITE(*,1041)
         DO I=1,2*MC,1
            T = -TF + (I - ONE)*DTF
            WRITE(*,1042)T, U1(I), V1(I)
         ENDDO
c     ..
c     .. Redirect certain types of calculations to end up the program
         IF(TPCLC(1:12).EQ.'.CORRELATION')GOTO 9999
c     ..
c     .. Doing the FFT of correlation function
         IF(TPCLC(1:9).EQ.'.SPECTRUM' .AND. TDI.EQ.'.TD')THEN 
c     ..
c     .. Appling the half-width at half-maximum
            IF(LHWHM)THEN
               DE = DE - E1I
               WRITE(*,1003)'Appling the half-width at half-maximum...'       
               CALL HWHM(MF, DTF, TF, CHWHM, DE, U1, V1) 
               WRITE(*,*)'Half-width at half-maximum applied.'
            ENDIF
            WRITE(*,1003)'   Processing spectrum...'       
            CALL SPECTRUMTD(TPWIND, MF, DTF, TF, WP, ADE, U1, V1)
            WRITE(*,*)'Spectrum done.'
         ENDIF
c
      ELSEIF(TPCLC(1:9).EQ.'.ONLYSPEC')THEN
c     ....................... 
c     .. Making just spectrum 
c     .......................
         MC = TWO**MF/TWO
         DTF = TF/MC
         WRITE(*,1003)'   Getting correlation function...'  
         CALL GETCORR(FILEAUX, MF, DTF, TF, U1, V1) 
         WRITE(*,*)'Correlation function got.'
         IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
            WRITE(*,*)
            WRITE(*,*)'Correlation Function:'
            WRITE(*,*)'"""""""""""""""""""""'
            WRITE(*,1041)
            DO I=1,2*MC,1
               T = -TF + (I - ONE)*DTF
               WRITE(*,1042)T, U1(I), V1(I)
            ENDDO
         ENDIF
c     ..
c     .. Appling the half-width at half-maximum
         IF(LHWHM)THEN
            DE = DE - E1I
            WRITE(*,1003)'Appling the half-width at half-maximum...'       
            CALL HWHM(MF, DTF, TF, CHWHM, DE, U1, V1) 
            WRITE(*,*)'Half-width at half-maximum applied.'
         ENDIF
         WRITE(*,1003)'   Processing spectrum...'  
         CALL SPECTRUMTD(TPWIND, MF, DTF, TF, WP, ADE, U1, V1)
         WRITE(*,*)'Spectrum done.'
      ENDIF
c     ..
 9999 CONTINUE
      WRITE(*,*)
      WRITE(*,*)'The eSPec program finished successfully!'
c     ..
 1001 FORMAT(10X,'****************************************************')
 1002 FORMAT(25X,'eSPec version 0.7 beta',//,14X,
     &     'Principal author: ',/,17X,6X,
     &     'Freddy Fernandes Guimaraes.',/,14X,
     &     'Contributors:',/,17X, 
     &     'Amary Cesar, Vinicius Vaz da Cruz,',/,17X,
     &     'Viviane Costa Felicissimo, Viktor Kimberg, ',/,17X,
     &     'Faris Gelmukhanov, and Yasen Velkov.') 
 1003 FORMAT(/,3X,A68)
 1004 FORMAT(1X,'    Number of points in the discretization:',
     &     1X,3(A1,I1,1X,'=',1X,I4,';',1X))
 1005 FORMAT(8X,'-------------------------------------------',
     &     '------------------')
 1006 FORMAT(8X,'|',3X,I3,6X,'|',2X,E12.6,2X,'|',3X,I3,6X,'|',2X,
     &     E12.6,2X,'|')
 1007 FORMAT(5X,'-----------------------------------------------',
     &     '-----------------')
 1008 FORMAT(5X,'|',4X,'R',I1,5X,'|',2X,E12.6,2X,'|',2X,E12.6,2X,'|',
     &     2X,E12.6,2X,'|')
 1009 FORMAT(1X,'Mass_',I1,'/a.u. = ',E12.6,';')
 1010 FORMAT(1X,'R',I1,'_i/a.u. =',1X,E12.6,';',1X,'R',I1,'_f/a.u. =',
     &     1X,E12.6,';',1X,'dR',I1,'/a.u. =',1X,E12.6,';')
 1011 FORMAT(1X,'A',I1,'_i/a.u. =',1X,E12.6,';',1X,'A',I1,'_f/a.u. =',
     &     1X,E12.6,';')     
 1013 FORMAT(3X,A68)
 1014 FORMAT(1X,'---------------------------------')
 1015 FORMAT(1X,'|',3X,I3,8X,'|',2X,E12.6,2X,'|')
 1031 FORMAT(/,33X,'Starting',1X,I1,A2,1X,
     &     'state to be propagated...') 
 1041 FORMAT(7X,'*t/fs*',7X,'*C(t)r*',7X,'*C(t)i*')
 1042 FORMAT(1X,F14.6,3X,F10.7,4X,F10.7)
 1045 FORMAT(3X,'Energy from state',1X,I3)
 1046 FORMAT(3X,E20.12)
 1047 FORMAT(3X,'Eigenvector from state',I3)
 1048 FORMAT(100000000(3X,E20.12))
 1901 FORMAT(/,/,A32)
      STOP
      END
Cxxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xx
C123456789012345678901234567890123456789012345678901234567890123456789012
