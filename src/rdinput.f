      SUBROUTINE RDINPUT(ABSORB, CHANGE, PRTEGVC, PRTPOT, PRTVEFF, 
     &     PRTEIGVC2, PRTPULS, PRTDIPL, LHWHM, FSTOP, SHIFTP, INFL,  
     &     TITLE, DIM, FILEAUX, POTCH, POTCHF, POTFILE, POTFILEF,  
     &     PRTCRL, TDI, TPABSOR, TPCLC, 
     &     TPDIAG, TPPROPG, TPTRANS, TPWIND, TPDIPL, EFC, INIEIGVC, 
     &     DMFILE, IIL, IIU, NIS, IFL, IFU, MF, MP, NFS, NSEED, NSHOT, 
     &     NPR, MXCST, MXDM, NISPG, NFSPG, NREORT, MAXINT, KP, KL, KLX, 
     &     NP, ABSTOL, DT, DTF, TI, TF, TOL, WP, E0, TP, TD, T0, SNI, 
     &     OMG, CHWHM, DE, E1I, ADE, GAMMA, OMEGA, DMX, TPX, T0X, AA,
     &     QC, DELQ, CSTI, CSTF, CSTFB, CSTDM, SM, SMF, XI, XF, AL, AR,
     &     rK, rA, X0, VOI, VAR)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
      LOGICAL       ABSORB, CHANGE, PRTEGVC, PRTPOT, PRTVEFF, PRTEIGVC2
      LOGICAL       PRTPULS, PRTDIPL, LHWHM, FSTOP, SHIFTP
      CHARACTER*(*) INFL, TITLE, DIM, FILEAUX, POTCH, POTCHF, POTFILE
      CHARACTER*(*) POTFILEF, PRTCRL, TDI, TPABSOR, TPCLC, TPDIAG 
      CHARACTER*(*) TPPROPG, TPTRANS, TPWIND, TPDIPL, EFC, INIEIGVC 
      CHARACTER*(*) DMFILE
      INTEGER       IIL, IIU, NIS, IFL, IFU, MF, MP, MXCST, MXDM 
      INTEGER       NISPG, NFS, NFSPG, NSEED, NSHOT, NPR, NREORT
      INTEGER       MAXINT, KP, KL, KLX
      REAL*8        ABSTOL, DT, DTF, TI, TF, TOL, WP, E0, TP, TD, T0
      REAL*8        SNI, OMG, CHWHM, DE, E1I, ADE
      REAL*8        GAMMA, OMEGA, DMX, TPX, T0X, AA, QC, DELQ
c     **
c     ** Array arguments
cdel      LOGICAL
      INTEGER       NP(*)
      REAL*8        CSTI(*), CSTF(*), CSTFB(*), SM(*), SMF(*), CSTDM(*)
      REAL*8        XI(*), XF(*), AL(*), AR(*), rK(*), rA(*), X0(*) 
      REAL*8        VOI(*), VAR(*)      
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
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (16/03/2003) First version written RDINPUT by Freddy. 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, A0A, FATEEVAU, SHBAR
      PARAMETER(
     &     ZERO     = +0.0D+0, 
     &     ONE      = +1.0D+0,     
     &     A0A      = +5.291772083D-1, !A/a.u.
     &     FATEEVAU = +27.2113834D+0, !eV/a.u. 
     &     SHBAR    = +6.58211889D-1 !ev*fs
     &     )
c     **
c     ** Local scalars 
cdel      LOGICAL
      CHARACTER*72  READAUX
      INTEGER       I, INFLGTH, ITEST, ND
      REAL*8        SKL, SKLX
c     **
c     ** Local arrays 
      REAL*8        WST(MXDM), VFY(10)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER
      INTEGER       ICHLENGTH, ICOMPAR, IOTEST     
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      DATA VFY / 10* ONE /
c
      ND = 0
      INFLGTH = ICHLENGTH(INFL, 0)
      OPEN(1, STATUS='OLD', FILE=INFL(1:INFLGTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.0)THEN
         WRITE(*,*)'<<<>>> Input file "', INFL(1:INFLGTH), 
     &        '" cannot be opened <<<>>>'
         IF(FSTOP) STOP
      ENDIF
c     
      READ(1,2,END=999,ERR=999)READAUX
      ITEST = ICOMPAR(READAUX, '*** eSPec input file ***        ')
      IF(ITEST.NE.0)THEN
         WRITE(*,1001)
         WRITE(*,*)'Stopping!!! Wrong input file'
         IF(FSTOP) STOP
      ENDIF
c     ..
c     .. Select group
 3    READ(1,'(A)',END=999,ERR=999)READAUX
 10   CONTINUE      
c      write(*,*)readaux
      IF(READAUX(1:2).EQ. '**')THEN
         IF(READAUX(1:6).EQ. '**MAIN')THEN
            GOTO 100
         ELSEIF(READAUX(1:4).EQ.'**TI')THEN
            GOTO 200
         ELSEIF(READAUX(1:4).EQ.'**TD')THEN
            GOTO 300
         ELSEIF(READAUX(1:10).EQ.'**SPECTRUM')THEN
            GOTO 400
         ELSEIF(READAUX(1:5).EQ.'**END')THEN
            GOTO 991
         ELSE
            WRITE(*,1001)
            WRITE(*,1002)READAUX
            IF(FSTOP) STOP
         ENDIF
      ENDIF
      GOTO 3
c     ..
c     .. MAIN GROUP 
 100  CONTINUE
      READ(1,'(A)',END=999,ERR=999)READAUX
c      write(*,*)'readaux',readaux
      IF(READAUX(1:2).EQ. '**')GOTO 10
      IF(READAUX(1:1).EQ. '*' .AND. READAUX(2:2).NE. '*')THEN
         IF(READAUX(1:6).EQ.'*TITLE')THEN
            READ(1,2,END=999,ERR=999)TITLE
         ELSEIF(READAUX(1:7).EQ.'*TPCALC')THEN
            READ(1,'(A)',END=999,ERR=999)TPCLC
            IF(TPCLC(1:7).EQ.'.ENERGY')THEN
               TDI = '.TI'
            ELSEIF(TPCLC(1:12).EQ.'.PROPAGATION' .OR.
     &              TPCLC(1:12).EQ.'.CORRELATION' .OR.
     &              TPCLC(1:10).EQ.'.COLLISION')THEN
               TDI = '.TD'
            ELSEIF(TPCLC(1:9).EQ.'.SPECTRUM')THEN
               READ(1,'(A)',END=999,ERR=999)TDI
            ELSEIF(TPCLC(1:9).EQ.'.ONLYSPEC')THEN
               READ(1,'(A)',END=999,ERR=999)FILEAUX
               TDI = '.NC'  
            ELSEIF(TPCLC(1:10).EQ.'.PES'  .OR.
     &              TPCLC(1:10).EQ.'.POTENTIAL')THEN
               TDI = '.NC'  
            ELSE
               WRITE(*,1001)
               WRITE(*,1002)TPCLC
               IF(FSTOP) STOP
            ENDIF
         ELSEIF(READAUX(1:10).EQ.'*DIMENSION')THEN
            READ(1,'(A)',END=999,ERR=999)DIM
            IF(DIM(1:3).EQ.'.1D')THEN
               READ(1,*,END=999,ERR=999)NP(1)
            ELSEIF(DIM(1:3).EQ.'.2D')THEN
               READ(1,*,END=999,ERR=999)NP(1),NP(2)
            ELSEIF(DIM(1:3).EQ.'.3D')THEN
               READ(1,*,END=999,ERR=999)NP(1),NP(2),NP(3)
            ELSE
               WRITE(*,1001)
               WRITE(*,1002)READAUX
               IF(FSTOP) STOP
            ENDIF
         ELSEIF(READAUX(1:9).EQ.'*GRID_CUT')THEN
            READ(1,*,END=999,ERR=999)(VAR(I), I=1,MXCST,1)
         ELSEIF(READAUX(1:8).EQ.'*TPTRANS')THEN
            READ(1,'(A)',END=999,ERR=999)TPTRANS
         ELSEIF(READAUX(1:8).EQ.'*NREORT')THEN
            READ(1,*,END=999,ERR=999)NREORT
         ELSEIF(READAUX(1:18).EQ.'*INITIAL_POTENTIAL')THEN 
            READ(1,'(A)',END=999,ERR=999)POTCH
            IF(POTCH(1:5).EQ.'.FILE')THEN
               READ(1,'(A)',END=999,ERR=999)POTFILE
            ELSE
               READ(1,*,END=999,ERR=999)(CSTI(I), I=1,MXCST,1)
            ENDIF 
         ELSEIF(READAUX(1:16).EQ.'*FINAL_POTENTIAL')THEN 
            READ(1,'(A)',END=999,ERR=999)POTCHF
            IF(POTCHF(1:5).EQ.'.FILE')THEN
               READ(1,'(A)',END=999,ERR=999)POTFILEF
            ELSE
               READ(1,*,END=999,ERR=999)(CSTF(I), I=1,MXCST,1)

            ENDIF 
         ELSEIF(READAUX(1:10).EQ.'*POTENTIAL')THEN 
            READ(1,'(A)',END=999,ERR=999)POTCH
            POTCHF = POTCH
            IF(POTCH(1:5).EQ.'.FILE')THEN
               READ(1,'(A)',END=999,ERR=999)POTFILE
               POTFILEF = POTFILE
            ELSE
               READ(1,*,END=999,ERR=999)(CSTI(I), I=1,MXCST,1)
               READ(1,*,END=999,ERR=999)(CSTF(I), I=1,MXCST,1)
c               IF(TPPROPG(1:7).EQ.'.P2PSOD')THEN
c                  READ(1,*,END=999,ERR=999)(CSTFB(I), I=1,MXCST,1)
c               ENDIF
c               READ(1,*,END=999,ERR=999)(XI(I), XF(I), I=1,MXDM,1)
            ENDIF 
         ELSEIF(READAUX(1:12).EQ.'*GRID_RANGES')THEN
            READ(1,*,END=999,ERR=999)(XI(I), XF(I), I=1,MXDM,1)
            IF(READAUX(14:15).EQ.'au')THEN 
               DO I=1,MXDM,1
                  XI(I) = XI(I)*A0A
                  XF(I) = XF(I)*A0A
c                  write(*,*)xi(i),xf(i)
               ENDDO   
            ENDIF
         ELSEIF(READAUX(1:9).EQ.'*INIEIGVC')THEN
            READ(1,'(A)',END=999,ERR=999)INIEIGVC
         ELSEIF(READAUX(1:12).EQ.'*MASS_FACTOR')THEN
            READ(1,*,END=999,ERR=999)(SMF(I), I=1,MXCST,1)
         ELSEIF(READAUX(1:5).EQ.'*MASS')THEN
            READ(1,*,END=999,ERR=999)(SM(I), I=1,MXCST,1)
c            write(*,*)'I',I,sm(1),sm(2),sm(3)
c            read(*,*)
c            DO J=1,I,1
c               SM(J)=SMAUX(I-J+1)
c            ENDDO

         ELSEIF(READAUX(1:6).EQ.'*SEED')THEN
            READ(1,*,END=999,ERR=999)NSEED
         ELSEIF(READAUX(1:6).EQ.'*NSHOT')THEN
            READ(1,*,END=999,ERR=999)NSHOT
         ELSEIF(READAUX(1:7).EQ.'*CHANGE')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:3).EQ.'.NO' .OR.
     &           READAUX(1:4).EQ.'.OFF')THEN
               CHANGE = .FALSE.
            ELSE
               CHANGE = .TRUE.
            ENDIF
         ELSEIF(READAUX(1:16).EQ.'*SHIFT_POTENTIAL')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:3).EQ.'.NO' .OR.
     &           READAUX(1:4).EQ.'.OFF')THEN
               SHIFTP = .FALSE.
            ELSE
               SHIFTP  = .TRUE.
            ENDIF   
         ELSEIF(READAUX(1:10).EQ.'*PRTEIGVC2')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               PRTEIGVC2 = .TRUE.
            ELSE
               PRTEIGVC2 = .FALSE.
            ENDIF
         ELSEIF(READAUX(1:9).EQ.'*PRTEIGVC')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               PRTEGVC = .TRUE.
            ELSE
               PRTEGVC = .FALSE.
            ENDIF 
         ELSEIF(READAUX(1:7).EQ.'*PRTPOT')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               PRTPOT = .TRUE.
            ELSE
               PRTPOT = .FALSE.
            ENDIF
c            write(*,*)PRTPOT
         ELSEIF(READAUX(1:8).EQ.'*PRTVEFF')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               PRTVEFF = .TRUE.
            ELSE
               PRTVEFF = .FALSE.
            ENDIF
         ELSEIF(READAUX(1:9).EQ.'*PRTPULSE')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               PRTPULS = .TRUE.
            ELSE
               PRTPULS = .FALSE.
            ENDIF
         ELSEIF(READAUX(1:10).EQ.'*PRTDIPOLE')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               PRTDIPL = .TRUE.
            ELSE
               PRTDIPL = .FALSE.
            ENDIF 
         ELSEIF(READAUX(1:10).EQ.'*DONOTSTOP')THEN
            READ(1,'(A)',END=999,ERR=999)READAUX
            IF(READAUX(1:4).EQ.'.YES' .OR.
     &           READAUX(1:3).EQ.'.ON')THEN
               FSTOP = .FALSE.
            ELSE
               FSTOP = .TRUE.
            ENDIF    
         ELSEIF(READAUX(1:7).EQ.'*PRTCRL')THEN
            READ(1,'(A)',END=999,ERR=999)PRTCRL
            IF(PRTCRL(1:5).EQ.'.FULL' .OR. PRTCRL(1:4).EQ.'.ALL')THEN
               PRTCRL='.YES'
            ELSEIF(PRTCRL(1:5).EQ.'.NONE' .OR. PRTCRL(1:5).EQ.'.NULL' 
     &              .OR.  PRTCRL(1:6).EQ.'.EMPTY' !.OR.
     &              )THEN   
               PRTCRL='.NO'
            ENDIF
         ELSE
            WRITE(*,1001)
            WRITE(*,1002)READAUX
            IF(FSTOP) STOP
         ENDIF
c      ELSE 
c         WRITE(*,1001)
c         WRITE(*,1002)READAUX
c         STOP
      ENDIF  
      GOTO 100
c     ..
c     .. TI GROUP 
 200  READ(1,'(A)',END=999,ERR=999)READAUX
c      write(*,*)readaux
      IF(READAUX(1:2).EQ. '**')GOTO 10
      IF(READAUX(1:1).EQ. '*' .AND. READAUX(2:2).NE. '*')THEN
         IF(READAUX(1:18).EQ.'*INITIAL_POTENTIAL' .OR.
     &        READAUX(1:18).EQ.'*INITIAL_CONDITION' .OR.
     &        READAUX(1:11).EQ.'*INITIAL_WP' .OR.
     &        READAUX(1:11).EQ.'*INITIAL_WF' .OR.
     &        READAUX(1:22).EQ.'*INITIAL_WAVE_FUNCTION' .OR.
     &        READAUX(1:14).EQ.'*WAVE_FUNCTION' .OR.
     &        READAUX(1:12).EQ.'*WAVE_PACKET' .OR.
     &        READAUX(1:10).EQ.'*POTENTIAL')THEN 
            READ(1,'(A)',END=999,ERR=999)POTCH
            IF(POTCH(1:5).EQ.'.FILE')THEN
               READ(1,'(A)',END=999,ERR=999)POTFILE
            ELSE
               READ(1,*,END=999,ERR=999)(CSTI(I), I=1,MXCST,1)
            ENDIF 
         ELSEIF(READAUX(1:7).EQ.'*TPDIAG')THEN
            READ(1,'(A)',END=999,ERR=999)TPDIAG
         ELSEIF(READAUX(1:5).EQ.'*NIST')THEN
            READ(1,*,END=999,ERR=999)NIS, IIL
            IIU = IIL + NIS - ONE
         ELSEIF(READAUX(1:8).EQ.'*NFST')THEN
            READ(1,*,END=999,ERR=999)NFS, IFL
            IFU = IFL + NFS - ONE
         ELSEIF(READAUX(1:7).EQ.'*ABSTOL')THEN
            READ(1,*,END=999,ERR=999)ABSTOL
         ELSEIF(READAUX(1:6).EQ.'*SEED')THEN
            READ(1,*,END=999,ERR=999)NSEED 
         ELSEIF(READAUX(1:7).EQ.'*MAXINT')THEN
            READ(1,*,END=999,ERR=999)MAXINT
         ELSEIF(READAUX(1:17).EQ.'*COLLISION_ENERGY')THEN
            READ(1,*,END=999,ERR=999)(rK(I), I=1,MXDM,1)
d            write(*,*)'rK(I)',rK(1),rK(2),rK(3)
d            read(*,*)
            IF(READAUX(19:20).NE.'au')THEN
               DO I=1,MXDM,1
                  rK(I) = rK(I)/FATEEVAU 
               ENDDO   
            ENDIF
         ELSEIF(READAUX(1:12).EQ.'*GRID_RANGES')THEN
            READ(1,*,END=999,ERR=999)(XI(I), XF(I), I=1,MXDM,1)
            IF(READAUX(14:15).NE.'au')THEN
               DO I=1,MXDM,1
                  XI(I) = XI(I)/A0A
                  XF(I) = XF(I)/A0A
               ENDDO   
            ENDIF
         ELSEIF(READAUX(1:15).EQ.'*INIT_POSSITION')THEN
            READ(1,*,END=999,ERR=999)(X0(I), I=1,MXDM,1)
            IF(READAUX(17:18).NE.'au')THEN
               DO I=1,MXDM,1
                  X0(I) = X0(I)/A0A
               ENDDO
            ENDIF
         ELSEIF(READAUX(1:17).EQ.'*SIZE_OF_GAUSSIAN')THEN
            READ(1,*,END=999,ERR=999)(rA(I), I=1,MXDM,1)
D            write(*,*)rA(1),rA(2),rA(3)
            IF(READAUX(19:20).NE.'au')THEN
               DO I=1,MXDM,1
                  rA(I) = rA(I)/A0A
               ENDDO
            ENDIF
         ELSE 
            WRITE(*,1001)
            WRITE(*,1002)READAUX
            IF(FSTOP) STOP
         ENDIF 
c      ELSE 
c         WRITE(*,1001)
c         WRITE(*,1002)READAUX
c         STOP
      ENDIF
      GOTO 200  
c     ..
c     .. TD GROUP 
 300  READ(1,'(A)',END=999,ERR=999)READAUX
c      write(*,*)readaux
 310  IF(READAUX(1:2).EQ. '**')GOTO 10
      IF(READAUX(1:1).EQ. '*' .AND. READAUX(2:2).NE. '*')THEN
         IF(READAUX(1:7).EQ.'*PROPAG')THEN
            READ(1,'(A)',END=999,ERR=999)TPPROPG
            IF(TPPROPG(1:5).EQ.'.PSIL' .OR.
     &           TPPROPG(1:6).EQ.'.PPSIL' .OR.
     &           TPPROPG(1:6).EQ.'.PPLNZ' .OR.
     &           TPPROPG(1:4).EQ.'.PSH' .OR.
     &           TPPROPG(1:5).EQ.'.PPSH')THEN
               READ(1,*,END=999,ERR=999)TI, TF, DT, MP
            ELSEIF(TPPROPG(1:5).EQ.'.PSOD' .OR. 
     &              TPPROPG(1:6).EQ.'.PFSOD' .OR. 
     &              TPPROPG(1:6).EQ.'.PWSOD' .OR. 
     &              TPPROPG(1:6).EQ.'.PPSOD' .OR.
     &              TPPROPG(1:6).EQ.'.PFSPO' .OR.
     &              TPPROPG(1:6).EQ.'.PWSPO')THEN
               READ(1,*,END=999,ERR=999)TI, TF, DT
            ELSEIF(TPPROPG(1:4).EQ.'.P2P')THEN
               READ(1,*,END=999,ERR=999)TI, TF, DT
               
 20            READ(1,'(A)',END=999,ERR=999)READAUX
               IF(READAUX(1:1).EQ.'*')GOTO 310
               IF(READAUX.EQ.'.PULSE')THEN
c TPX = duration of the pulse (fs)
c T0X = maximum may be shifted to t = t0 (fs)
c KLX = shape of the pulse
                  READ(1,*,END=999,ERR=999)TPX, T0X, SKLX
                  KLX = SKLX
                  VFY(1) = ZERO
                  GOTO 20
               ELSEIF(READAUX(1:22).EQ.'.GAMMA_OMEG_TRNSDIPOL')THEN
c GAMMA = Decay rate of the core excited potential
c OMEGA = Xray frequency ()
c DMX   = Transition dipole moment ()
                  READ(1,*,END=999,ERR=999)GAMMA, OMEGA, DMX
                  GAMMA = GAMMA/FATEEVAU
                  OMEGA = OMEGA/SHBAR
                  DMX   = DMX/FATEEVAU
                  VFY(2) = ZERO
                  GOTO 20
               ELSEIF(READAUX(1:11).EQ.'.AA_QC_DELQ')THEN
                  READ(1,*,END=999,ERR=999)AA, QC, DELQ
                  IF(READAUX(12:15).EQ.'(au)')THEN
                     CONTINUE
                  ELSE
                     AA = AA/FATEEVAU
                     QC = QC/A0A
                     DELQ = DELQ/A0A
                  ENDIF
                  VFY(3) = ZERO
                  GOTO 20
c                  WRITE(*,*)'It is necessary specify the pulse ',
c     &                 'parameters!'
c                  STOP
c               ENDIF
c               READ(1,'(A)',END=999,ERR=999)READAUX
c               IF(READAUX(1:22).EQ.'.GAMMA_OMEG_TRNSDIPOL')THEN
c GAMMA = Decay rate of the core excited potential
c OMEGA = Xray frequency ()
c DMX   = Transition dipole moment ()
c                  READ(1,*,END=999,ERR=999)GAMMA, OMEGA, DMX
c                  GAMMA = GAMMA/FATEEVAU
c                  OMEGA = OMEGA/SHBAR
c                  DMX   = DMX/FATEEVAU
               ELSE
                  WRITE(*,*)'<<<>>> It is necessary specify the decay ra
     &te, frequency, pulse parameters and transition dipole moment! <<<>
     &>>'
                  IF(FSTOP) STOP
               ENDIF  
            ELSE
               WRITE(*,1001)
               WRITE(*,1002)TPPROPG
               IF(FSTOP) STOP
            ENDIF
         ELSEIF(READAUX(1:16).EQ.'*FINAL_POTENTIAL' .OR.
     &           READAUX(1:10).EQ.'*POTENTIAL' )THEN 
            READ(1,'(A)',END=999,ERR=999)POTCHF
            IF(POTCHF(1:5).EQ.'.FILE')THEN
               READ(1,'(A)',END=999,ERR=999)POTFILEF
            ELSE
               READ(1,*,END=999,ERR=999)(CSTF(I), I=1,MXCST,1)
            ENDIF
         ELSEIF(READAUX(1:12).EQ.'*GRID_RANGES')THEN
            READ(1,*,END=999,ERR=999)(XI(I), XF(I), I=1,MXDM,1)
            IF(READAUX(14:15).NE.'au')THEN
               DO I=1,MXDM,1
                  XI(I) = XI(I)/A0A
                  XF(I) = XF(I)/A0A
               ENDDO   
            ENDIF
         ELSEIF(READAUX(1:10).EQ.'*ABSORBING')THEN
            READ(1,'(A)',END=999,ERR=999)TPABSOR
c            write(*,*)TPABSOR
           IF(TPABSOR(1:8).EQ.'.SMOOTHW' .OR. 
     &           TPABSOR(1:5).EQ.'.SMRS')THEN
               ABSORB = .TRUE.
               READ(1,*,END=999,ERR=999)(VOI(I), I=1,MXDM,1)
c               write(*,*)VOI(1),AL(1), AR(1), MXDM
               READ(1,*,END=999,ERR=999)(AL(I), AR(I), I=1,MXDM,1)
c               write(*,*)VOI(1),AL(1), AR(1)
            ELSEIF(TPABSOR(1:7).EQ.'.VOPTIC')THEN
               ABSORB = .TRUE.
               READ(1,*,END=999,ERR=999)(VOI(I), I=1,MXDM,1)
               READ(1,*,END=999,ERR=999)(AL(I), AR(I), I=1,MXDM,1)
            ELSEIF(TPABSOR(1:5).EQ.'.NULL')THEN
               ABSORB = .FALSE.
               READ(1,*,END=999,ERR=999)(WST(I), I=1,MXDM,1)
               READ(1,*,END=999,ERR=999)(WST(I), WST(I), I=1,MXDM,1)
            ELSE
               WRITE(*,1001)
               WRITE(*,1002)TPABSOR
               IF(FSTOP) STOP
            ENDIF
         ELSEIF(READAUX(1:17).EQ.'*PULSE_AND_DIPOLE')THEN
            READ(1,*,END=999,ERR=999)EFC
            IF(EFC(1:5).EQ.'.GAUS')THEN
c E0 = maximum E0 of gaussian envelop (w/cm^2)
c TP = pulse duration (fs)
c OMG = the pulse frequency (eV)
c T0 = maximum may be shifted to t = t0 (fs)
c TD = delay time of the pulse (fs)
c SNI = the phase in general is set iqual to 0
               READ(1,*,END=999,ERR=999)E0, TP, OMG, T0, SNI
            ELSEIF(EFC(1:6).EQ.'.GGAUS')THEN
               READ(1,*,END=999,ERR=999)E0, TP, OMG, T0, SNI, SKL
               KL = SKL
            ELSEIF(EFC(1:5).EQ.'.SIN2')THEN
               READ(1,*,END=999,ERR=999)E0, TP, OMG, TD, SNI
            ELSEIF(EFC(1:5).EQ.'.NONE')THEN
               READ(1,*,END=999,ERR=999)(WST(I), I=1,MXDM,1)
               CONTINUE
            ELSE
               WRITE(*,1001)
               WRITE(*,1002)EFC
            ENDIF
            READ(1,*,END=999,ERR=999)TPDIPL
            IF(TPDIPL(1:5).EQ.'.READ' .OR. TPDIPL(1:4).EQ.'.GET')THEN
               READ(1,*,END=999,ERR=999)DMFILE
            ELSEIF(TPDIPL(1:5).EQ.'.NULL')THEN
               CSTDM(1) = ZERO
            ELSE
               READ(1,*,END=999,ERR=999)(CSTDM(I), I=1,MXCST,1)
            ENDIF
         ELSEIF(READAUX(1:10).EQ.'*PRPGSTATE')THEN
            READ(1,*,END=999,ERR=999)NISPG
         ELSEIF(READAUX(1:8).EQ.'*TPTRANS')THEN
            READ(1,'(A)',END=999,ERR=999)TPTRANS
         ELSEIF(READAUX(1:7).EQ.'*PRPTOL')THEN
            READ(1,*,END=999,ERR=999)TOL
         ELSEIF(READAUX(1:6).EQ.'*NSHOT')THEN
            READ(1,*,END=999,ERR=999)NSHOT
         ELSEIF(READAUX(1:13).EQ.'*NPROJECTIONS')THEN
            READ(1,*,END=999,ERR=999)NPR, KP
         ELSEIF(READAUX(1:8).EQ.'*FOURIER')THEN
            READ(1,*,END=999,ERR=999)MF
         ELSE
            WRITE(*,1001)
            WRITE(*,1002)READAUX
            IF(FSTOP) STOP
         ENDIF 
c      ELSE
c         WRITE(*,1001)
c         WRITE(*,1002)READAUX
c         STOP
      ENDIF
      GOTO 300
c     ..
c     .. SPECTRUM GROUP 
 400  READ(1,'(A)',END=999,ERR=999)READAUX
      IF(READAUX(1:2).EQ. '**')GOTO 10
      IF(READAUX(1:1).EQ. '*' .AND. READAUX(2:2).NE. '*')THEN
         IF(READAUX(1:8).EQ.'*FOURIER')THEN
            READ(1,*,END=999,ERR=999)MF
         ELSEIF(READAUX(1:10).EQ.'*WINDOWING')THEN
            READ(1,'(A)',END=999,ERR=999)TPWIND 
            READ(1,*,END=999,ERR=999)WP
         ELSEIF(READAUX(1:5).EQ.'*HWHM')THEN
c     CHWHM (eV)
c     DE (eV)
c     E1I (eV)
            READ(1,*,END=999,ERR=999)CHWHM, DE, E1I, ADE
            LHWHM = .TRUE.
         ELSE
            WRITE(*,1001)
            WRITE(*,1002)READAUX
            IF(FSTOP) STOP
         ENDIF 
c      ELSE
c         WRITE(*,1001)
c         WRITE(*,1002)READAUX
c         STOP
      ENDIF
      GOTO 400
c
 991  CONTINUE
c
      IF(DIM(1:3).EQ.'.1D')THEN
         DO I=2,MXDM,1
            NP(I) = ONE
         ENDDO
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
         DO I=3,MXDM,1
            NP(I) = ONE
         ENDDO
      ELSEIF(TPCLC(1:9).EQ.'.ONLYSPEC')THEN
         CONTINUE
      ELSE
         WRITE(*,1001)
         WRITE(*,1002)DIM
         IF(FSTOP) STOP
      ENDIF
c
      IF(INIEIGVC(1:4).EQ.'.GET' .OR. INIEIGVC(1:5).EQ.'.READ'
     &     .OR. TPCLC(1:9).EQ.'.ONLYSPEC')THEN
         NIS = 1
      ELSEIF(INIEIGVC(1:5).EQ.'.CALC')THEN
         CONTINUE
      ELSE
         WRITE(*,1001)
         WRITE(*,1002)INIEIGVC
         IF(FSTOP) STOP
      ENDIF
c 
      IF(NPR+KP.GT.NIS)NPR = NPR - NIS
      IF(NPR.GT.6)NPR = 6
c
      IF(TPPROPG(1:7).EQ.'.P2PABM' .AND. VFY(1) .NE. ZERO .AND.
     &     TPPROPG(1:8).NE.'.P2PABM2')THEN
         WRITE(*,*)'It is necessary specify the pulse parameters!'
         IF(FSTOP) STOP
      ELSEIF(TPPROPG(1:7).EQ.'.P2PABM' .AND. VFY(2) .NE. ZERO .AND.
     &        TPPROPG(1:8).NE.'.P2PABM2')THEN
         WRITE(*,*)'<<<>>> It is necessary specify the decay rate, frequ
     &ency and transition the dipole moment! <<<>>>'
         IF(FSTOP) STOP
      ENDIF
c
      IF(TPPROPG(1:7).EQ.'.P2PSOD' .AND. VFY(3) .NE. ZERO .OR.
     &     TPPROPG(1:8).EQ.'.P2PABM2' .AND. VFY(3) .NE. ZERO)THEN
         WRITE(*,*)'<<<>>> It is necessary specify the coupling constant
     &s! <<<>>>'
         IF(FSTOP) STOP
      ENDIF
c
      CLOSE(1)
      WRITE(*,992)
      RETURN
c
 999  WRITE(*,1001)
 2    FORMAT(A72)
 992  FORMAT(1X,'End of input file. Reading finish!')
 1001 FORMAT(1X,'<<<>>> Input file error, check and resubmit! <<<>>>') 
 1002 FORMAT(1X,'It wasn`t possible find',1X,A20,'.')
      IF(FSTOP)THEN
         STOP
      ELSE
         GOTO 3
      ENDIF
      END
