      PROGRAM FLUXP
      IMPLICIT NONE
*     ..
*     Purpose
*     =======
*     To compute the probability of a reaction through the flux expression.
*
*     ..
*     Argumens
*     =========
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes - email:freddy@ufmg.br
*
*     ..
*     Historic
*     ========
*     (01/10/2004) First version of FLUXP written by Freddy 
*     This first version the program is based in the NORMF also 
*     written by Freddy 
*     
c     ** 
c     ** Parameters 
      INTEGER       NDIM, MXDCT, MXSPL, LMTREORT, KXORD, KYORD
      PARAMETER     (NDIM = 3.0D+0, 
     &     MXDCT = 1.0D+5,
     &     MXSPL = 3.5D+2,
     &     LMTREORT = 2.0D+1,
     &     KXORD = 3.0D+0,
     &     KYORD = 3.0D+0)
      REAL*8        ZERO, ONE, TWO, AUFS
      PARAMETER     (ZERO = 0.0D+0, 
     &     ONE = +1.0D+0, 
     &     TWO = +2.0D+0,
     &     AUFS = +41.34145D+0)
c     **
c     ** Scalar arguments
      LOGICAL       PRTPOT, LEVELP, PRTEGVC
      CHARACTER*4   CHNUM 
      CHARACTER*16  TPDIAG, EIGST
      CHARACTER*72  INFILEAUX, INPESFILE, INFILE, CHLX, RDAUX, INFLEIG
      INTEGER       NC, NINI, NFIN, ND, IOTEST, I, J, K, L, NT, NPP
      INTEGER       IIL, IIU, INFO, IFLTH, IFLTHT, NPTU, MP1, LAUX
      INTEGER       NPF, KM
      REAL*8        DT, DT1, DT2, RXA, RYA, PROB1, PROB2, PROB3, Q1A 
      REAL*8        Q2A, RX, RY, ABSTOL, Q01, Q02, XAI, XAF, BETA, XSS
      REAL*8        UR, VI, FLUXPT, Q1LIM, X, Y, X0, Y0, SSM, YAI, YAF
      REAL*8        YSS, FLUXPTI, DUR, DVI, DY, DUR1, DVI1, DUR2, DVI2
      REAL*8        tot1, tot2, DUR11, DUR01, DUR10, DUR00, DVI11, DVI01 
      REAL*8        DVI10, DVI00, ABSTOL2
c     **
c     ** Array arguments
      LOGICAL       CHSTEP(NDIM)
      INTEGER       NP(NDIM), IWK(MXDCT), NPE(NDIM)
      REAL*8        XI(NDIM), XF(NDIM), SH(NDIM), U(MXDCT), V(MXDCT)
      REAL*8        EIGVL(LMTREORT), SHM(NDIM), Q1(MXDCT), Q2(MXDCT)
      REAL*8        SM(NDIM), SHY(NDIM), FLUXPP(LMTREORT), WK1(MXDCT)
      REAL*8        WK2(MXDCT), XKNOT(MXSPL), YKNOT(MXSPL), XP(NDIM)
      REAL*8        SR(LMTREORT), SI(LMTREORT), SDR(LMTREORT), XIE(NDIM)
      REAL*8        SDI(LMTREORT), FLUXPPI(LMTREORT), SHE(NDIM) 
      REAL*8        BCOEFR(MXSPL,MXSPL), BCOEFI(MXSPL,MXSPL)
      REAL*8        GST(MXDCT,LMTREORT), FQI(MXSPL,MXSPL)
      REAL*8        FQR(MXSPL,MXSPL)
c     **
c     ** External functions 
      INTEGER       ICHLENGTH
      REAL*8        DBS2VL 
c     **
c     ** External subroutines 
      EXTERNAL      DBS2IN, DBSNAK, RDPT2
c     ..
c     .. Starting program
      WRITE(*,1001)
      WRITE(*,1002)
      WRITE(*,1001)
      WRITE(*,1003)'Initializing...'
c     .. 
c     .. Start logical variables 
      PRTPOT = .FALSE.
      LEVELP = .FALSE.
      PRTEGVC = .FALSE.
c     .. 
c     .. Start character variables
      CHNUM     = '1   ' 
      INFILEAUX = 'ReIm_   '
      INPESFILE = 'veff_0001.dat  '
      INFILE    = 'Null    '
      CHLX      = 'Null    '
      EIGST     = '.CALC   '
      TPDIAG    = '.MTRXDIAG   '
c     .. 
c     .. Start integer variables
      ND   = ZERO
      NT   = ZERO
      INFO = ZERO
      IIL  = ONE
      IIU  = ONE
      NPP  = ONE
      NINI = ONE
      ND   = 2.0D+0
      NFIN = 1.0D+4
c     .. 
c     .. Start real variables 
      DT  = ZERO 
      DT1 = ZERO
      DT2 = ZERO
      RXA = ZERO
      RYA = ZERO
      RX  = ZERO
      RY  = ZERO
      DY  = 1.0D-4
      ABSTOL  = 1.0D-12
d      ABSTOL2 = 1.0D-4
c     ..
c     .. Start logical arrays
      DATA CHSTEP  / NDIM* .TRUE. /
c     ..
c     .. Starting integer arrays 
      DATA NP / NDIM* ONE /
c     .. 
c     .. Starting real arrays
      DATA U  / MXDCT* ZERO /
      DATA V  / MXDCT* ZERO /
      DATA XI / NDIM* ZERO /
      DATA XF / NDIM* ZERO /
      DATA SH / NDIM* ZERO /
      DATA Q1 / MXDCT* ZERO /
      DATA Q2 / MXDCT* ZERO /
      DATA EIGVL / LMTREORT* ZERO /
c 
c     .. Getting input parameters
      WRITE(*,1003)'Getting input parameters...'
 1    READ(*,*)RDAUX 
d      print *,rdaux
      IF(RDAUX(1:16).EQ.'*INPUT_FILE_NAME')THEN
         READ(*,*)INFILEAUX
         IFLTH =  ICHLENGTH(INFILEAUX, 1) 
         IFLTHT = IFLTH + 8
      ELSEIF(RDAUX(1:14).EQ.'*PES_FILE_NAME')THEN
         READ(*,*)INPESFILE
      ELSEIF(RDAUX(1:22).EQ.'*NUMBER_OF_POINTS_FLUX')THEN
         READ(*,*)NPTU
      ELSEIF(RDAUX(1:21).EQ.'*NUMBER_OF_DIMENSIONS')THEN
         READ(*,*)ND
      ELSEIF(RDAUX(1:20).EQ.'*INITIAL_FILE_NUMBER')THEN
         READ(*,*)NINI, NFIN
      ELSEIF(RDAUX(1:23).EQ.'*STATIONARY_EIGENSTATES')THEN
         READ(*,*)EIGST
         print *,'eigst',eigst
      ELSEIF(RDAUX(1:14).EQ.'*FLUX_POSITION')THEN
         READ(*,*)Q01,Q02          
      ELSEIF(RDAUX(1:16).EQ.'*SPLINE_STARTING')THEN
         READ(*,*)Q1LIM
      ELSEIF(RDAUX(1:6).EQ.'*XI_XF')THEN
         READ(*,*)XAI,XAF   
      ELSEIF(RDAUX(1:3).EQ.'*DY')THEN
         READ(*,*)DY   
      ELSEIF(RDAUX(1:12).EQ.'*BETA_FACTOR')THEN
         READ(*,*)BETA
      ELSEIF(RDAUX(1:12).EQ.'*PROB_PROJEC')THEN
         READ(*,*)RDAUX
         IF(RDAUX(1:3).EQ.'.ON' .OR. RDAUX(1:4).EQ.'.YES')THEN
            LEVELP  = .TRUE.
         ELSE
            LEVELP  = .FALSE.
         ENDIF 
      ELSEIF(RDAUX(1:16).EQ.'*PRINT_POTENTIAL')THEN
         READ(*,*)RDAUX
         IF(RDAUX(1:3).EQ.'.ON' .OR. RDAUX(1:4).EQ.'.YES')THEN
            PRTPOT = .TRUE.
         ELSE
            PRTPOT = .FALSE.
         ENDIF
      ELSEIF(RDAUX(1:18).EQ.'*PRINT_EIGENVECTOR')THEN
         READ(*,*)RDAUX
         IF(RDAUX(1:3).EQ.'.ON' .OR. RDAUX(1:4).EQ.'.YES')THEN
            PRTEGVC = .TRUE.
         ELSE
            PRTEGVC = .FALSE.
         ENDIF
      ELSEIF(RDAUX(1:8).EQ.'*IIL_IIU')THEN
         READ(*,*)IIL, IIU 
         IIL = IIL + ONE
         IIU = IIU + ONE
      ELSEIF(RDAUX(1:6).EQ.'*MMASS')THEN
         READ(*,*)SSM 
      ELSEIF(RDAUX(1:5).EQ.'*MASS')THEN
         READ(*,*)(SM(I), I=1,NDIM,1)
      ELSEIF(RDAUX(1:4).EQ.'*END')THEN
         GOTO 2
      ELSE 
         WRITE(*,*)'<<<>>> Keyword,', RDAUX, 
     &        ' was not recognized. <<<>>>'
         WRITE(*,*)'Please, check the imput file and resubmit!'
         STOP
      ENDIF
      GOTO 1
 2    CONTINUE
      NPP = IIU - IIL + ONE
      WRITE(*,*)'Input parameters got.' 
c     ..
c     .. Write input description
      WRITE(*,*)
      WRITE(*,*)'Input description:'
      WRITE(*,*)'""""""""""""""""""'
      WRITE(*,*)'   The wave packet input file name is: ',
     &     INFILEAUX(1:IFLTH)//'????.dat'
c     ..
c     .. Read potential if asked and check the size of the files
      IF(EIGST(1:5).EQ.'.CALC')THEN
         NC = ICHLENGTH(INPESFILE, 0)
         OPEN(2, STATUS='OLD', FILE=INPESFILE(1:NC), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)THEN
            WRITE(*,*)'<<<>>> Potential file "', INFILE(1:NC), 
     &           '" cannot be opened <<<>>>'
            STOP
         ENDIF 
         READ(2,*,END=999,ERR=999)CHLX
         READ(2,*, END=10,ERR=10)(XI(J), J=ND,1,-1), U(1)
         WRITE(*,*)'Initial grid points: ',('D',J,' = ',XI(J),';   ',
     &        J=ND,1,-1)
c     
         DO I=2,MXDCT,1
            READ(2,*, END=10,ERR=10)(XF(J), J=ND,1,-1), U(I)
d            print *,i,xf(1),xf(2),u(i)
c     .. Get the step in space
            DO J=1,ND,1
               IF(CHSTEP(J) .AND. XF(J).NE.XI(J))THEN
                  CHSTEP(J) = .FALSE.
                  SH(J) = XF(J) - XI(J)
                  NP(J) = ONE
               ELSEIF(CHSTEP(J))THEN
                  NP(J-ONE) = NP(J-ONE) + ONE
               ENDIF
            ENDDO
         ENDDO
      ENDIF
c     ..
c     .. Check the size of the files when the potential is not 
c     given.
      IF(EIGST(1:5).NE.'.CALC')THEN
         NC = ICHLENGTH(INFILEAUX, 0)
         WRITE(CHNUM,'(I4)')NINI
         IF(NINI.LT.10)THEN
            INFILE = INFILEAUX(1:IFLTH)//'000'//CHNUM(4:4)//'.dat'
         ELSEIF(NINI.LT.100)THEN
            INFILE = INFILEAUX(1:IFLTH)//'00'//CHNUM(3:4)//'.dat'
         ELSEIF(NINI.LT.1000)THEN
            INFILE = INFILEAUX(1:IFLTH)//'0'//CHNUM(2:4)//'.dat'
         ELSEIF(NINI.LT.10000)THEN
            INFILE = INFILEAUX(1:IFLTH)//CHNUM(1:4)//'.dat'
         ENDIF
         NC = ICHLENGTH(INFILE, 0)
         OPEN(2, STATUS='OLD', FILE=INFILE(1:NC), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)THEN
            WRITE(*,*)'<<<>>> Input file "', INFILE(1:NC), 
     &           '" cannot be opened <<<>>>'
            STOP
         ENDIF
         READ(2,*,END=999,ERR=999)CHLX, CHLX, DT1
         READ(2,*, END=10,ERR=10)(XI(J), J=ND,1,-1), U(1), V(1)
         WRITE(*,*)'Initial grid points:',('D_',J,' = ',XI(J),';   ',
     &        J=ND,1,-1)
c     
         DO I=2,MXDCT,1
            READ(2,*, END=10,ERR=10)(XF(J), J=ND,1,-1), U(I), V(I)
c     .. Get the step in space
            DO J=1,ND,1
               IF(CHSTEP(J) .AND. XF(J).NE.XI(J))THEN
                  CHSTEP(J) = .FALSE.
                  SH(J) = XF(J) - XI(J)
                  NP(J) = ONE
               ELSEIF(CHSTEP(J))THEN
                  NP(J-ONE) = NP(J-ONE) + ONE
               ENDIF
            ENDDO
         ENDDO
      ENDIF
c     ..
c     .. Exit from the loop indicating the end of the file.
 10   CONTINUE
c     ..
c     .. Compute the total number of points and number of points per dimension
      NT = I - ONE
      NP(ND) = NT
      DO I=ND-ONE,1,-1
         NP(ND) = NP(ND)/NP(I)
      ENDDO 
c
      WRITE(*,*)'Final grid points: ',('D',J,' = ',XF(J),';   ',
     &     J=ND,1,-1)
      WRITE(*,*)'Step between points in the grid:'
     &     ,('D',J,'= ',SH(J),';   ',J=ND,1,-1)
      DO I=1,ND,1
         NP(I) = NINT((XF(I) - XI(I))/SH(I) + ONE)
      ENDDO
      WRITE(*,*)'There are ',NT,' points in each file.'
      WRITE(*,*)'   with ',(NP(J),' in the dimension', J,';   ',
     &     J=ND,1,-1)
c     .. Compute step in the direction pependicular to the reaction path 
c      XAI = XI(I) - BETA*(-XI(I)*BETA + Q01*BETA + Q02)
c      XAF = (XI(I) + (NP(1) - ONE)*SH(1))*(ONE + BETA**2) 
c     &     - Q01*BETA**2 - Q02*BETA
      X0 = Q01 - BETA*Q02
      Y0 = 5.0D-1*(Q01 + ONE/BETA*Q02)
d      print *, 'X0,Y0',X0,Y0
d      XAI = TWO*(XI(1) -   Y0)
d      IF(TWO*Y0 - TWO*BETA*XI(2).GT.XAI)XAI = TWO*(Y0 - BETA*XI(2))  
d      XAF = TWO*(XF(1) -   Y0)
d      IF(TWO*Y0 - TWO*BETA*XF(2).LT.XAF)XAF = TWO*(Y0 - BETA*XF(2))
d      print *,  'XAI, XAF',TWO*(XI(1) -   Y0),TWO*(Y0 - BETA*XI(2)),
d     &     TWO*(XF(1) -   Y0), TWO*(Y0 - BETA*XF(2)),XAI, XAF
c
      XAI = X0 - XAI
      XAF = X0 + XAF
      XSS = (XAF - XAI)/(NPTU - ONE)
c
c      YAI = Y0 - 5.58
c      YAF = Y0 + 3.1
d      YSS = (YAF - YAI)/(NPTU - ONE)
c
d      Q1LIM = 0d0!5.0D-1*X0 + Y0
      WRITE(*,*)'X0 = ',X0,';  Y0 = ',Y0
      WRITE(*,1010)1,XAI,1,XAF,1,XSS
      IF(XSS.GT.4.0D-2)THEN
         WRITE(*,*)'<<<>>> Step bigger than', 4.0D-2,
     &        ' a.u. <<<>>>'
         STOP
      ELSEIF(XSS.GT.1.0D-2)THEN
         WRITE(*,*)'-> Warning step bigger than', 1.0D-2,
     &        ' a.u. !'
      ENDIF  
c     ..
c     .. Generate Q vectors
      KM = ONE
      DO I=1,NP(1),1
         Q1(I) = XI(1) + (I - ONE)*SH(1)
         IF(Q1(I).LT.Q1LIM)KM = KM + 1
      ENDDO
      DO I=1,NP(2),1
         Q2(I) = XI(2) + (I - ONE)*SH(2)
      ENDDO
c     ..
c     .. Compute K_0s referent to q_1^0 and q_2^0 positions
d      K01 = NINT((Q01 - XI(1))/SH(1) + ONE)
d      K02 = NINT((Q02 - XI(2))/SH(2) + ONE)
c     ..
c     .. Generate KNOTs for the spline procedure
      MP1 = NP(1) - KM + ONE
      CALL DBSNAK(MP1, Q1(KM), KXORD, XKNOT)
      CALL DBSNAK(NP(2), Q2(1), KYORD, YKNOT)
c     ..
c     .. Reading stationary eigenstates 
      IF(LEVELP .AND. EIGST(1:5).NE.'.CALC')THEN
         WRITE(*,1003)'Getting eigenvetors ...'
         NC = ICHLENGTH(INFLEIG, 0)
         OPEN(11, STATUS='OLD', FILE=INFLEIG(1:NC), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)THEN
            WRITE(*,*)'<<<>>> Input file "', INFILE(1:NC), 
     &           '" cannot be opened <<<>>>'
            STOP
         ENDIF 
         DO J=1,MXDCT,1
            READ(11,*, END=20,ERR=20)(GST(J,K), K=1,NPP,1)
         ENDDO
 20      CONTINUE
         NPF = J - ONE
         IF(NPF.NE.NPTU)THEN
            WRITE(*,*)'<<<>>> Number of data in the file differ from the
     &inputed one <<<>>>'
            STOP
         ENDIF
      ENDIF
      IF(EIGST(1:5).EQ.'.CALC')THEN
c     .. Loop to collect the data for the spline
         DO J=1,NP(2),1
            LAUX = (J - ONE)*NP(1)
            DO K=KM,NP(1),1
               L = K + LAUX
               FQR(K,J) = U(L)
            ENDDO
         ENDDO
c     .. Compute the matrix of coefficients for the spline
         CALL DBS2IN(MP1, Q1(KM), NP(2), Q2(1), FQR(KM,1), MXSPL, KXORD, 
     &        KYORD, XKNOT, YKNOT, BCOEFR)
c     .. Collect data in the direction pependicular to the reaction path at y0
         IF(PRTPOT)WRITE(*,*)'Potential points:'
         IF(PRTPOT)WRITE(*,*)'"""""""""""""""""'
         DO I=1,NPTU,1
            X = XAI + (I - ONE)*XSS 
            Y = Y0
d            Y = YAI + (I - ONE)*YSS
d            X = X0
            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
            Q1A = X + BETA*Q2A
c     .. Generate througth spline the pontential in the direction 
c     pependicular to the reaction path at y0            
            V(I) = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFR)
c     .. Print the potential points when asked
            IF(PRTPOT)WRITE(*,*)X, V(I)
            IF(PRTPOT)WRITE(3,*)Q2A, Q1A, V(I)
         ENDDO
c     .. Compute eigenvalues and eigenstates in the direction pependicular 
c     to the reaction path
         WRITE(*,*)
         WRITE(*,*)'Computing the eigenvalues and eigenvectors of the pote
     &ntial,'
         IF(TPDIAG(1:9).EQ.'.MTRXDIAG')THEN
            WRITE(*,*)'Using complete matrix diagonalization procedure,'
            SHM(1) = ONE/(2.4D+1*SSM*1.836152667539D+3*XSS*XSS)
            print *, shm(1),SSM,XSS,iil, iiu
            CALL MTRXDIAG('.1D', IIL, IIU, INFO, MXDCT, NPTU, ABSTOL,  
     &            IWK, NP, EIGVL, SHM, V, WK1, WK2, GST)
D*12         ELSEIF(TPDIAG(1:8).EQ.'.LANCZSG')THEN
D*12            WRITE(*,*)'Using interative lanczos tri-diagonalization ',
D*12     &           'procedure,'
D*12            CALL LANCZSG(CHANGE, DIM, iil, iiu, INFO, NREORT, LMTREORT, 
D*12     &           NT, NSEED, MXDCT, ABSTOL, IWK, MAXINT, np, U1, V1, 
D*12     &           EIGVL, SHM, WK1, WK2, v, WK3, WK4, gst, LANCZ)
D*12         ELSEIF(TPDIAG(1:7).EQ.'.LANCZS')THEN
D*12            WRITE(*,*)'Using lanczos tri-diagonalization procedure,'
D*12            CALL LANCZS(CHANGE, DIM, iil, iiu, INFO, NREORT, LMTREORT, 
D*12     &           nt, NSEED, MXDCT, ABSTOL, IWK, np, U1, V1, EIGVL, 
D*12     &           SHM, WK1, WK2, v, WK3, WK4, gst, LANCZ)
         ELSE
            WRITE(*,*)'<<<>>> Initial eigenvalue(s) and eigenvector(s)', 
     &           'error <<<>>>'
            WRITE(*,*)TPDIAG
            STOP
         ENDIF
      ENDIF
c     .. Normalize the eigenvectors
      DO I=1,NPP,1
         DO J=1,NPTU,1
           GST(J,I) = GST(J,I)/SQRT(XSS)
d           IF(ABS(GST(J,I)).LT.ABSTOL2)GST(J,I) = ZERO
         ENDDO
      ENDDO
c
      IF(INFO.NE.0)THEN 
         WRITE(*,*)'<<<>>> It wasn`t possible obtain eignvalues and ',
     &        'eignvectors correctly <<<>>>'
         WRITE(*,*)'Exiting.'
         STOP
      ENDIF
      IF(EIGST(1:5).EQ.'.CALC')THEN
         WRITE(*,*)
         WRITE(*,*)'Eigenvalues:'
         WRITE(*,*)'""""""""""""'
         WRITE(*,1014)XSS
         WRITE(*,*)'|  Eigenstate  |  Eigenvalues   |'
         WRITE(*,*)'================================='
         DO I = IIL,IIU,1
            WRITE(*,1015)I-1,EIGVL(I-IIL+1)
            WRITE(*,1014)
         ENDDO
      ENDIF
      IF(PRTEGVC)THEN
         NPE(1) = NPTU
         NPE(2) = ONE
         NPE(3) = ONE
         XIE(1) = XAI
         SHE(1) = XSS
         WRITE(*,*)
         WRITE(*,*)'Eigenvectors:'
         WRITE(*,*)'"""""""""""""'
         CALL PRTEIGVC(TPDIAG, 6, MXDCT, 1, NPP, NPE, XP, XIE, SHE, 
     &        GST)
      ENDIF
c
c     .. Check the time step of the files
      IF(LEVELP .AND. EIGST(1:5).EQ.'.CALC')THEN
         NC = ICHLENGTH(INFILEAUX, 0)
         WRITE(CHNUM,'(I4)')NINI
c         print *,'chnum',nini,'   ',chnum
         IF(NINI.LT.10)THEN
            INFILE = INFILEAUX(1:NC)//'000'//CHNUM(4:4)//'.dat'
         ELSEIF(NINI.LT.100)THEN
            INFILE = INFILEAUX(1:NC)//'00'//CHNUM(3:4)//'.dat'
         ELSEIF(NINI.LT.1000)THEN
            INFILE = INFILEAUX(1:NC)//'0'//CHNUM(2:4)//'.dat'
         ELSEIF(NINI.LT.10000)THEN
            INFILE = INFILEAUX(1:NC)//CHNUM(1:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF
         NC = ICHLENGTH(INFILE, 0)
         OPEN(2, STATUS='OLD', FILE=INFILE(1:NC), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)THEN
            WRITE(*,*)'<<<>>> Input file "', INFILE(1:NC), 
     &           '" cannot be opened <<<>>>1'
            STOP
         ENDIF
         READ(2,*,END=999,ERR=999)CHLX, CHLX, DT1
         CLOSE(2)
      ENDIF
      NC = ICHLENGTH(INFILEAUX, 0)
      WRITE(CHNUM,'(I4)')NINI + 1
      IF(NINI + 1.LT.10)THEN
         INFILE = INFILEAUX(1:NC)//'000'//CHNUM(4:4)//'.dat'
      ELSEIF(NINI + 1.LT.100)THEN
         INFILE = INFILEAUX(1:NC)//'00'//CHNUM(3:4)//'.dat'
      ELSEIF(NINI + 1.LT.1000)THEN
         INFILE = INFILEAUX(1:NC)//'0'//CHNUM(2:4)//'.dat'
      ELSEIF(NINI + 1.LT.10000)THEN
         INFILE = INFILEAUX(1:NC)//CHNUM(1:4)//'.dat'
      ELSE
         WRITE(*,*)'Too many files to be read. Stopping!'
         STOP
      ENDIF
      NC = ICHLENGTH(INFILE, 0)
      OPEN(2, STATUS='OLD', FILE=INFILE(1:NC), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:NC), 
     &        '" cannot be opened <<<>>>2'
         STOP
      ENDIF
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT2
      DT = ABS(DT2 - DT1)
      CLOSE(2)
      WRITE(*,*)'The time step between files is:',DT,'.'
c
      WRITE(*,1011)
      FLUXPTI = ZERO
      DO K=1,NPP,1
         FLUXPPI(K) = ZERO
      ENDDO
c
      DO I=NINI,NFIN,1
d         print *, 'i',i
c     ..
c     .. Generate the input file name to be read
         NC = ICHLENGTH(INFILEAUX, 0)
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = INFILEAUX(1:NC)//'000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE = INFILEAUX(1:NC)//'00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = INFILEAUX(1:NC)//'0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = INFILEAUX(1:NC)//CHNUM(1:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF
c     ..
c     .. Read data from "INFILE"
         CALL RDPT2(INFILE, NT, ND, XI, XF, U, V)
c     ..
c     .. Loop to collect the data for the spline
         tot1 = zero
         tot2 = zero
         DO J=1,NP(2),1
            LAUX = (J - ONE)*NP(1)
            DO K=KM,NP(1),1
               L = K + LAUX
               FQR(K,J) = U(L)/SQRT(SH(1)*SH(2))
               FQI(K,J) = V(L)/SQRT(SH(1)*SH(2))
d               tot2 = tot2 + (fqr(k,j)*fqr(k,j) + fqi(k,j)*fqi(k,j))
d     &              *sh(1)*sh(2)
d               tot1 = tot1 + (u(l)*u(l) + v(l)*v(l))*sh(1)*sh(2)
            ENDDO
         ENDDO
d         write(*,*)'tot1,sqrt(sh(1)*sh(2))',tot1,tot2,sh(1)*sh(2)
c     ..
c     .. Perform the spline procedure
         CALL DBS2IN(MP1, Q1(KM), NP(2), Q2(1), FQR(KM,1), MXSPL, KXORD, 
     &        KYORD, XKNOT, YKNOT, BCOEFR)
         CALL DBS2IN(MP1, Q1(KM), NP(2), Q2(1), FQI(KM,1), MXSPL, KXORD, 
     &        KYORD, XKNOT, YKNOT, BCOEFI)

d         x0 = 0
d         y0 = 0
d         xai = x0 - 10
d         xaf = x0 + 10
d         yai = y0 - 10
d         yaf = y0 + 10
d         xss = (xaf - xai)/(nptu - one)
d         yss = (yaf - yai)/(nptu - one)
d         
d         write(*,*)(yknot(k),k=1,nptu+5,1)
d         tot1 = zero
d         do j=1,nptu,1
d            y = yai + (j - one)*yss
d            do k=1,nptu,1
d               x = xai + (k - one)*xss 
d               q2a = beta*(two*y - x)/(one + beta**2)
d               q1a = x + beta*q2a
c               print *,x,y,q1a,q2a,xf(1),xf(2),xi(1),xi(2)
d               if(q1a.ge.xf(1) .or. q2a.ge.xf(2) .or. 
d     &              q1a.le.xi(1) .or. q2a.le.xi(2))goto 1212
d               print *,'oi'
d               ur = dbs2vl(q1a, q2a, kxord, kyord, xknot, yknot, 
d     &              mp1, np(2), bcoefr)
d               vi = dbs2vl(q1a, q2a, kxord, kyord, xknot, yknot, 
d     &           mp1, np(2), bcoefi)
d               write(23,*)y, x, ur, vi
d               tot1 = tot1 + (ur*ur + vi*vi)*yss*xss
d 1212          continue
d            enddo
d         enddo
d         print *,'norm =',tot1
c     ..
c     .. Perform the probability calculation via the flux expression
c         DATA FLUXPP / LMTREORT* ZERO /
c         DATA SR  / LMTREORT* ZERO /
c         DATA SI  / LMTREORT* ZERO /
c         DATA SDR / LMTREORT* ZERO /
c         DATA SDI / LMTREORT* ZERO /
         FLUXPT = ZERO
         DO K=1,NPP,1
            FLUXPP(K) = ZERO
            SR(K) = ZERO  
            SI(K) = ZERO 
            SDR(K) = ZERO 
            SDI(K) = ZERO
         ENDDO
c     .. Loop in the eigenstates to perform the projections
c         print *, 'sr',sr(1),sr(2),sr(3),sr(4)
c         print *, 'si',si(1),si(2),si(3),si(4)
c         print *, 'sdr',sdr(1),sdr(2),sdr(3),sdr(4)
c         print *, 'sdi',sdi(1),sdi(2),sdi(3),sdi(4)
         DO J=1,NPTU,1
c     .. Collect data in the direction pependicular to the reaction path at y0
            X = XAI + (J - ONE)*XSS 
            Y = Y0
            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
            Q1A = X + BETA*Q2A
c     .. Generate through spline the digenvector in the direction 
c     pependicular to the reaction path at y0
            UR = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFR)
            VI = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFI)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
d            IF(ABS(UR).LT.ABSTOL2)UR = ZERO
d            IF(ABS(VI).LT.ABSTOL2)VI = ZERO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            Y = Y0 + DY
            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
            Q1A = X + BETA*Q2A
            DUR1 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFR)
            DVI1 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFI)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            IF(ABS(DUR1).LT.ABSTOL)DUR1 = ZERO
c            IF(ABS(DVI1).LT.ABSTOL)DVI1 = ZERO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            Y = Y0 - DY
            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
            Q1A = X + BETA*Q2A
            DUR2 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFR)
            DVI2 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
     &           MP1, NP(2), BCOEFI)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            IF(ABS(DUR2).LT.ABSTOL)DUR2 = ZERO
c            IF(ABS(DVI2).LT.ABSTOL)DVI2 = ZERO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            DUR = (DUR1 - DUR2)/(TWO*DY)
            DVI = (DVI1 - DVI2)/(TWO*DY)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
d            IF(ABS(DUR).LT.ABSTOL2)DUR = ZERO
d            IF(ABS(DVI).LT.ABSTOL2)DVI = ZERO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            X = XAI + (J - ONE)*XSS + DY
c            Y = Y0 + DY
c            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
c            Q1A = X + BETA*Q2A
c            DUR11 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFR)
c            DVI11 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFI)
c
c            X = XAI + (J - ONE)*XSS - DY
c            Y = Y0 + DY
c            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
c            Q1A = X + BETA*Q2A
c            DUR10 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFR)
c            DVI10 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFI)
c
c            X = XAI + (J - ONE)*XSS + DY
c            Y = Y0 - DY
c            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
c            Q1A = X + BETA*Q2A
c            DUR01 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFR)
c            DVI01 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFI)
c
c            X = XAI + (J - ONE)*XSS - DY
c            Y = Y0 - DY
c            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
c            Q1A = X + BETA*Q2A
c            DUR00 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFR)
c            DVI00 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
c     &           MP1, NP(2), BCOEFI)
c
c            DUR1 = (DUR11 - DUR01 + DUR10 - DUR00)/(4.0D0*DY)
c            DVI1 = (DVI11 - DVI01 + DVI10 - DVI00)/(4.0D0*DY)
c            X = XAI + (J - ONE)*XSS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
d            Y = Y0 + TWO*DY
d            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
d            Q1A = X + BETA*Q2A
d            DUR11 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
d     &           MP1, NP(2), BCOEFR)
d            DVI11 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
d     &           MP1, NP(2), BCOEFI)
c
d            Y = Y0 - TWO*DY
d            Q2A = BETA*(TWO*Y - X)/(ONE + BETA**2)
d            Q1A = X + BETA*Q2A
d            DUR00 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
d     &           MP1, NP(2), BCOEFR)
d            DVI00 = DBS2VL(Q1A, Q2A, KXORD, KYORD, XKNOT, YKNOT, 
d     &           MP1, NP(2), BCOEFI)
c
d            DUR = (-DUR11 + 8.0D0*DUR1 - 8.0D0*DUR2 + DUR00)
d     &           /(1.2D1*DY)
d            DVI = (-DVI11 + 8.0D0*DVI1 - 8.0D0*DVI2 + DVI00)
d     &           /(1.2D1*DY)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     .. Compute the projections in the estationary wave functions
            IF(LEVELP)THEN
               DO K=1,NPP,1
                  SR(K) = SR(K) + UR*GST(J,K)*XSS
                  SI(K) = SI(K) + VI*GST(J,K)*XSS
                  SDR(K) = SDR(K) + DUR*GST(J,K)*XSS
                  SDI(K) = SDI(K) + DVI*GST(J,K)*XSS
d                  FLUXPP(K) =  FLUXPP(K) +
d     &                (UR*GST(J,K)*DVI*GST(J,K) 
d     &                 - VI*GST(J,K)*DUR*GST(J,K))*XSS*XSS 
               ENDDO
            ENDIF 
c     .. 
d            FLUXPT = FLUXPT + (UR*UR + VI*VI)*XSS
            FLUXPT = FLUXPT + UR*DVI - VI*DUR
c            write(13,1030)x,ur,vi,dur,dvi,dur1,dvi1!,dvi-dvi1!(gst(j,k), k=1,npp,1)
c            WRITE(12,*)X, (UR*UR + VI*VI)*XSS, UR, VI
c            write(11,*)q2a,q1a,ur**2+vi**2
         ENDDO
         write(25,*)dt1 + (i - 1)*dt, (ur*ur + vi*vi)*xss,
     &        (sr(k)*sr(k) + si(k)*si(k), k=1,npp,1)
         FLUXPT = FLUXPT*XSS
c         close(13)
c         print *,'oi'
c         stop
c         read(*,*)
c         write(12,*)
c         close(11)
         IF(LEVELP)THEN
            DO K=1,NPP,1
               FLUXPP(K) = (SR(K)*SDI(K) - SI(K)*SDR(K))!*XSS*XSS
               FLUXPPI(K) = FLUXPPI(K) + FLUXPP(K)*DT*AUFS
            ENDDO
         ENDIF


         FLUXPTI = FLUXPTI + FLUXPT*DT*AUFS
         WRITE(1,*) 
         IF(LEVELP)THEN
            WRITE(*,1030)DT1 + (I - 1)*DT, FLUXPT, FLUXPTI, 
     &           (FLUXPP(K), K=1,NPP,1),(FLUXPPI(K), K=1,NPP,1)
         ELSE
            WRITE(*,*)DT1 + (I - 1)*DT, FLUXPT, FLUXPTI 
         ENDIF
      ENDDO
c
C vinicius 27/09/2011-----------
C      RETURN
C--------------------
c

 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
c     ..
 1011 FORMAT(3X,'#w/a.u.#',5X,'#Int/arb. units#')
 1012 FORMAT(1X,E12.6,5X,E12.6)
 1014 FORMAT(1X,'---------------------------------')
 1015 FORMAT(1X,'|',3X,I3,8X,'|',2X,E12.6,2X,'|')
 1030 FORMAT(1X,F12.6,3X,100(E12.6,4X))
 1010 FORMAT(1X,'X',I1,'_i/a.u. =',1X,E12.6,';',1X,'X',I1,'_f/a.u. =',
     &     1X,E12.6,';',1X,'dX',I1,'/a.u. =',1X,E12.6,';')  
 1001 FORMAT(10X,'****************************************************')
 1002 FORMAT(19X,'FLUX PROBABILITY version 0.1 beta',//,14X,
     &     'Principal author: ',/,17X,6X,
     &     'Freddy Fernandes Guimaraes.',/,14X,
     &     'Contributors:',/,17X, 
     &     'Amary Cesar, Viviane Costa Felicissimo,',/,17X,
     &     'Faris Gelmukhanov, and Yasen Velkov.')
 1003 FORMAT(/,3X,A68)
      
      END



