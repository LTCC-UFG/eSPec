      PROGRAM PUMPPROBE
      IMPLICIT NONE
*     ..
*     Purpose
*     =======
*     The RAMAN can compute the wave packet to be propagated in the final potential 
*     in simulations of x-ray Raman experiments. 
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
*     (01/10/2004) First version of NORMF written by Freddy
*     In this first version the program is able to canstruct IR-X-ray Pump-probe spectra.
*     (02/01/2006) Small modifications in the NORMF to include an external input file and 
*     the name of the program was changed for PUMPPROBE.
*     
c     **
c     ** Scalar arguments
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILEAUX, INFILE, TPPULS
      CHARACTER*72  CHLX
      INTEGER       I, J, M, NI, NJ, IFLTH, IFLTHT, IOTEST
      INTEGER       NIAUX, MAUX, IAUX, KXORD
      REAL*8        XL, SUM, T, DT, DT1, DT2, DW, EMF, DTAUX
      REAL*8        TF, AUX, TMAX, UAUX, VAUX, CS, SN, DEM, DE
c     ** 
c     ** Parameters 
      INTEGER       KX, NPRV, NPR, MNPR
      REAL*8        ZERO, ONE, PI, FATE 
      PARAMETER     (NPRV = +4D+0*8192D+0, 
                     NPR  = +1025D+0, 
                     KX   = +3D+0, 
                     MNPR = +NPRV*NPR, 
                     ZERO = +0.0D+0, 
                     ONE  = +1.0D+0, 
                     FATE = +4.13566727D+0, 
                     PI   = +3.14159265358979323846D+0)
c     **
c     ** Array arguments
      REAL*8        U(NPRV), V(NPRV), EUC(NPRV), EUCK(NPRV)
      REAL*8        TV(NPRV), XKNOT(NPRV+KX), BCOEFU(NPRV), BCOEFV(NPRV)
      REAL*8        UM(NPRV,NPR), VM(NPRV,NPR)
c     **
c     ** External functions 
      INTEGER       ICHLENGTH
      REAL*8        THETA, ECNORM, DBSVAL
c     **
c     ** External subroutines 
      EXTERNAL      FFTNT, DBSINT, DBSNAK
c     **
c     ** Intrinsic functions 
      INTRINSIC     DSIN, DCOS
c     ..
c     .. Starting program
      WRITE(*,1001)
      WRITE(*,1002)
      WRITE(*,1001)
      WRITE(*,1003)'Initializing...' 
c     .. 
c     .. Starting default logical values
c     .. 
c     .. Starting character default values
      INFILEAUX = 'wp_000'
      FLTH =  ICHLENGTH(INFILEAUX, 1) 
      IFLTHT = IFLTH + 4
c     .. 
c     .. Starting scalar default values
      KXORD = NINT(3.0D+0)
c     ..
c     .. Starting integer arrays values
c     .. 
c     .. Starting real arrays values  
      DATA U        / NPRV* ZERO /
      DATA V        / NPRV* ZERO /
      DATA EUC      / NPRV* ZERO /
      DATA EUCK     / NPRV* ZERO /
c     .. Getting input parameters
      WRITE(*,1003)'Getting input parameters...'
 1    CONTINUE
      READ(*,*)RDAUX 
      IF(RDAUX(1:10).EQ.'*INPUTFILE')THEN
         READ(*,*)INFILEAUX
         IFLTH =  ICHLENGTH(INFILEAUX, 1) 
         IFLTHT = IFLTH + 4
         GOTO 1
      ELSEIF(RDAUX(1:11).EQ.'*OUTPUTFILE')THEN
         READ(*,*)OUTFILE
         GOTO 1
      ELSEIF(RDAUX(1:6).EQ.'*GAMMA')THEN
         READ(*,*)GAMMA
         GOTO 1
      ELSEIF(RDAUX(1:6).EQ.'*EXC_ENERGS')THEN 
         READ(*,*,END=135,ERR=135)(WC(I), I=1,MXNXE,1)
 135     CONTINUE
         NXERG = I - ONE
         GOTO 1
      ELSEIF(RDAUX(1:4).EQ.'*END')THEN
         CONTINUE
      ELSE
         WRITE(*,*)'<<<>>> Key word,', RDAUX, 
     &        ' was not recognized. <<<>>>'
         WRITE(*,*)'Please, check the imput file and resubmit!'
         STOP
      ENDIF
      WRITE(*,*)'Input parameters got.'
c 
      WRITE(*,*)
      WRITE(*,*)'Input description:'
      WRITE(*,*)'""""""""""""""""""'
      WRITE(*,*)'   The wave packet input file name is: ',
     &     INFILEAUX(1:IFLTH)//'?.dat'
      WRITE(*,*)'   The life time of the core excited state is:',
     &     GAMMA,' eV or ',GAMMA,' 1/fs'
      WRITE(*,*)'   It was chosen', NXERG,' distinct excitation energies
     &:'
      DO I=1,NXERG,3
         WRITE(*,*)'      ',(J, WC(J), J=I,I+3,1)
      ENDDO
c     ..
c     .. Reading wave packets
      WRITE(*,1003)'Acquiring wave packets...'
      I = 1
      INFILE = INFILEAUX(1:IFLTH)//'1.dat'
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTHT), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,'(A)',END=999,ERR=999)CHLX
      READ(2,*, END=10,ERR=999)XI, UM(I,J), VM(I,J)
      READ(2,*, END=10,ERR=999)X2, UM(I,J), VM(I,J)
      DX = X2 - XI
      WRITE(*,*)'   The space between points is:', DX
      DO J=3,NPR,1
         READ(2,*, END=10,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
 10   CONTINUE
c
      NJ = J - 1
      WRITE(*,*)'   Each wave packet file contains',NJ ,' points.' 
c
      I = 2
      INFILE = INFILEAUX(1:IFLTH)//'2.dat'
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT1
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      INFILE = INFILEAUX(1:IFLTH)//'3.dat'
      I = 3
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX,CHLX,DT2
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      DT = ABS(DT2 - DT1)
      WRITE(*,*)'   Delta time between the wave packets acquired =',DT            
c
      DO I=4,NPRV,1
         write(*,*)'i',i,CHNUM(1:4)
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = INFILEAUX(1:IFLTH)//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE = INFILEAUX(1:IFLTH-1)//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = INFILEAUX(1:IFLTH-2)//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = INFILEAUX(1:IFLTH-3)//CHNUM(1:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF  
         OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)GOTO 20
         READ(2,'(A)',END=999,ERR=999)CHLX
         DO J=1,NJ,1
            READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
            IF(J.GT.NJ-5)THEN
               IF((UM(I,J) .OR. VM(I,J)).GT.1D-2)THEN                 
                  WRITE(*,*)'<<<>>> The wave packet touch the artifitial
     & barrier! Stopping. <<<>>>'
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDDO
 20   CONTINUE
      WRITE(*,*)'   Total number of wave packets acquired =',i - 1
      WRITE(*,*)'Wave packets acquired.'
c     ..
c     .. Compute the auxiliary parameters to perform the FFT
      WRITE(*,1003)'Determining parameters to compute the FFT...'
      M = INT(LOG(FLOAT(I))/LOG(2.0D0))
      NI = 2.0D+0**M
      TF = (NI - 1)*DT
      TMAX = TF
      MAUX = M + 1
      IF((THETA(TPPULS,T,TP,TF,TAL,KL).OR.THETA(TPPULS,T,TP,TF,TAL,KL)) 
     &     .GT.1D-6)THEN
         WRITE(*,*)'<<<>>> The FFT is probably wrong <<<>>> '
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(*,*)'log2(2*ni)',M
      write(*,*)'#ni =',ni,'; #nj =',nj
c      read(*,*)
c      do i=1,ni,1
c         write(*,*)'i',i
c         do j=1,nj,1
c            write(31,*)j, um(i,j)*um(i,j)+vm(i,j)*vm(i,j)
c         enddo
c         write(31,*)
c         read(*,*)
c         if(i.eq.1000)close(31)
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      TAUX = TMAX**2/(LOG(ONE/WP)*2.302585093)
      NIAUX = 2.0D+0**MAUX
      DTAUX = TF/NIAUX
      DW = 1.0D0/(NIAUX*DTAUX)
c      DW = 1.0D0/(2.0D+0*NIAUX*DTAUX)
      write(*,*)'#dw =',dw
      DEM = DE/6.58211889D-1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc 

      DO I=1,NI,1
         TV(I) = (I - ONE)*DT
      ENDDO
      CALL DBSNAK(NI, TV, KX, XKNOT) 
cccccccccccccccccccccccccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      DO J=1,NJ,1
         write(*,*)'j',j
         IF(M.NE.MAUX)THEN
            DO I=1,NI,1
               TV(I) = (I - ONE)*DT
               U(I) = UM(I,J)
               V(I) = VM(I,J)
            ENDDO 
            CALL DBSINT(NI, TV, U, KX, XKNOT, BCOEFU)
            CALL DBSINT(NI, TV, V, KX, XKNOT, BCOEFV)
         ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
         DO I=1,NIAUX,1
            IAUX=2*NIAUX-I+1
            T = (I - ONE)*DTAUX
            IF(J.EQ.1)write(11,*)t,theta(TPPULS,t,tp,tf, TAL, KL)
            IF(DEM.EQ.0)THEN
               CS = 1.0D+0
               SN = 0.0D+0
            ELSE
               CS = DCOS(DEM*T)
               SN = DSIN(DEM*T)
            ENDIF
            IF(M.EQ.MAUX)THEN
               UAUX = UM(I,J)
               VAUX = VM(I,J)
            ELSE
               UAUX = DBSVAL(T, KX, XKNOT, NI, BCOEFU)
               VAUX = DBSVAL(T, KX, XKNOT, NI, BCOEFV)
            ENDIF               
            U(I) = (UAUX*CS - VAUX*SN)*THETA(TPPULS,T,TP,TF,TAL,KL)
            V(I) = (VAUX*CS + UAUX*SN)*THETA(TPPULS,T,TP,TF,TAL,KL)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccU(NIAUX+I) = (UAUX*CS - VAUX*SN)*SGAUX*THETA(T,TP,TF)
cccccU(NIAUX-I+2)= (UAUX*CS + VAUX*SN)*SGAUX*THETA(T,TP,TF)
cccccV(NIAUX+I) = (UAUX*CS + VAUX*SN)*SGAUX*THETA(T,TP,TF)
cccccV(NIAUX-I+2) = (UAUX*CS - VAUX*SN)*SGAUX*THETA(T,TP,TF)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ENDDO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cccccU(1) = U(2)
cccccV(1) = V(2)  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         write(*,*)'Final'
c         read(*,*)
c         DO I=1,2*NIAUX,1
c            write(3,*)i,u(i),v(i),(u(i)*u(i)+v(i)*v(i))
c         ENDDO
c         close(3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
         CALL FFTNT(U, V, NIAUX, MAUX, 1)
c
         DO I=1,NIAUX,1
            UM(I,J) = U(I)*DT
            VM(I,J) = V(I)*DT
c
            EUCK(I) = EUCK(I) + DT*DT*(U(I)*U(I) + V(I)*V(I))
c
            write(22,*)i,(u(i)*u(i)+v(i)*v(i)),euck(i)
         ENDDO
         write(22,*)
c         close(2)
c         read(*,*)
      ENDDO
c     
      WRITE(*,1001)
      write(*,*)dw,fate
      DO I=1,NIAUX,1
         EUC(I) = ZERO 
         DO J=1,NJ,1
            EUC(I) = EUC(I) + UM(I,J)*UM(I,J) + VM(I,J)*VM(I,J)
         ENDDO
         WRITE(12,*)(I-1)*DW*FATE-DE, EUC(I), EUCK(I)
      ENDDO
c
      STOP
c
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
c     ..
 1011 FORMAT(3X,'#w/a.u.#',5X,'#Int/arb. units#')
 1012 FORMAT(1X,E12.6,5X,E12.6)
 1001 FORMAT(10X,'****************************************************')
 1002 FORMAT(25X,'PUMPPROBE version 0.1 beta',//,14X,
     &     'Principal author: ',/,17X,6X,
     &     'Freddy Fernandes Guimaraes.',/,14X,
     &     'Contributors:',/,17X, 
     &     'Amary Cesar, Viviane Costa Felicissimo,',/,17X,
     &     'and  Faris Gelmukhanov.')
 1003 FORMAT(/,3X,A68)
      END
      
      FUNCTION THETA(TPPULS, T, TP, TF, TAL, KL)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
      CHARACTER*(*)  TPPULS
      REAL*8         T, TP, TF, TAL, KL
*     ..
*     Purpose
*     =======
*     Compute the theta function, which represents the life time decay of the 
*     core excited state.
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
*     (01/09/2004) First version THETA written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        PI
      PARAMETER     (PI = 3.14159265358979323846)
c     **
c     ** Local scalars 
      REAL*8        THETA
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
      INTRINSIC     DEXP, ATAN
c     .. Start program      
      IF(TPPULS(1:6).EQ.'.GAUSS')THEN
         THETA = DEXP(-((T - TP)/TAL)**(2*KL))
      ELSEIF(TPPULS(1:5).EQ.'.ATAN')THEN 
         THETA = (2.0D0/PI*ATAN((T - TP)/TAL) + 1.0D+0)/2.0D+0
      ELSEIF(TPPULS(1:8).EQ.'.COSTANT')THEN
         THETA = 1.0D+0
      ELSEIF(TPPULS(1:5).EQ.'.STEP')THEN
         IF(T.GE.TP .AND. T.LE.TF)THEN
            THETA = 1.0D+0
         ELSE
            THETA = 0.0D+0
         ENDIF
      ELSE
         WRITE(*,*)'The specified type of pulse was not found!'
         STOP
      ENDIF
c     ..
      RETURN
      END
