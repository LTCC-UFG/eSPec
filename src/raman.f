      PROGRAM RAMAN
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
*     (02/01/2006) Modifications in the NORMF for the the calculation of the initial Wave 
*     packet to computed the X-ray RAMAN spectrum by Freddy. The name of this new program 
*     is RAMAN.
*     
c     **
c     ** Scalar arguments
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILEAUX, INFILE
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
      INTRINSIC     SIN, COS
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
 10   READ(*,*)RDAUX 
      IF(RDAUX(1:10).EQ.'*INPUTFILE')THEN
         READ(*,*)INFILEAUX
         IFLTH =  ICHLENGTH(INFILEAUX, 1) 
         IFLTHT = IFLTH + 4
         GOTO 1
      ELSEIF(RDAUX(1:11).EQ.'*OUTPUTFILE')THEN
         READ(*,*)OUTFILE
      ELSEIF(RDAUX(1:6).EQ.'*GAMMA')THEN
         READ(*,*)GAMMA
      ELSEIF(RDAUX(1:6).EQ.'*EXC_ENERGS')THEN 
         READ(*,*,END=135,ERR=135)(WC(I), I=1,MXNXE,1)
 135     CONTINUE
         NXERG = I - ONE
      ELSEIF(RDAUX(1:4).EQ.'*END')THEN
         GOTO 20
      ELSE
         WRITE(*,*)'<<<>>> Key word,', RDAUX, 
     &        ' was not recognized. <<<>>>'
         WRITE(*,*)'Please, check the imput file and resubmit!'
         STOP
      ENDIF
      GOTO 10
 20   CONTINUE
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
     & barrier! <<<>>>'
                  TF = (I - 2)*DT
                  IF(THETA(T,GAMMA).GT.1D-6)THEN
                     WRITE(*,*)'<<<>>> The decay is not enough to avoid 
     & this. Stopping!. <<<>>>'
                     
                     STOP
                  ENDIF
                  WRITE(*,*)'The decay rate allow the program to continue
     &!'
                  I = I - 1
                  GOTO 20
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
      IF(THETA(TF,GAMMA).GT.1D-6)THEN
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
 30   CONTINUE
      NIAUX = 2.0D+0**MAUX
c      IF(THETA(TF,GAMMA).GT.1D-4)THEN
c         WRITE(*,*)'<<<>>> The FFT is probabily wrong <<<>>> '
c         WRITE(*,*)'Increasing M value...'
c         MAUX = MAUX + 1
c         GOTO 30
c      ENDIF
      DTAUX = TF/NIAUX
      DW = 1.0D0/(NIAUX*DTAUX)
      write(*,*)'#dw =',dw
      DEM = DE/6.58211889D-1
c
      DO I=1,NI,1
         TV(I) = (I - ONE)*DT
      ENDDO
c
      DO I=1,NIAUX,1
         WV(I) = (I - ONE)*DW
      ENDDO   
      WRITE(*,*)'Parameters to compute the FFT determined.'
c     ..
c     .. Output title
      WRITE(*,*)(WC(I),'        ', I=1,NXERG,1) 
c     ..
c     .. Compute the FFT for the x_i coordinate 
      DO J=1,NJ,1
c     ..
c     .. Perform a spline procedure in case the number of points do not be enough
         IF(M.NE.MAUX)THEN
            CALL DBSNAK(NI, TV, KX, XKNOT)
            DO I=1,NI,1
               TV(I) = (I - ONE)*DT
               U(I) = UM(I,J)
               V(I) = VM(I,J)
            ENDDO 
            CALL DBSINT(NI, TV, U, KX, XKNOT, BCOEFU)
            CALL DBSINT(NI, TV, V, KX, XKNOT, BCOEFV)
         ENDIF
c     ..
c     .. Build the functions to be Fourier transformed
         DO I=1,NIAUX,1
            IAUX1 = NIAUX - I + ONE
            IAUX2 = NIAUX + I - ONE
            T = (I - ONE)*DTAUX
            IF(J.EQ.1)write(11,*)t,theta(t,gamma)
            IF(M.EQ.MAUX)THEN
               UAUX = UM(I,J)
               VAUX = VM(I,J)
            ELSE
               UAUX = DBSVAL(T, KX, XKNOT, NI, BCOEFU)
               VAUX = DBSVAL(T, KX, XKNOT, NI, BCOEFV)
            ENDIF               
            U(IAUX2) =  UAUX*CS*THETA(T,GAMMA)
            U(IAUX1) =  UAUX*CS*THETA(T,GAMMA)
            V(IAUX2) =  VAUX*CS*THETA(T,GAMMA)
            U(IAUX1) = -UAUX*CS*THETA(T,GAMMA)
         ENDDO
c     ..
c     .. Perform the Fast Fourier Transform in the x coordinate, which correspond to J 
         CALL FFTNT(U, V, 2*NIAUX, MAUX+1, 1)



c     ..
c     .. Spline the the Fourier Transform in the x coordinate
         CALL DBSNAK(NIAUX, WV, KX, XKNOT)
         CALL DBSINT(NIAUX, WV, U, KX, XKNOT, BCOEFU)
         CALL DBSINT(NIAUX, WV, V, KX, XKNOT, BCOEFV)
c     ..
c     .. Compute the initial wave packet in the desired excitation energies
         DO I=1,NXERG,1
            UA(I) = DBSVAL(WC(I), KX, XKNOT, NIAUX, BCOEFU)
            VA(I) = DBSVAL(WC(I), KX, XKNOT, NIAUX, BCOEFV)
         ENDDO
         WRITE(*,*)(J - ONE)*DX + XI, 
     &        (UA(I), VA(I),'        ' ,I=1,NXERG,1)
      ENDDO
      STOP
c
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
c     ..
 1011 FORMAT(3X,'#w/a.u.#',5X,'#Int/arb. units#')
 1012 FORMAT(1X,E12.6,5X,E12.6)
 1001 FORMAT(10X,'****************************************************')
 1002 FORMAT(25X,'RAMAN version 0.1 beta',//,14X,
     &     'Principal author: ',/,17X,6X,
     &     'Freddy Fernandes Guimaraes.',/,14X,
     &     'Contributors:',/,17X, 
     &     'Amary Cesar, Viviane Costa Felicissimo,',/,17X,
     &     'and  Faris Gelmukhanov.')
 1003 FORMAT(/,3X,A68)
      END

      
      FUNCTION THETA(T, GAMMA)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
      REAL*8         T, GAMMA
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
c      REAL*8        PI
c      PARAMETER     (PI = 3.14159265358979323846)
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
      INTRINSIC     EXP
c     .. Start program      
      THETA = EXP(-T*GAMMA)
c     ..
      RETURN
      END
