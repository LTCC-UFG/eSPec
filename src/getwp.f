      SUBROUTINE GETWPS(INFILEAUX, NPR, NF, NJ, DT, XI, DX, UM, VM)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      CHARACTER*(*) INFILEAUX
      INTEGER       NPR, NF, NJ
      REAL*8        DT, XI, DX
c     **
c     ** Array arguments
      REAL*8        UM(NPR,*), VM(NPR,*)
c     **
*     ..
*     Purpose
*     =======
*     Read the WPs from the eSPec output.
*
*     ..
*     Arguments
*     =========
*
*
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (01/03/2006) First version GETWPS written by Freddy.   
*
c     **
c     ** Parameters 
      REAL*8        ZERO
      PARAMETER     (ZERO = +0.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILE
      CHARACTER*72  CHLX
      INTEGER       I, J, IFLTHT, IOTEST, NPRV, IFLTH
      REAL*8        DT1, DT2, XL, X2
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER       
c      REAL*8        U0(N), V0(N)  
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       ICHLENGTH      
c      REAL*8        FCORREL, EIGENERG    
c     **
c     ** External subroutines 
c      EXTERNAL      PSODAUX
c     **
c     ** Intrinsic functions 
c      INTRINSIC     SQRT 
c     .. Start program      
c     ..
c     .. Starting values
c     .. Starting integer local values
      NPRV = NPR*1000
c     ..
c     .. Computing the size of the input file name
      IFLTH = ICHLENGTH(INFILEAUX,0)
c     .. 
c     .. Reading wave packets
      WRITE(*,1003)'Acquiring wave packets...'
      I = 1
      INFILE = INFILEAUX(1:IFLTH)//'_0001.dat'
      IFLTHT = IFLTH + 9
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTHT), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,'(A)',END=999,ERR=999)CHLX
      READ(2,*, END=10,ERR=999)XI, UM(I,1), VM(I,1)
      READ(2,*, END=10,ERR=999)X2, UM(I,2), VM(I,2)
      DX = X2 - XI 
      WRITE(*,*)'   The space between points is', DX,' a.u.'
      DO J=3,NPR,1
         READ(2,*, END=10,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
 10   CONTINUE
c
      NJ = J - 1
      WRITE(*,*)'   Each wave packet file contains',NJ ,' points.' 
c
      I = 2
      INFILE = INFILEAUX(1:IFLTH)//'_0002.dat'
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT1
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      INFILE = INFILEAUX(1:IFLTH)//'_0003.dat'
      I = 3
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTHT), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX,CHLX,DT2
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      DT = ABS(DT2 - DT1)
      WRITE(*,*)'   Delta time between the acquired wave packets =',DT,
     &     ' fs'
c
      DO I=4,NPRV,1
c         write(*,*)'i',i,CHNUM(1:4)
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = INFILEAUX(1:IFLTH)//'_000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE = INFILEAUX(1:IFLTH)//'_00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = INFILEAUX(1:IFLTH)//'_0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = INFILEAUX(1:IFLTH)//'_'//CHNUM(1:4)//'.dat'
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
               IF(UM(I,J).GT.1D-2 .OR. VM(I,J).GT.1D-2)THEN                 
                  WRITE(*,*)'<<<>>> The wave packet touch the artifitial
     & barrier! <<<>>>'
                  WRITE(*,*)'Stopping in file:', I
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDDO
 20   CONTINUE
      NF = I - 1
      WRITE(*,*)'   Total number of wave packets acquired =',NF
      WRITE(*,*)'Wave packets acquired.'
c
 1003 FORMAT(/,3X,A68)   
c
      RETURN
c
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP 
c
      END
