c Este programa calula o valor esperado
      PROGRAM project
      IMPLICIT NONE
c
      REAL*8        NPRV, NPR, ZERO, MNPR, ONE
      PARAMETER     (NPRV = 2*8192, NPR = 1025, MNPR = NPRV*NPR)
      PARAMETER     (ZERO = +0.0D+0, ONE = 1.0D+0)
c 
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILE
      CHARACTER*72  CHLX
      INTEGER       I, J, NJ, NI, NII, NIF, NIS, IFLTH, IOTEST

      REAL*8        DT, TI, TF, TS, XMAX, UMAX, XAV
c
      REAL*8        XL(NPRV), UM(NPRV,NPR), VM(NPRV,NPR)
c 
c      DATA UM / MNPR * ZERO /
      DATA XL / NPRV* ZERO /
c  
c      INFILE = 'eigvc_0001.dat'
      IFLTH = 13
c 
      TI = 700
      TF = 900
      TS = 5
      XAV = ZERO
      DT = ZERO
c
      
c
      I = 1
      INFILE = 'ReIm_0001.dat'
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,'(A)',END=999,ERR=999)CHLX
      XMAX = ZERO
      UMAX = ZERO
      DO J=1,NPR,1
         READ(2,*, END=10,ERR=10)XL(J), UM(I,J), VM(I,J)
      ENDDO
      WRITE(1,*)XMAX,DT
 10   CONTINUE
      TS = DT
c
      NJ = J - 1
      WRITE(*,*)'Number of lines in the first file: ',NJ
c
      DO J=1,NJ,1
         XAV = XAV + UM(I,J)*XL(J)*UM(I,J) + VM(I,J)*XL(J)*VM(I,J)
      ENDDO
      WRITE(1,*)DT, XAV
c
      I = 2
      INFILE = 'ReIm_0002.dat'
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT
      XMAX = ZERO
      UMAX = ZERO
      DO J=1,NJ,1
         READ(2,*, END=999,ERR=999)XL(J), UM(I,J), VM(I,J)
      ENDDO
c
      DT = DT - 0.5D-4
      WRITE(*,*)'Time step: ',DT
      XAV = ZERO
      DO J=1,NJ,1
         XAV = XAV + UM(I,J)*XL(J)*UM(I,J) + VM(I,J)*XL(J)*VM(I,J)
      ENDDO
      WRITE(1,*)DT, XAV
c
      DO I=3,NPRV,1
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = 'ReIm_000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE  = 'ReIm_00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = 'ReIm_0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = 'ReIm_'//CHNUM(1:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF  
         OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)GOTO 20
         READ(2,'(A)',END=999,ERR=999)CHLX
         XMAX = ZERO
         UMAX = ZERO
         DO J=1,NJ,1
            READ(2,*,END=999,ERR=999)XL(J), UM(I,J), VM(I,J)
         ENDDO
c
         XAV = ZERO
         DO J=1,NJ,1
            XAV = XAV + UM(I,J)*XL(J)*UM(I,J) + VM(I,J)*XL(J)*VM(I,J)
         ENDDO
         WRITE(1,*)(I - ONE)*DT, XAV
      ENDDO
 20   CONTINUE  
      NI = I - 1
      WRITE(*,*)'Number of files: ',NI
c
      NII = INT(TI/DT) + ONE
      NIF = INT(TF/DT) + ONE
      NIS = INT(TS/DT) 
      IF(NIS.LT.1)NIS = ONE
      WRITE(*,*)'File steps: ',NIS
c      write(*,*)NII, NIF, NIS
      DO I=NII,NIF,NIS
         DO J=1,NJ,1
c            WRITE(3,1001)(I - ONE)*DT, XL(J), UM(I,J)+0.5
            WRITE(3,1001)XL(J), UM(I,J)+0.5
         ENDDO
         WRITE(3,*)
      ENDDO
c
      STOP
c
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
 1001 FORMAT(1X,G12.5,3X,E12.6,3X,E12.6)
c 1002 FORMAT(A14,G12.6)
c     ..
      END
