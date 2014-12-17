      PROGRAM NORMF
      IMPLICIT NONE
c
      REAL*8        NPRV, NPR, ZERO, ONE
      PARAMETER     (NPRV = 4*8192, NPR = 10025, ZERO = +0.0D+0)
      PARAMETER     (ONE = 1.0D+0)
c 
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILE, OUTFILE
      CHARACTER*72  CHLX
      INTEGER       I, J, NJ, IFLTH, OFLTH, IOTEST, IST, M
      REAL*8        XL, TAU, DT
c
      REAL*8        U(NPRV)
c 
      DATA U / NPRV* ZERO /

c  
      INFILE = 'eigvc_0001.dat'
      IFLTH = 14
      OUTFILE ='lteigvc_0001.dat'
      OFLTH = 16
      TAU = 20.56  ! fs
      DT = 0.366   ! Step from the files eigvc_????.dat 
      IST = 5      ! Step to output files

c
      I = 1
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      OPEN(UNIT=3, STATUS='UNKNOWN', FILE=OUTFILE(1:OFLTH))
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,'(A)',END=999,ERR=999)CHLX
      DO J=1,NPR,1
         READ(2,*, END=10,ERR=10)XL, U(J)
         WRITE(3,*)XL, U(J)
      ENDDO
      CLOSE(3)
 10   CONTINUE
c
      NJ = J - 1
      write(*,*)'#nj =',nj
c
      M = 1
      DO I=2,NPRV,IST
         M = M + 1
         WRITE(CHNUM,'(I4)')I
         write(*,*)'i',m,i,CHNUM(1:4)
         IF(I.LT.10)THEN
            INFILE = 'eigvc_000'//CHNUM(4:4)//'.dat'
            OUTFILE = 'lteigvc_000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE = 'eigvc_00'//CHNUM(3:4)//'.dat'
            OUTFILE = 'lteigvc_00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = 'eigvc_0'//CHNUM(2:4)//'.dat'
            OUTFILE = 'lteigvc_0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = 'eigvc_'//CHNUM(1:4)//'.dat'
            OUTFILE = 'lteigvc_0'//CHNUM(2:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF  
         OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
         OPEN(UNIT=3, STATUS='UNKNOWN', FILE=OUTFILE(1:OFLTH))
         IF(IOTEST.NE.ZERO)GOTO 20
         READ(2,'(A)',END=999,ERR=999)CHLX
         DO J=1,NJ,1
            READ(2,*,END=999,ERR=999)XL, U(J)
            WRITE(3,*)XL, U(J)*EXP(-(I - ONE)*DT/TAU)
         ENDDO
         CLOSE(3)
      ENDDO
 20   CONTINUE
      write(*,*)'#ni =',i - 1
c
      RETURN
c
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP 
c
      END
