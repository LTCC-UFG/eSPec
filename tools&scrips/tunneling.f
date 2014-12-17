      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILE
      CHARACTER*72  LX
C
      INFILE = 'ReIm_0001.dat'
      IFLTH = 13
      IOTEST = 0
      OPEN(1, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.0)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTH),
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(1,*)LX,LX,T,LX
      TLD = 0.0D+0
      DO J=1,10000,1
         READ(1,*,END=10)CR, SR, U, V
         IF(SR.GT.CR/2)TLD = TLD+ U*U + V*V
      ENDDO
 10   CONTINUE
      CLOSE(1) 
      WRITE(*,*)T, TLD
      NT = J - 1
      WRITE(*,*)'NT =',nt
C
      DO I=1,1000,1
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
         OPEN(1, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
         IF(IOTEST.NE.0)GOTO 20
         TLD = 0.0D+0
         READ(1,*)LX,LX,T,LX
         DO J=1,NT,1
            READ(1,*)CR, SR, U, V
            IF(SR.GT.CR/2)TLD = TLD+ U*U + V*V
         ENDDO
         CLOSE(1)
         WRITE(*,*)T, TLD
      ENDDO
 20   CONTINUE
      WRITE(*,*)'Finish I=',I
      END
C
