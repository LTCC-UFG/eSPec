      PROGRAM PROBABILITY
      IMPLICIT NONE 
c
      INTEGER       NDIM, MXDCT
      PARAMETER     (NDIM = 3.0D+0, MXDCT = 1.0D+5)
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = 0.0D+0, ONE = +1.0D+0)
c
      CHARACTER*4   CHNUM  
      CHARACTER*72  INFILEI, INFILE, CHLX
      INTEGER       NC, NINI, NFIN, ND, IOTEST, I, J, K, L, NT
      REAL*8        DT, DT1, DT2, RXA, RYA, PROB1, PROB2, PROB3
      REAL*8        RX, RY
c
      LOGICAL       CHSTEP(NDIM)
      INTEGER       NP(NDIM)
      REAL*8        XI(NDIM), XF(NDIM), SH(NDIM), U(MXDCT), V(MXDCT)
c
      INTEGER       ICHLENGTH      
c
      CHNUM   = '   ' 
      INFILEI = '   '
      INFILE  = '   '
      CHLX    = '   '
c
      NINI = ZERO
      NFIN = ZERO
      ND   = ZERO
      NT   = ZERO
c
      DT  = ZERO 
      DT1 = ZERO
      DT2 = ZERO
      RXA = ZERO
      RYA = ZERO
      RX  = ZERO
      RY  = ZERO
c     ..
      DATA CHSTEP  / NDIM* .TRUE. /
c
      DATA NP / NDIM* ZERO /
c
      DATA XI / NDIM* ZERO /
      DATA XF / NDIM* ZERO /
      DATA SH / NDIM* ZERO /
      DATA U  / MXDCT* ZERO /
      DATA V  / MXDCT* ZERO /
c 
      READ(*,*)INFILEI
      WRITE(*,*)'The input file names start with', INFILEI
      READ(*,*)NINI, NFIN
      WRITE(*,*)'The enumeration starts in', NINI,
     &     ' and finish in', NFIN
      READ(*,*)ND
      WRITE(*,*)'The files contain',ND,' dimensional data.'
      READ(*,*)RXA, RYA
      WRITE(*,*)'Probability division in dimension 1',RXA
      WRITE(*,*)'   and in dimension 2',RYA,'.'
c
c     .. Check the size of the files
      NC = ICHLENGTH(INFILEI, 0)
      WRITE(CHNUM,'(I4)')NINI
      INFILE = INFILEI(1:NC)//'_000'//CHNUM(4:4)//'.dat   '
      NC = ICHLENGTH(INFILE, 0)
      OPEN(2, STATUS='OLD', FILE=INFILE(1:NC), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:NC), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT1
      READ(2,*, END=10,ERR=10)(XI(J), J=ND,1,-1), U(1), V(1)
      WRITE(*,*)'Initial grid points:',('D',J,'= ',XI(J),';   ',
     &     J=ND,1,-1)
c
      DO I=2,MXDCT,1
         READ(2,*, END=10,ERR=10)(XF(J), J=ND,1,-1), U(I), V(I)
c     .. Get the step in space
         DO J=1,ND,1
            IF(CHSTEP(J) .AND. XF(J).NE.XI(J))THEN
               CHSTEP(J) = .FALSE.
               SH(J) = XF(J) - XI(J)
            ENDIF
         ENDDO
      ENDDO
 10   CONTINUE
      WRITE(*,*)'Final grid points:',('D',J,'= ',XF(J),';   ',
     &     J=ND,1,-1)
      WRITE(*,*)'Step between points in the grid:'
     &     ,('D',J,'= ',SH(J),';   ',J=ND,1,-1)
      NT = I - 1
      DO I=1,ND,1
         NP(I) = (XF(I) - XI(I))/SH(I) + ONE
      ENDDO
      WRITE(*,*)'There are ',NT,' points in each file.'
      WRITE(*,*)'   with ',(NP(J),' in the dimension', J,';   ',
     &     J=ND,1,-1)
c
c     .. Check the time step of the files
      NC = ICHLENGTH(INFILEI, 0)
      WRITE(CHNUM,'(I4)')NINI + 1
      INFILE = INFILEI(1:NC)//'_000'//CHNUM(4:4)//'.dat   '
      NC = ICHLENGTH(INFILE, 0)
      OPEN(2, STATUS='OLD', FILE=INFILE(1:NC), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:NC), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT2
      DT = ABS(DT2 - DT1)
      WRITE(*,*)'The time step between files is:',DT,'.'
c
      DO I=NINI,NFIN,1
         NC = ICHLENGTH(INFILEI, 0)
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = INFILEI(1:NC)//'_000'//CHNUM(4:4)//'.dat   '
         ELSEIF(I.LT.100)THEN
            INFILE = INFILEI(1:NC)//'_00'//CHNUM(3:4)//'.dat   '
         ELSEIF(I.LT.1000)THEN
            INFILE = INFILEI(1:NC)//'_0'//CHNUM(2:4)//'.dat   '
         ELSEIF(I.LT.10000)THEN
            INFILE = INFILEI(1:NC)//'_'//CHNUM(1:4)//'.dat   '
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF
c     ..
c     .. Read data from file
         CALL RDPT2(INFILE, NT, ND, XI, XF, U, V)
c     ..
c     .. Perform the probability calculation
         PROB1 = ZERO
         PROB2 = ZERO
         PROB3 = ZERO 
         DO J=1,NP(1),1
            RY = XI(2) + (J - 1)*SH(2)
            DO K=1,NP(2),1
               RX = XI(1) + (K - 1)*SH(1)
               L = K + (J - 1)*NP(1)              
c               IF(RX.GT.RXA)THEN
c                  PROB3 = PROB3 + U(L)*U(L) + V(L)*V(L)
c               ELSE
               IF(RY.GT.RYA)THEN
                  PROB2 = PROB2 + U(L)*U(L) + V(L)*V(L)
               ELSE
                  PROB1 = PROB1 + U(L)*U(L) + V(L)*V(L)
               ENDIF
            ENDDO
         ENDDO
         WRITE(*,*)DT1 + (I - 1)*DT, PROB1, PROB2, PROB3
      ENDDO
c
      RETURN
c
 999  WRITE(*,*)'Input file error!'
c
      END


!  f77 /home/freddy/espec_v0.5/src/probil.f /home/freddy/espec_v0.5/src/chlength.f /home/freddy/espec_v0.5/src/rdpt2.f
!  a.out < inp > out2.log
