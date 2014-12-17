      PROGRAM NORMF
      IMPLICIT NONE
c
      INTEGER       KX
      PARAMETER     (KX = 3D0)
      REAL*8        NPRV, NPR, ZERO, MNPR, ONE
      PARAMETER     (NPRV = 8192, NPR = 2048, ZERO = +0.0D+0)
      PARAMETER     (MNPR = NPRV*NPR, ONE = 1.0D+0)
c 
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILE
      CHARACTER*72  CHLX
      INTEGER       I, J, M, NI, NJ, IFLTH, IOTEST
      INTEGER       NIAUX, MAUX, IAUX, KXORD
      REAL*8        XL, SUM, T, DT, DT1, DT2, DW, EMF, FATE
      REAL*8        SGAUX, TF, WP, AUX, TMAX, TAUX, TP, LBDA, TT
c
      REAL*8        U(NPRV), V(NPRV), EUC(NPRV), EUCK(NPRV)
      REAL*8        UM(NPRV,NPR), VM(NPRV,NPR)
      REAL*8        TV(NPRV), XKNOT(NPRV+KX), BCOEFU(NPRV), BCOEFV(NPRV)
c 
      REAL*8        SG, ECNORM, THETA
c      DATA UM / MNPR * ZERO /
c      DATA VM / MNPR * ZERO /
      DATA U / NPRV* ZERO /
      DATA V / NPRV* ZERO /
      DATA EUC / NPRV* ZERO /
      DATA EUCK / NPRV* ZERO /
c    
      INFILE = 'wp_0001.dat'
      IFLTH = 11
      FATE = 4.13566727D0 !1.519829846003D-1*27.2114
      WP = 1D-4
      TP = 700
      LBDA = 0.055
      KXORD = NINT(3.0D+0)
c 
      I = 1
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:IFLTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(2,'(A)',END=999,ERR=999)CHLX
      DO J=1,NPR,1
         READ(2,*, END=10,ERR=10)XL, UM(I,J), VM(I,J)
      ENDDO
 10   CONTINUE
c
      NJ = J - 1
c
      INFILE = 'wp_0002.dat'
      I = 2
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX,CHLX,DT1
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      INFILE = 'wp_0003.dat'
      I = 3
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX,CHLX,DT2
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      DT = ABS(DT2 - DT1)
      write(*,*)'#dt =',dt
      
c
      DO I=4,NPRV,1
         write(*,*)'i',i,CHNUM(1:4)
         WRITE(CHNUM,'(I4)')I
         IF(I.LT.10)THEN
            INFILE = 'wp_000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE  = 'wp_00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = 'wp_0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = 'wp_'//CHNUM(1:4)//'.dat'
         ELSE
            WRITE(*,*)'Too many files to be read. Stopping!'
            STOP
         ENDIF  
         OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
         IF(IOTEST.NE.ZERO)GOTO 20
         READ(2,'(A)',END=999,ERR=999)CHLX
         DO J=1,NJ,1
            READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
         ENDDO
      ENDDO
 20   CONTINUE
      write(*,*)'#ni =',i - 1
      M = INT(LOG(FLOAT(I))/LOG(2.0D0))
      NI = 2.0D+0**M
      write(*,*)'#dw =',dw
      TF = (NI - 1)*DT
      TMAX = NI*DT
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
c
      TAUX = TMAX**2/(LOG(ONE/WP)*2.302585093)
      TT = 700
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
         DO I=1,NI,1
            TV(I) = (I - ONE)*DT
            U(I) = UM(I,J)
            V(I) = VM(I,J)
         ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
         CALL DBSINT(NI, TV, U, KX, XKNOT, BCOEFU)
         CALL DBSINT(NI, TV, V, KX, XKNOT, BCOEFV)
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         MAUX = M + 1
         NIAUX = 2.0D+0**MAUX
         DTAUX = 2*TF/NIAUX
         DW = 1.0D0/(NIAUX*DTAUX)
         DO I=1,NIAUX,1
            IAUX=2*NI-I+1
            T = (I - ONE)*DTAUX
c            SGAUX = DEXP(-T*LBDA*1.51928)
            SGAUX = DEXP(-(T-TT)**2/TAUX)
            if(T.LT.TT)SGAUX = 1
            write(11,*)t,theta(t,tp,tf),sgaux
            U(I) = DBSVAL(T, KX, XKNOT, NI, BCOEFU)
     &           *SGAUX*THETA(T,TP,TF)
            V(I) = DBSVAL(T, KX, XKNOT, NI, BCOEFV)
     &           *SGAUX*THETA(T,TP,TF)
            U(IAUX) = U(I)
            V(IAUX) = V(I)
         ENDDO
         close(11)
c         write(*,*)'Final'
c         read(*,*)
         DO I=1,NIAUX,1
            write(3,*)i,u(i),v(i),(u(i)*u(i)+v(i)*v(i))
         ENDDO
         close(3)
c
         CALL FFTNT(U, V, NIAUX, MAUX, 1)
c
         DO I=1,NI,1
            UM(I,J) = U(I)*DT
            VM(I,J) = V(I)*DT
c
            EUCK(I) = EUCK(I) + DT*DT*(U(I)*U(I) + V(I)*V(I))
            write(2,*)i,u(i),v(i),(u(i)*u(i)+v(i)*v(i)),euck(i)
         ENDDO
         close(2)
c         read(*,*)
      ENDDO
c
      DO I=1,NI,1
         EUC(I) = ZERO 
         DO J=1,NJ,1
            EUC(I) = EUC(I) + UM(I,J)*UM(I,J) + VM(I,J)*VM(I,J)
         ENDDO
      ENDDO
c
      WRITE(*,1001)
      write(*,*)dw,fate
c      read(*,*)
      AUX = ECNORM(NI, EUC)
c
      DO I=1,NI,1
         WRITE(12,*)(I-1)*DW*FATE, EUC(I)/AUX, EUCK(I)
      ENDDO
c
      STOP
c
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
c     ..
 1001 FORMAT(3X,'#w/a.u.#',5X,'#Int/arb. units#')
 1002 FORMAT(1X,E12.6,5X,E12.6)
      END

      FUNCTION SG(T,P)
      IMPLICIT NONE
      REAL*8 T, P, SG
c
      SG = DEXP(-T*T/P)
c
      RETURN
      END     
      
      FUNCTION THETA(T,TP,TF)
      IMPLICIT NONE
      REAL*8 T, TP, TF, THETA, TAL, PI, KL
c
      PI = 3.1415
      TAL = 10
      KL = 2
c      THETA = 2/PI*ATAN((T-TP)/TAL)
      THETA = EXP(-((T - TP)/TAL)**(2*KL))
c      IF(T.GE.TP .AND. T.LE.TF)THEN
c         THETA = 1.0D+0
c      ELSE
c         THETA = 0.0D+0
c      ENDIF
c
      RETURN
      END
