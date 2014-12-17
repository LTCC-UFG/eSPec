      PROGRAM NORMF
      IMPLICIT NONE
c
      INTEGER       KX
      PARAMETER     (KX = 3D0)
      INTEGER        NPRV, NPR, ZERO, MNPR, ONE
      PARAMETER     (NPRV = 4*8192, NPR = 1025, ZERO = +0.0D+0)
      PARAMETER     (MNPR = NPRV*NPR, ONE = 1.0D+0)
c 
      CHARACTER*4   CHNUM
      CHARACTER*16  INFILE
      CHARACTER*72  CHLX
      INTEGER       I, J, M, NI, NJ, IFLTH, IOTEST
      INTEGER       NIAUX, MAUX, IAUX, KXORD
      REAL*8        XL, SUM, T, DT, DT1, DT2, DW, EMF, FATE, DTAUX
      REAL*8        SGAUX, TF, WP, AUX, TMAX, TAUX, TP, LBDA, TT
      REAL*8        UAUX, VAUX, CS, SN, DEM, DE, PI, TAL
c
      REAL*8        U(NPRV), V(NPRV), EUC(NPRV), EUCK(NPRV)
      REAL*8        UM(NPRV,NPR), VM(NPRV,NPR)
      REAL*8        TV(NPRV), XKNOT(NPRV+KX), BCOEFU(NPRV), BCOEFV(NPRV)
c 
      REAL*8        SG, ECNORM, THETA, DBSVAL
c      DATA UM / MNPR * ZERO /
c      DATA VM / MNPR * ZERO /
      DATA U / NPRV* ZERO /
      DATA V / NPRV* ZERO /
      DATA EUC / NPRV* ZERO /
      DATA EUCK / NPRV* ZERO /
c  
      PI = 3.1415  
      INFILE = 'wpA_0001.dat'
      IFLTH = 12
      FATE = 4.13566727D0 !1.519829846003D-1*27.2114
      WP = 1D-1
      TP = 1550
      WRITE(*,*)'type the time for the maximum of the x-ray pulse (fs)'
      READ(*,*)TP
      WRITE(*,*)'type the the HWHM of the x-ray pulse (fs)'
      READ(*,*)TAL
      LBDA = 0.055
      KXORD = NINT(3.0D+0)
      DE = 3
      TT = 200
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
      INFILE = 'wpA_0002.dat'
      I = 2
      OPEN(2, STATUS='OLD', FILE=INFILE(1:IFLTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.ZERO)GOTO 20
      READ(2,*,END=999,ERR=999)CHLX, CHLX, DT1
      DO J=1,NJ,1
         READ(2,*,END=999,ERR=999)XL, UM(I,J), VM(I,J)
      ENDDO
c
      INFILE = 'wpA_0003.dat'
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
            INFILE = 'wpA_000'//CHNUM(4:4)//'.dat'
         ELSEIF(I.LT.100)THEN
            INFILE  = 'wpA_00'//CHNUM(3:4)//'.dat'
         ELSEIF(I.LT.1000)THEN
            INFILE = 'wpA_0'//CHNUM(2:4)//'.dat'
         ELSEIF(I.LT.10000)THEN
            INFILE = 'wpA_'//CHNUM(1:4)//'.dat'
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
      TF = (NI - 1)*DT
      TMAX = TF
      MAUX = M + 1
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
c            SGAUX = DEXP(-T*LBDA*1.51928)
c            SGAUX = DEXP(-T**2/TAUX)
            SGAUX = 1.0D+0
c            SGAUX = (2.0D0/PI*ATAN((TP - T)/TT) + 1.0D+0)/2.0D+0
            IF(J.EQ.1)write(11,*)t,theta(t,tp,tf,tal), sgaux
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
            U(I) = (UAUX*CS - VAUX*SN)*SGAUX*THETA(T,TP,TF,TAL)
            V(I) = (VAUX*CS + UAUX*SN)*SGAUX*THETA(T,TP,TF,TAL)
c            U(IAUX) = 0d0
c            V(IAUX) = 0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccU(NIAUX+I) = (UAUX*CS - VAUX*SN)*SGAUX*THETA(T,TP,TF,TAL)
cccccU(NIAUX-I+2)= (UAUX*CS + VAUX*SN)*SGAUX*THETA(T,TP,TF,TAL)
cccccV(NIAUX+I) = (UAUX*CS + VAUX*SN)*SGAUX*THETA(T,TP,TF,TAL)
cccccV(NIAUX-I+2) = (UAUX*CS - VAUX*SN)*SGAUX*THETA(T,TP,TF,TAL)
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
            UM(I,J) = U(I)*DTAUX
            VM(I,J) = V(I)*DTAUX
c
            EUCK(I) = EUCK(I) + DTAUX*DTAUX*(U(I)*U(I) + V(I)*V(I))
c
            write(22,*)i,(u(i)*u(i)+v(i)*v(i)),euck(i)
         ENDDO
         write(22,*)
c         close(2)
c         read(*,*)
      ENDDO
c
      DO I=1,NIAUX,1
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
      WRITE(12,*)'# TP=',TP
      WRITE(12,*)'# TAL=',TAL
      DO I=1,NIAUX,1
         WRITE(12,*)(I-1)*DW*FATE-DE, EUC(I), EUCK(I)
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

      FUNCTION SG(T, P)
      IMPLICIT NONE
      REAL*8 T, P, SG
c
      SG = DEXP(-T*T/P)
c
      RETURN
      END     
      
      FUNCTION THETA(T, TP, TF, TAL)
      IMPLICIT NONE
      REAL*8 T, TP, TF, THETA, TAL, PI, KL
c
      PI = 3.1416
c      TAL = 10
      KL = 1
c      THETA = (2.0D0/PI*ATAN((T - TP)/TAL) + 1.0D+0)/2.0D+0
c      THETA = 1.0D+0
      THETA = EXP(-((T - TP)/TAL)**(2*KL))
c      IF(T.GE.TP .AND. T.LE.TF)THEN
c         THETA = 1.0D+0
c      ELSE
c         THETA = 0.0D+0
c      ENDIF
c
      RETURN
      END
