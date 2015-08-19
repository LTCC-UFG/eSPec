      PROGRAM TESTMTRX
      implicit none
      INTEGER I,J,K,NP(3),N,L,M
      REAL*8 APAUX,AP(100000000),SHT,SHM(3)
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1) 
CCC TESTFUNCTION
      REAL*8 US(10000),AU(10000),VPOT(10000),A
      REAL*8 AU2(10000),ANS(10000),APP(10000,10000)
      REAL*8 X,XI,XF,Y,YI,YF,STEPX,STEPY,DIFF,B
      REAL*8 NFXX,NFYY,NFXY,FXX,FYY,FXY,GAUSS
      REAL*8 NUMDERIV,NANS
CC-------------------
      XI=-0.5D+1
      XF=0.5D+1
      YI=-0.3D+1
      YF=0.3D+1
cc----------------

      NP(1)=10
      NP(2)=10
      N=NP(1)*NP(2)
      SHM(1)=-ONE
      SHM(2)=-ONE
      SHM(3)=-ONE
      STEPX = (XF - XI)/(NP(1) -1)
      STEPY = (YF - YI)/(NP(2) -1)
      SHM(1) = -1.0/(12*STEPX**2) !0.0D+0!1.0/(12*STEPX**2)
      SHM(2) = -1.0/(12*STEPY**2) !0.0D+0!1.0/(12*STEPY**2)
      SHM(3) = -1.0/(48*STEPX*STEPY)

c-------

c      DO K = 1,N*N/2

CCC GENERATE FUNCTION VECTOR
      DO L=1,NP(1),1
         X = XI + (L-1)*STEPX
         DO M=1,NP(2),1
            Y = YI + (M-1)*STEPY
            I = L + (M-1)*NP(1)
            US(I) = GAUSS(X,Y,0) !A*DEXP(-B*x**2)*DEXP(-B*y**2)
            VPOT(I) = 0.0D+0
            FXY = GAUSS(X,Y,3)  
            FXX = GAUSS(X,Y,1)  
            FYY = GAUSS(X,Y,2)  
            ANS(I) =  FXX + FYY + FXY
         ENDDO
      ENDDO

c------  Generate Hamiltonian matrix
!     CALL AP2D(NP,AP,SHM,VPOT)
      CALL AP2DCT(NP,AP,SHM,VPOT)

cccccccccccccccccccccccccccccc PROPAGATION ROUTINE

      !TOTAL DXX + DYY + DXY  
      SHM(1) = -1.0/(12*STEPX**2)!0.0D+0!1.0/(12*STEPX**2)
      SHM(2) = -1.0/(12*STEPY**2)!0.0D+0!1.0/(12*STEPY**2)
      SHM(3) = -1.0/(48*STEPX*STEPY)!0.0D+0!-1.0/(48*STEPX*STEPY)
      CALL AU_2DCT(NP, SHM, VPOT, US, AU)
      
ccccccccccccccccccccccccccc

CCC   APPLY DERIVATIVE MATRIX

      DO J=1,N,1
         DO I=1,J,1
            APP(I,J) = AP(I+(J-1)*J/2)
            APP(J,I) = AP(I+(J-1)*J/2)
         ENDDO
      ENDDO


      DO I=1,N,1
         AU2(I)=0.0D+0
         DO J=1,N,1
            AU2(I) = AU2(I)+APP(I,J)*US(J)
            IF(ABS(APP(I,J)).GT.ZERO)THEN
               !AU2(I) = AU2(I)+APP(I,J)*US(J)
               WRITE(2,'(I6,I6,ES26.18,ES26.18,ES26.18)')I,J,APP(I,J)
     &              ,US(J),APP(I,J)*US(J)
            ENDIF
            !AU2(I) = AU2(I)+APP(I,J)*US(J)
         ENDDO
         WRITE(2,'(I6,I6,ES26.18)')I,J,AU2(I)
         WRITE(*,'(I6,I6,ES26.18,ES26.18,ES26.18,ES20.9)')I,J,
     &      AU2(I),AU(I),ANS(I),(AU(I)-AU2(I))/(AU(I))
         IF(ABS((AU(I)-AU2(I))/(AU(I))).EQ.0.0D+0)WRITE(3,*)I,
     &        ABS((AU(I)-AU2(I))/(AU(I)))
         IF(ABS((AU(I)-AU2(I))/(AU(I))).GT.5.0D-14) WRITE(1,*)I,
     &        ABS((AU(I)-AU2(I))/(AU(I)))
      ENDDO
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 1111 FORMAT(100(I3,1X))
      STOP
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
CCCCCCCCCCCCCCCCc TEST FUNCTION

      FUNCTION GAUSS(X,Y,I)
      implicit none
      INTEGER I
      REAL*8 A,B,X,Y,GAUSS
      
      A = 1.0D+0
      B = 1.5D+0

      IF(I.EQ.0)THEN
         IF(X.GT.0.5D+1 .OR. X.LT.-0.5D+1)THEN 
            GAUSS = 0.0D+0
         ELSEIF(Y.GT.0.3D+1 .OR. Y.LT.-0.3D+1)THEN 
            GAUSS = 0.0D+0
         ELSE
            GAUSS = A*DEXP(-B*x**2)*DEXP(-B*y**2)
         ENDIF
      ELSEIF(I.EQ.1)THEN !FXX
         GAUSS = (4*X*X*B -2)*B*A*DEXP(-B*x**2 - B*y**2)
      ELSEIF(I.EQ.2)THEN !FYY
         GAUSS = (4*Y*Y*B -2)*B*A*DEXP(-B*x**2 - B*y**2)
      ELSEIF(I.EQ.3)THEN !FXY
         GAUSS = 4*X*Y*A*B*B*DEXP(-B*x**2 - B*y**2)
      ENDIF

      RETURN
      END

      FUNCTION NUMDERIV(X,Y,I,STEPX,STEPY)
      implicit none
      INTEGER I
      REAL*8 X,Y,XX(10),YY(10),STEPX,STEPY,NUMDERIV,GAUSS

      IF(I.EQ.1)THEN !FXX
         XX(1) = X + 2.0*STEPX
         XX(2) = X + STEPX
         XX(3) = X 
         XX(4) = X - STEPX
         XX(5) = X - 2.0*STEPX
         NUMDERIV = (-GAUSS(XX(1),Y,0) +16.0*GAUSS(XX(2),Y,0)
     &      -30.0*GAUSS(XX(3),Y,0) +16.0*GAUSS(XX(4),Y,0) 
     &      -GAUSS(XX(5),Y,0))/(12.0*STEPX**2)
      ELSEIF(I.EQ.2)THEN !FYY
         YY(1) = Y + 2.0*STEPY
         YY(2) = Y + STEPY
         YY(3) = Y 
         YY(4) = Y - STEPY
         YY(5) = Y - 2.0*STEPY
         NUMDERIV = (-GAUSS(X,YY(1),0) +16.0*GAUSS(X,YY(2),0)
     &      -30.0*GAUSS(X,YY(3),0) +16.0*GAUSS(X,YY(4),0) 
     &      -GAUSS(X,YY(5),0))/(12.0*STEPY**2)
      ELSEIF(I.EQ.3)THEN !FXY
         XX(1) = X + 2.0*STEPX
         XX(2) = X + STEPX
         XX(3) = X - STEPX
         XX(4) = X - 2.0*STEPX
         YY(1) = Y + 2.0*STEPY
         YY(2) = Y + STEPY
         YY(3) = Y - STEPY
         YY(4) = Y - 2.0*STEPY
         NUMDERIV =(16*GAUSS(XX(2),YY(2),0) 
     &      + 16*GAUSS(XX(3),YY(3),0) -GAUSS(XX(1),YY(1),0)
     &      -GAUSS(XX(4),YY(4),0) -16*GAUSS(XX(2),YY(3),0)
     &      -16*GAUSS(XX(3),YY(2),0) + GAUSS(XX(1),YY(4),0)
     &      + GAUSS(XX(4),YY(1),0))/(48.0*STEPX*STEPY)
      ENDIF

      RETURN
      END
