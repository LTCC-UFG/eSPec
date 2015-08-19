      SUBROUTINE INIT_COND(DIM, INIEIGVC, NP, rK, rA, SH, XI, X0, U, V, 
     &     WORK1, WORK2, VAR)
      
c     Purpose:
c     Obtaining the initial conditions for the scattering wave function at t=0.
      
      IMPLICIT NONE
c     **
c     ** Scalar arguments
      CHARACTER*(*) DIM, INIEIGVC
c     **
c     **  Array arguments
      INTEGER       NP(*)
      REAL*8        SH(*), XI(*), U(*), V(*), WORK1(*), WORK2(*), X0(*)
      REAL*8        rK(*), rA(*), VAR(*)
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
c      REAL*8        TWOPI, FFREQ
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0) 
c     **
c     ** Local scalars
      INTEGER       I, J, K, iTEMP, IK
      REAL*8        COORD, COS_T, SIN_T, GAUSS_T, rNORM, COAUX
      REAL*8        rNORM1, rNORM2, TOL, SHIG
c     **
c     ** External functions
      REAL*8        GET_GAUSS
c      
      IF (DIM(1:3).EQ.'.1D' .AND. INIEIGVC(1:5).EQ.'.CALC') THEN
         PRINT*, ">>>>> vinicius"
         PRINT*, NP(1), X0(1), rK(1), rA(1)
         rNORM = ZERO    
         DO I = 1, NP(1), 1
            COORD = XI(1) + (I - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, ONE)
            rNORM = rNORM + GAUSS_T**2
         ENDDO
         rNORM = ONE/SQRT(rNORM)
c
         DO I = 1, NP(1), 1
            COORD = XI(1) + (I - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, rNORM)
c     
            COAUX = rK(1)*COORD
            U(I) =  DCOS(COAUX)*GAUSS_T
c           I think the sign is wrong. 
c            V(I) = -DSIN(COAUX)*GAUSS_T
c
            V(I) = DSIN(COAUX)*GAUSS_T
         ENDDO
c           
      ELSEIF(DIM(1:3).EQ.'.2D' .AND. INIEIGVC(1:6).EQ.'.GETR2')THEN
c         write(*,*)'tes',ra(2),rk(2),x0(2)
         DO I = 1, NP(2), 1    
            WORK1(I) = U(I)
         ENDDO
C
         rNORM = ZERO    
         DO I = 1, NP(1), 1
            COORD = XI(1) + (I - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, ONE)
            rNORM = rNORM + GAUSS_T**2
         ENDDO
         rNORM = ONE/SQRT(rNORM)
C
         DO J = 1, NP(1), 1   
            COORD = XI(1) + (J - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, rnorm)
            
            COAUX = rK(1)*COORD
            COS_T = DCOS(COAUX)*GAUSS_T
            SIN_T = DSIN(COAUX)*GAUSS_T
c            
            DO I = 1, NP(2), 1
               K = (I - 1)*NP(1) + J
               U(K) =  WORK1(I)*COS_T
               V(K) = -WORK1(I)*SIN_T
            ENDDO          
         ENDDO    
      ELSEIF(DIM(1:3).EQ.'.2D' .AND. 
     &        (INIEIGVC(1:5).EQ.'.GETR' .OR. 
     &        INIEIGVC(1:6).EQ.'.GETR1')) THEN
c         write(*,*)'tes',ra(2),rk(2),x0(2)
         DO I = 1, NP(1), 1    
            WORK1(I) = U(I)
         ENDDO
C
         rNORM = ZERO    
         DO I = 1, NP(2), 1
            COORD = XI(2) + (I - 1)*SH(2) - X0(2)
            GAUSS_T = GET_GAUSS(rA(2), COORD, ONE)
            rNORM = rNORM + GAUSS_T**2
         ENDDO
         rNORM = ONE/SQRT(rNORM)
C
         DO J = 1, NP(2), 1   
            COORD = XI(2) + (J - 1)*SH(2) - X0(2)
            GAUSS_T = GET_GAUSS(rA(2), COORD, rnorm)
            
            COAUX = rK(2) * COORD
            COS_T = DCOS(COAUX)*GAUSS_T
            SIN_T = DSIN(COAUX)*GAUSS_T
            
            iTEMP = (J - 1) * NP(1)
            DO I = 1, NP(1), 1
               K = iTEMP + I
               U(K) =  WORK1(I)*COS_T
               V(K) = -WORK1(I)*SIN_T
            ENDDO
            
         ENDDO
      ELSEIF(DIM(1:3).EQ.'.2D' .AND. INIEIGVC(1:6).EQ.'.GETC2') THEN   
         DO I = 1, NP(2), 1    
            WORK1(I) = U(I)
            WORK2(I) = V(I)
         ENDDO
c
         rNORM = ZERO    
         DO I = 1, NP(1), 1
            COORD = XI(1) + (I - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, ONE)
            rNORM = rNORM + GAUSS_T**2
         ENDDO
         rNORM = 1.0D0/SQRT(rNORM)
c
         DO J = 1, NP(1), 1   
            COORD = XI(1) + (J - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, rNORM)
c            
            COAUX = rK(1)*COORD
            COS_T = DCOS(COAUX)*GAUSS_T
            SIN_T = DSIN(COAUX)*GAUSS_T
c            
            DO I = 1, NP(2), 1
               K = (I - 1)*NP(1) + J
               U(K) =  WORK1(I)*COS_T + WORK2(I)*SIN_T
               V(K) = -WORK1(I)*SIN_T + WORK2(I)*COS_T 
            ENDDO    
         ENDDO
      ELSEIF(DIM(1:3).EQ.'.2D' .AND. 
     &        (INIEIGVC(1:5).EQ.'.GETC' .OR.
     &        INIEIGVC(1:6).EQ.'.GETC1'))THEN   
         DO I = 1, NP(1), 1    
            WORK1(I) = U(I)
            WORK2(I) = V(I)
         ENDDO
C
         rNORM = ZERO    
         DO I = 1, NP(2), 1
            COORD = XI(2) + (I - 1)*SH(2) - X0(2)
            GAUSS_T = GET_GAUSS(rA(2), COORD, ONE)
            rNORM = rNORM + GAUSS_T**2
         ENDDO
         rNORM = 1.0D0/SQRT(rNORM)
C
         DO J = 1, NP(2), 1   
            COORD = XI(2) + (J - 1)*SH(2) - X0(2)
            GAUSS_T = GET_GAUSS(rA(2), COORD, rNORM)
            
            COAUX = rK(2)*COORD
            COS_T = DCOS(COAUX)*GAUSS_T
            SIN_T = DSIN(COAUX)*GAUSS_T
            
            iTEMP = (J - 1) * NP(1)
            DO I = 1, NP(1), 1
               K = iTEMP + I
               U(K) =  WORK1(I)*COS_T + WORK2(I)*SIN_T
               V(K) = -WORK1(I)*SIN_T + WORK2(I)*COS_T 
            ENDDO    
         ENDDO
      ELSEIF(DIM(1:3).EQ.'.2D' .AND. INIEIGVC(1:5).EQ.'.CALC')THEN
d         write(*,*)'np1,np2,xi,xf',NP(1),np(2),xi(1),xi(2),sh(1),sh(2)
d         write(*,*)'rk',rk(1),rk(2)
d         read(*,*)
c     .. 
c     .. Compute the normalization constant of the gaussian function 1
         rNORM1 = ZERO    
         DO I = 1, NP(1), 1
            COORD = XI(1) + (I - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, ONE)
            rNORM1 = rNORM1 + GAUSS_T**2
         ENDDO
         rNORM1 = ONE/SQRT(rNORM1)
c     .. 
c     .. Compute the normalization constant of the gaussian function 2
         rNORM2 = ZERO   
         DO I = 1, NP(2), 1
            COORD = XI(2) + (I - 1)*SH(2) - X0(2)
            GAUSS_T = GET_GAUSS(rA(2), COORD, ONE)
            rNORM2 = rNORM2 + GAUSS_T**2
         ENDDO
         rNORM2 = ONE/SQRT(rNORM2)
c     ..
c     .. Compute the two dimensional gaussian wave packet
         DO I = 1, NP(1), 1 
            COORD = XI(1) + (I - 1)*SH(1) - X0(1)
            GAUSS_T = GET_GAUSS(rA(1), COORD, rNORM1)
c
            COAUX = rK(1)*COORD
            WORK1(I) =  DCOS(COAUX)*GAUSS_T
            WORK2(I) = -DSIN(COAUX)*GAUSS_T
         ENDDO
c      
         DO J = 1, NP(2), 1   
            COORD = XI(2) + (J - 1)*SH(2) - X0(2)
            GAUSS_T = GET_GAUSS(rA(2), COORD, rNORM2)
c            
            COAUX = rK(2)*COORD
            COS_T = DCOS(COAUX)*GAUSS_T
            SIN_T = DSIN(COAUX)*GAUSS_T
c            
            iTEMP = (J - 1) * NP(1)
            DO I = 1, NP(1), 1
               K = iTEMP + I
               U(K) =  WORK1(I)*COS_T + WORK2(I)*SIN_T
               V(K) = -WORK1(I)*SIN_T + WORK2(I)*COS_T 
            ENDDO
         ENDDO  
c

      ELSE
         WRITE(*,*) "<<<>>> Input file error! Please check the input fil
     &e and resubmit. <<<>>>"
         STOP
      ENDIF
c     ..
      IF(DIM(1:4).EQ.'.2DT')THEN
d      IF(DIM(1:3).EQ.'.2D')THEN
         TOL = 1.0D-99
         I = ZERO
 10      I = I + 1
         IK = VAR(1) + (I - ONE)*NP(1)
         write(*,*)'I,IK',I,IK,I*NP(1)
         DO K=IK,I*NP(1),1
            IF(ABS(U(K)).GT.TOL)TOL =  U(K)
            IF(ABS(V(K)).GT.TOL)TOL =  V(K)
         ENDDO
d         write(*,*)tol
         IF(K.LT.NP(1)*NP(2))GOTO 10
         WRITE(*,*)'   For the triangular grid is necessary to use a cut
     &off =', TOL, I,K
         WRITE(*,*)'   This means that all points of the initial conditi
     &on smaller' 
         WRITE(*,*)'      than',TOL,' in absolute value will be set iqua
     &l to zero'
         rNORM1 = ZERO
         DO K=1,NP(1)*NP(2),1
            IF(ABS(U(K)).LT.TOL)U(K) =  ZERO
            IF(U(K).GT.SHIG)SHIG =  U(K)
            IF(ABS(V(K)).LT.TOL)V(K) =  ZERO
            IF(V(K).GT.SHIG)SHIG =  V(K)
            rNORM1 = rNORM1 + U(K)*U(K) + V(K)*V(K)
         ENDDO
         WRITE(*,*)'   The biggest value in the grid is iqual to ',SHIG
c         IF(SHIG-TOL)
         rNORM1 = ONE/SQRT(rNORM1)
         DO K=1,NP(1)*NP(2),1
            U(K) = U(K)/rNORM1 
            V(K) = V(K)/rNORM1 
         ENDDO
      ENDIF
c     ..
c 1001 FORMAT(A167)
      RETURN
      END

