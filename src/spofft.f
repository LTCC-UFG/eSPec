c     subroutine PSPOFFTAUX
c     =====================
      SUBROUTINE PSPOFFTAUX(ABSORB, TPABSOR, NDIM, N, B, B1, B2, CST1, 
     &     CST2, NP, U, V, TC, VPOT, VABC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
      LOGICAL       ABSORB
      CHARACTER*(*) TPABSOR
      INTEGER       NDIM, N
      REAL*8        B, B1, B2, CST1, CST2
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*)
      REAL*8        U(*), V(*), TC(*), VPOT(*), VABC(*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     |p^(2) = exp(-i*dt/hbar*V/2)*exp(-i*dt/hbar*T)*exp(-i*dt/hbar*V/2)|p^(1).
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
*     (02/2004) First version PSPOAUX written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, NF, NH
      REAL*8        UAUX, VAUX, TH1, TH2, CTH1, CTH2, STH1, STH2
      REAL*8        TH11, TH22

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL      FFT
c     **
c     ** Intrinsic functions 
      INTRINSIC     DSIN, DCOS, DEXP
c     .. Start program 
      NH = N/TWO
      IF(ABSORB .AND. TPABSOR(1:7).EQ.'.VOPTIC')THEN
         DO I=1,NH,1
            TH1 = CST2*VPOT(I)
            TH11 = DEXP(CST2*VABC(I))
            CTH1 = DCOS(TH1)*TH11
            STH1 = DSIN(TH1)*TH11
            UAUX = CTH1*U(I) + STH1*V(I)
            VAUX = CTH1*V(I) - STH1*U(I)
            TH2 = CST2*VPOT(I+NH)
            TH22 = DEXP(CST2*VABC(I+NH))
            CTH2 = DCOS(TH2)*TH22
            STH2 = DSIN(TH2)*TH22
            U(I) = CTH2*U(I+NH) + STH2*V(I+NH)
            V(I) = CTH2*V(I+NH) - STH2*U(I+NH)
            U(I+NH) = UAUX
            V(I+NH) = VAUX
         ENDDO
c
         NF = ONE
         DO I=1,NDIM,1
            NF = NF*NP(I)
            CALL FFT(U, V, N, NP(I), NF, 1)
         ENDDO
c
         DO I=1,N,1
            TH1 = CST1*TC(I)
            CTH1 = DCOS(TH1)
            STH1 = DSIN(TH1)
            U(I) = CTH1*U(I) + STH1*V(I)
            V(I) = CTH1*V(I) - STH1*U(I)
         ENDDO
c
         NF = ONE
         DO I=1,NDIM,1
            NF = NF*NP(I)
            CALL FFT(U, V, N, NP(I), NF, -1)
         ENDDO
c
         DO I=1,NH,1
            TH1 = CST2*VPOT(I+NH)
            TH11 = DEXP(CST2*VABC(I+NH))
            CTH1 = DCOS(TH1)*TH11
            STH1 = DSIN(TH1)*TH11
            UAUX = CTH1*B*U(I) + STH1*B*V(I)
            VAUX = CTH1*B*V(I) - STH1*B*U(I)
            TH2 = CST2*VPOT(I)
            TH22 = DEXP(CST2*VABC(I))
            CTH2 = DCOS(TH2)*TH22
            STH2 = DSIN(TH2)*TH22
            U(I) = CTH2*B*U(I+NH) + STH2*B*V(I+NH)
            V(I) = CTH2*B*V(I+NH) - STH2*B*U(I+NH)
            U(I+NH) = UAUX
            V(I+NH) = VAUX
         ENDDO
      ELSE
         DO I=1,NH,1
            TH1 = CST2*VPOT(I)
            CTH1 = DCOS(TH1)
            STH1 = DSIN(TH1)
            UAUX = CTH1*U(I) + STH1*V(I)
            VAUX = CTH1*V(I) - STH1*U(I)
            TH2 = CST2*VPOT(I+NH)
            CTH2 = DCOS(TH2)
            STH2 = DSIN(TH2)
            U(I) = CTH2*U(I+NH) + STH2*V(I+NH)
            V(I) = CTH2*V(I+NH) - STH2*U(I+NH)
            U(I+NH) = UAUX
            V(I+NH) = VAUX
         ENDDO
c         DO I=1,N,1
c            TH1 = CST2*VPOT(I)
c            CTH1 = DCOS(TH1)
c            STH1 = DSIN(TH1)
c            U(I) = CTH1*U(I) + STH1*V(I)
c            V(I) = CTH1*V(I) - STH1*U(I)
c         ENDDO         
c     
         NF = ONE
         DO I=1,NDIM,1
            NF = NF*NP(I)
            CALL FFT(U, V, N, NP(I), NF, 1)
         ENDDO
c     
         DO I=1,N,1
            TH1 = CST1*TC(I)
            CTH1 = DCOS(TH1)
            STH1 = DSIN(TH1)
            U(I) = CTH1*U(I) + STH1*V(I)
            V(I) = CTH1*V(I) - STH1*U(I)
         ENDDO
c     
         NF = ONE
         DO I=1,NDIM,1
            NF = NF*NP(I)
            CALL FFT(U, V, N, NP(I), NF, -1)
         ENDDO
c     
         IF(ABSORB .AND. TPABSOR(1:8).EQ.'.SMOOTHW')THEN
            DO I=1,NH,1
               TH1 = CST2*VPOT(I+NH)
               CTH1 = DCOS(TH1)
               STH1 = DSIN(TH1)
               UAUX = CTH1*B*U(I) + STH1*B*V(I)
               VAUX = CTH1*B*V(I) - STH1*B*U(I)
               TH2 = CST2*VPOT(I)
               CTH2 = DCOS(TH2)
               STH2 = DSIN(TH2)
               U(I) = VABC(I)*(CTH2*B*U(I+NH) + STH2*B*V(I+NH))
               V(I) = VABC(I)*(CTH2*B*V(I+NH) - STH2*B*U(I+NH))
               U(I+NH) = VABC(I+NH)*UAUX
               V(I+NH) = VABC(I+NH)*VAUX
            ENDDO
         ELSE
            DO I=1,NH,1
               TH1 = CST2*VPOT(I+NH)
               CTH1 = DCOS(TH1)
               STH1 = DSIN(TH1)
               UAUX = CTH1*B*U(I) + STH1*B*V(I)
               VAUX = CTH1*B*V(I) - STH1*B*U(I)
               TH2 = CST2*VPOT(I)
               CTH2 = DCOS(TH2)
               STH2 = DSIN(TH2)
               U(I) = CTH2*B*U(I+NH) + STH2*B*V(I+NH)
               V(I) = CTH2*B*V(I+NH) - STH2*B*U(I+NH)
               U(I+NH) = UAUX
               V(I+NH) = VAUX
            ENDDO
c            DO I=1,N,1
c               TH1 = CST2*VPOT(I+NH)
c               CTH1 = DCOS(TH1)
c               STH1 = DSIN(TH1)
c               U(I) = CTH1*B*U(I) + STH1*B*V(I)
c               V(I) = CTH1*B*V(I) - STH1*B*U(I)
c            ENDDO           
         ENDIF
      ENDIF
c     ..
      RETURN
      END
