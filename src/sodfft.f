c     subroutine PSODFFTAUX
c     =====================
      SUBROUTINE PSODFFTAUX(ABSORB, TPABSOR, NDIM, N, B, CST, NP, U1, 
     &     U2, V1,  V2, TC, VPOT, VABC, WRKR, WRKI)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
      LOGICAL       ABSORB
      CHARACTER*(*) TPABSOR
      INTEGER       NDIM, N
      REAL*8        B, CST
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*)
      REAL*8        U1(*), U2(*), V1(*), V2(*), TC(*), VPOT(*), VABC(*)
      REAL*8        WRKR(*), WRKI(*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     |p^(2) = |p^(0) - 2*i*dt*|H* |p^(1)/hbar.
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
*     (21/03/2003) First version PSODAUX written by Freddy 
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
      REAL*8        UAUX, VAUX 

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
cdel      INTRINSIC     
c     .. Start program 
      NH = N/TWO
c
      IF(ABSORB .AND. TPABSOR(1:7).EQ.'.VOPTIC')THEN
         DO I=1,N,1  
            WRKR(I) = U1(I) + CST*VPOT(I)*V2(I) + CST*VABC(I)*U2(I)
            WRKI(I) = V1(I) - CST*VPOT(I)*U2(I) + CST*VABC(I)*V2(I)
            U1(I) = U2(I)
            V1(I) = V2(I)
         ENDDO
      ELSE
         DO I=1,N,1  
            WRKR(I) = U1(I) + CST*VPOT(I)*V2(I) ! WRKI = V*Psi
            WRKI(I) = V1(I) - CST*VPOT(I)*U2(I) ! WRKR = U*Psi
            U1(I) = U2(I)
            V1(I) = V2(I)
         ENDDO
      ENDIF
c
      DO I=1,NH,1
         UAUX = U2(I)
         VAUX = V2(I)
         U2(I) = U2(I+NH)
         V2(I) = V2(I+NH)
         U2(I+NH) = UAUX
         V2(I+NH) = VAUX
      ENDDO
c
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(U2, V2, N, NP(I), NF, 1)
      ENDDO
c
      DO I=1,N,1
         U2(I) = TC(I)*U2(I)
         V2(I) = TC(I)*V2(I)
      ENDDO
c
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(U2, V2, N, NP(I), NF, -1)
      ENDDO
c
      IF(ABSORB .AND. TPABSOR(1:7).EQ.'.SMOOTHW')THEN
         DO I=1,NH,1
            UAUX = WRKR(I+NH) + CST*B*V2(I)
            VAUX = WRKI(I+NH) - CST*B*U2(I)
            U2(I) = VABC(I)*(WRKR(I) + CST*B*V2(I+NH))
            V2(I) = VABC(I)*(WRKI(I) - CST*B*U2(I+NH))
            U2(I+NH) = VABC(I)*UAUX
            V2(I+NH) = VABC(I)*VAUX
         ENDDO
      ELSE
         DO I=1,NH,1
            UAUX = WRKR(I+NH) + CST*B*V2(I)
            VAUX = WRKI(I+NH) - CST*B*U2(I)
            U2(I) = WRKR(I) + CST*B*V2(I+NH)
            V2(I) = WRKI(I) - CST*B*U2(I+NH)
            U2(I+NH) = UAUX
            V2(I+NH) = VAUX
         ENDDO
      ENDIF
c     ..
      RETURN
      END
