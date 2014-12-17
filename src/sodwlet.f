c     subroutine PSODFFTAUX
c     =====================
      SUBROUTINE PSODWLTAUX(NDIM, N, B, CST, NP, U1, U2, V1,  
     &     V2, TC, VPOT, WRKR, WRKI)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
cdel      LOGICAL
cdel      CHARACTER*(*) DIM
      INTEGER       NDIM, N
      REAL*8        B, CST
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*)
      REAL*8        U1(*), U2(*), V1(*), V2(*), TC(*), VPOT(*), WRKR(*)
      REAL*8        WRKI(*)
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
*     (21/03/2003) First version PSODWLTAUX written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, NF
      REAL*8        UAUX, VAUX 

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL      NWLET
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      DO I=1,N,1  
         WRKR(I) = U1(I) + CST*VPOT(I)*V2(I) ! WRKI = V*Psi
         WRKI(I) = V1(I) - CST*VPOT(I)*U2(I) ! WRKR = U*Psi
         U1(I) = U2(I)
         V1(I) = V2(I)
      ENDDO    
c
      DO I=1,N/2,1
         UAUX = U2(I)
         VAUX = V2(I)
         U2(I) = U2(I+N/2)
         V2(I) = V2(I+N/2)
         U2(I+N/2) = UAUX
         V2(I+N/2) = VAUX
      ENDDO
c
c      NF = ONE
c      DO I=1,NDIM,1
c         NF = NF*NP(I)
c         CALL NWLET(U2, V2, N, NP(I), NF, 1)
c      ENDDO
c
      DO I=1,N,2
         U2(I) = TC(I)*U2(I)
         U2(I+1) = TC(I+1)*U2(I+1)
         V2(I) = TC(I)*V2(I)
         V2(I+1) = TC(I+1)*V2(I+1)
      ENDDO
c
c      NF = ONE
c      DO I=1,NDIM,1
c         NF = NF*NP(I)
c         CALL FFT(U2, V2, N, NP(I), NF, -1)
c      ENDDO
c
      DO I=1,N/2,1
         UAUX = WRKR(I+N/2) + CST*B*V2(I)
         VAUX = WRKI(I+N/2) - CST*B*U2(I)
         U2(I) = WRKR(I) + CST*B*V2(I+N/2)
         V2(I) = WRKI(I) - CST*B*U2(I+N/2)
         U2(I+N/2) = UAUX
         V2(I+N/2) = VAUX
      ENDDO
c
      RETURN
      END
