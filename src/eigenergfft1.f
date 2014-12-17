      FUNCTION EFFT2(NDIM, N, NP, B, B1, U, V, TC, VPOT)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
c      CHARACTER*(*) DIM
      INTEGER       NDIM, N
      REAL*8        B, B1
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
      INTEGER       NP(*)
      REAL*8        U(*), V(*), TC(*), VPOT(*)
c     **
*     ..
*     Purpose
*     =======
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
*     () First version written by 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
cdel      REAL*8 

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       I, NH, NF
      REAL*8        SUMU1, SUMU2, SUMV1, SUMV2, EFFT2, UAUX, VAUX
c     **
c     ** External subroutines 
      EXTERNAL      FFT
c     **
c     ** Intrinsic functions 
      INTRINSIC     NINT
c     .. Start program  
      NH = N/TWO
      SUMU1 = ZERO
      SUMV1 = ZERO
      DO I=1,N,1
         SUMU1 = SUMU1 + U(I)*VPOT(I)*U(I)
         SUMV1 = SUMV1 + V(I)*VPOT(I)*V(I)
      ENDDO
c
      DO I=1,NH,1
         UAUX = U(I)
         VAUX = V(I)
         U(I) = U(I+NH)
         V(I) = V(I+NH)
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
 
      SUMU2 = ZERO
      SUMV2 = ZERO
      DO I=1,N,1
         SUMU2 = SUMU2 + B1*U(I)*TC(I)*B1*U(I) 
         SUMV2 = SUMV2 + B1*V(I)*TC(I)*B1*V(I) 
      ENDDO
c
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(U, V, N, NP(I), NF, -1)
      ENDDO
c
      DO I=1,NH,1
         UAUX = B*U(I)
         VAUX = B*V(I)
         U(I) = B*U(I+NH)
         V(I) = B*V(I+NH)
         U(I+NH) = UAUX
         V(I+NH) = VAUX
      ENDDO
c
D      write(*,*)'oi2',SUMU1, SUMU2, SUMV1, SUMV2
      EFFT2 = SUMU1 + SUMU2 + SUMV1 + SUMV2
c
      RETURN
      END
