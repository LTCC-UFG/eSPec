      FUNCTION EFFT4(NDIM, N, NP, B, U, V, WU, WV, TC, VPOT)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*(*) DIM
      INTEGER       NDIM, N
      REAL*8        B
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
      INTEGER       NP(*)
      REAL*8        U(*), V(*), WU(*), WV(*), TC(*), VPOT(*)
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
      REAL*8        SUMU1, SUMU2, SUMV1, SUMV2, EFFT4
c     **
c     ** External subroutines 
      EXTERNAL      FFT
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
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
         WU(I) = U(I+NH)
         WV(I) = V(I+NH)
         WU(I+NH) = U(I)
         WV(I+NH) = V(I)
      ENDDO
c
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(WU, WV, N, NP(I), NF, 1)
      ENDDO
c      
      SUMU2 = ZERO
      SUMV2 = ZERO
      DO I=1,N,1
         WU(I) = TC(I)*WU(I) 
         WV(I) = TC(I)*WV(I) 
      ENDDO
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(WU, WV, N, NP(I), NF, -1)
      ENDDO
c
      SUMU2 = ZERO
      SUMV2 = ZERO
      DO I=1,NH,1
         SUMU2 = SUMU2 + U(I)*B*WU(I+NH)
         SUMV2 = SUMV2 + V(I)*B*WV(I+NH)
         SUMU2 = SUMU2 + U(I+NH)*B*WU(I)
         SUMV2 = SUMV2 + V(I+NH)*B*WV(I)
      ENDDO
c
      EFFT4 = SUMU1 + SUMU2 + SUMV1 + SUMV2
c
      RETURN
      END

