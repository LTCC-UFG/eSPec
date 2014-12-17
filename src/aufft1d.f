      program teste
      implicit real*8(a-h,o-z)
      integer       np(2048)
      real*8        tc(2048), vpot(2048), us(2048), au(2048)
      real*8        wk1(2048), wk2(2048)
c
      n = 1024
      ndim = 1
      np(1) = n
      np(2) = 1d0
      np(3) = 1d0
c
      do i=1,n,1
         tc(i) = i**2/1d5
         vpot(i) = 0d0
         us(i) = i
         au(i) = 0d0
         wk1(i)  = 0d0
         wk2(i) = 0d0
      enddo     
c
      call aufft_nd(n, ndim, np, tc, vpot, us, au, wk1, wk2)
c
      end

      SUBROUTINE AUFFT_ND(N, NDIM, NP, TC, VPOT, US, AU, WK1, WK2)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      INTEGER       N, NDIM
cdel      REAL*8       
c     **
c     ** Array arguments
      INTEGER       NP(*)
      REAL*8        TC(*), VPOT(*), US(*), AU(*), WK1(*), WK2(*)
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
*     (03/02/2003) First version AU_1D written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8          ZERO, ONE
      PARAMETER       (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** local scalars **
      INTEGER        I, NF
cdel      REAL*8         SHAUX
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program      
c paralelizar (estudar numero de processadores para otimizar a paral.)
      DO I=1,N,1
         WK1(I) = US(I)
         WK2(I) = ZERO
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,n,1
         write(1,*)i,wk1(i), wk2(i), au(i)
      enddo
cccccccccccccccccccccccccccccccccccccccccccc
c     
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(WK1, WK2, 2*N, 2*NP(I), 2*NF, 1)
      ENDDO
c
cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,2*n,1
         write(2,*)i,wk1(i), wk2(i), au(i)
      enddo
cccccccccccccccccccccccccccccccccccccccccccc
c     
      DO I=1,N,1
         WK1(I) = TC(I)*WK1(I) 
         WK2(I) = TC(I)*WK2(I)
         WK1(2*N-I+1) = TC(I)*WK1(2*N-I+1)
         WK2(2*N-I+1) = TC(I)*WK2(2*N-I+1)
      ENDDO
c
cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,2*n,1
         write(3,*)i,wk1(i), wk2(i), au(i)
      enddo
cccccccccccccccccccccccccccccccccccccccccccc      
c     
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(WK1, WK2, 2*N, 2*NP(I), 2*NF, -1)
      ENDDO
c     
cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,2*n,1
         write(4,*)i,wk1(i)/(2*N), wk2(i)/(2*N), au(i)
      enddo
cccccccccccccccccccccccccccccccccccccccccccc
c
      DO I=1,N,1
         AU(I) = WK1(I)/(2*N) + VPOT(I)*US(I)
c     AV(I) = 
      ENDDO
c
cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,2*n,1
         write(11,*)i, au(i), wk2(i)
      enddo
cccccccccccccccccccccccccccccccccccccccccccc
c     
      RETURN
      END

