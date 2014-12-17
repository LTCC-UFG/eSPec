      SUBROUTINE HWHM(M, DT, TF, TAL, DE, U, V) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*(*) 
      INTEGER       M
      REAL*8        DT, TF, TAL, DE
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
cdel      INTEGER   
      REAL*8        U(*), V(*)       
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
*     (21/03/2003) First version SPECTRUM written by Freddy. 
*
c     **
c     ** Parameters 
      REAL*8        ONE, TWO, PI 
      PARAMETER     (ONE = +1.0D+0, TWO = +2.0D+0, 
     &     PI = +3.14159265358979323846D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, N
      REAL*8        CSTL, CS, SN, T, DEM, UI

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC
c     .. Start program      
      N = TWO**M
      DEM = - DE*1.5192676D0
      TAL = TAL*1.5192676D0
c      write(*,*)'***********DEM',DEM, DE, tal, n
c
      DO I=1,N,1
         T = -TF + (I - ONE)*DT
c         write(11,*)t,U(I),V(I)
         CS = DCOS(DEM*T)
         SN = DSIN(DEM*T)
         CSTL = DEXP(-ABS(T)*TAL)
         UI = U(I)
         U(I) = CS*U(I) - SN*V(I)
         V(I) = CS*V(I) + SN*UI
         U(I) = CSTL*U(I)
         V(I) = CSTL*V(I)
c         write(*,*)'oi2'
c         write(11,*)t,U(I),V(I)
      ENDDO
c      write(1,*)'#CSTL',cstl
c      DO I=N/2+1,N,1
c         T = -TF + (I - ONE)*DT
c         CS = DCOS(DEM*T)
c         SN = DSIN(DEM*T)
c         CSTL = DEXP(-PI*ABS(T)*TAL)
c         UI = U(I)
c         U(I) = CS*U(I) + SN*V(I)
c         V(I) = CS*V(I) - SN*UI
c         U(I) = CSTL*U(I)
c         V(I) = CSTL*V(I)
c      ENDDO
c
      RETURN
      END
