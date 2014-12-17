      SUBROUTINE DVTN(N, B, V1, V2)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       N
      REAL*8        B
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        V1(*), V2(*)
c     **
*     ..
*     Purpose
*     =======
*     v1 = A*v2
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
*     (01/2004) First version written by Freddy. 
*
c     **
c     ** Parameters 
cdel      REAL*8        ZERO, ONE
cdel      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I
cdel      REAL*8        
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
      DO I=1,N,1
         V2(I) = B*V1(I)
      ENDDO  
c
      RETURN
      END
