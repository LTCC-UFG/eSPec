      FUNCTION ECNORM(N, US)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*1      
      INTEGER       N
cdel      REAL*8      
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)  
cdel      INTEGER  
      REAL*8        US(*)  
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
*     (03/02/2003) First version ECNORM written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO
      PARAMETER     (ZERO = +0.0E+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I  
      REAL*8        SUM, ECNORM

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
      INTRINSIC     SQRT
c     .. Start program      
c     .. init and test parameters 
      SUM = ZERO
c 
      DO I=1,N,1
         SUM = SUM + US(I)*US(I)
      ENDDO
c     
      ECNORM = SQRT(SUM)
c     
      RETURN
      END
