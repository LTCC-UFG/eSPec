      FUNCTION FCORREL(N, U1, U2)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       N  
      REAL*8        U1(*), U2(*)
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
cdel      INTEGER  
cdel      REAL*8       
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
      REAL*8        ZERO
      PARAMETER     (ZERO = +0.0D+0)
c     **
c     ** Local scalars 
cdel       LOGICAL
cdel       CHARACTER*1
      INTEGER       I 
      REAL*8        SUM, FCORREL

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
c     ..
c     .. init and test lanczos parameters 
      SUM = ZERO
C
      DO I=1,N,1
         SUM = SUM + U1(I)*U2(I)
      ENDDO 
      FCORREL = SUM
C
      RETURN
      END
