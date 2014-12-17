c     .. random number generator
c     ==========================
      FUNCTION RAN(NSEED) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      INTEGER       NSEED
c     **
*     ..
*     Purpose
*     =======
*     Generate random numbers.
*
*     ..
*     Arguments
*     =========
*     NSEED 
*
*     ..
*     Authors
*     =======
*     
*
*     ..
*     Historic
*     ========
*     
c     **
c     ** Parameters 
      INTEGER       L, NC, M
      REAL*8        ONE
      PARAMETER     (L = 1029, NC = 221591, M = 1048576, ONE = +1.0E+0)
c     **
c     ** Local scalars 
      REAL*8        RAN
c
c     .. Start program            
      NSEED = MOD(NSEED*L + NC, M)
      RAN = ONE*NSEED/M
c
      RETURN
      END
