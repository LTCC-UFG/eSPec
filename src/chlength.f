      FUNCTION ICHLENGTH(NAME, MSP)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      CHARACTER*(*) NAME     
      INTEGER       MSP      
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
*     (10/03/2003) First version ICHLENGTH written by Freddy 
**
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, FOUR
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, FOUR =  +4.0D+0)
c     **
c     ** Local scalars 
      CHARACTER*1   LGTH
      INTEGER       I, J, ICHLENGTH
c     .. Start program
c     ..
      I = ZERO
 10   I = I + ONE
      LGTH = NAME(I:I)
      IF(LGTH.NE.' ')GOTO 10
      DO J=1,MSP,1
         I = I + ONE
         LGTH = NAME(I:I)
         IF(LGTH.NE.' ')GOTO 10
      ENDDO
c
      ICHLENGTH = I - (MSP + ONE)
c     .. End function
      RETURN
      END
