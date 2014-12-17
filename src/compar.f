      FUNCTION ICOMPAR(NAME1, NAME2)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL   
      CHARACTER*(*) NAME1, NAME2       
      INTEGER       ICOMPAR  
cdel      REAL*8      
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
*     (12/03/2003) First version ICOMPAR written by Freddy
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER
      INTEGER       N1, N2  
cdel      REAL*8 

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER
      INTEGER       ICHLENGTH     
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program      
      N1 = ICHLENGTH(NAME1, 5)
      N2 = ICHLENGTH(NAME2, 5)
      IF(NAME1(1:N1).EQ.NAME2(1:N2))THEN
         ICOMPAR = ZERO
      ELSE
         ICOMPAR = - ONE
      ENDIF
c
      RETURN
      END
