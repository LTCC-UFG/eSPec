      FUNCTION EIGENERG(DIM, N, NP, SH, U1, VPOT, WORK, VAR)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) DIM
      INTEGER       N
cdel      REAL*8        AC
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
      INTEGER       NP(*)
      REAL*8        SH(*), U1(*), VPOT(*), WORK(*), VAR(*)
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
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
cdel      REAL*8 

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       I
      REAL*8        SUM, EIGENERG
c     **
c     ** External subroutines 
      EXTERNAL       AU
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
D      print *, var(1), DIM, N
D      print *, VPOT(1), U1(90), WORK(75)
D      read(*,*)
      CALL AU(DIM, N, NP, SH, VPOT, U1, WORK, VAR)
c
      SUM = ZERO
      DO I=1,N,1
         SUM = SUM + U1(I)*WORK(I)
      ENDDO 
      EIGENERG = SUM
c     ..
      RETURN
      END
