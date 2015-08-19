c     subroutine AU
c     ==================
      SUBROUTINE AU(DIM, N, NP, SHM, VPOT, U, WORK, VAR)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      CHARACTER*(*) DIM
      INTEGER       N
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
      INTEGER       NP(*)
      REAL*8        SHM(*), VPOT(*), U(*), WORK(*), VAR(*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     |p^(2) = |p^(0) + 2*i*dt*|H* |p^(1)/hbar.
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
*     (21/03/2003) First version PSODAUX written by Freddy 
*
c     **
c     ** Parameters 
c      REAL*8        ZERO, ONE
c      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
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
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL      AU_1D, AU_2D, AU_2DC, AU_3D    
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      
D      print *, var(1)
D      read(*,*)
      IF(DIM(1:3).EQ.'.1D')THEN
         CALL AU_1D(N, SHM, VPOT, U, WORK)
      ELSEIF(DIM(1:5).EQ.'.2DCT')THEN
         CALL AU_2DCT(NP, SHM, VPOT, U, WORK)
      ELSEIF(DIM(1:4).EQ.'.2DC')THEN
         CALL AU_2DC(NP, SHM, VPOT, U, WORK)
      ELSEIF(DIM(1:4).EQ.'.2DT')THEN
         CALL AU_2DT(NP, VAR, SHM, VPOT, U, WORK)
      ELSEIF(DIM(1:4).EQ.'.2DL')THEN
         CALL AU_2DL(NP, VAR, SHM, VPOT, U, WORK)
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
         CALL AU_2D(NP, SHM, VPOT, U, WORK)
      ELSEIF(DIM(1:3).EQ.'.3D')THEN
         CALL AU_3D(NP, SHM, VPOT, U, WORK)
      ELSE
         WRITE(*,1001)
         STOP
      ENDIF
c     
 1001 FORMAT('<<<>>> Dimension error <<<>>>')
c     ..
      RETURN
      END
