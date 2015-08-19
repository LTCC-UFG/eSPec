      SUBROUTINE MTRXDIAG(DIM, IL, IU, INFO, MXDCT, N, ABSTOL, IWORK, 
     &     NP, EIGVL, SHM, VPOT, WORK, WK, EIGVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
c     del      LOGICAL       
      CHARACTER*(*) DIM
      INTEGER       IL, IU, INFO, MXDCT, N
      REAL*8        ABSTOL
c     **
c     ** Array arguments
c     del      LOGICAL       
c     del      CHARACTER*(*) 
      INTEGER       IWORK(*), NP(*)
      REAL*8        EIGVL(*), SHM(*), VPOT(*), WORK(*), WK(*)
      REAL*8        EIGVC(MXDCT,20)
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
*     (14/04/2003) First version MTRXDIAG written by Freddy.
*     (21/09/2014) Modified by Vinicius, created ap subroutines 
*      for organization sake, and fixed a small bug.

c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1)
c     **
c     ** Local scalars 
c     del      LOGICAL       
c     del      CHARACTER*1   
      INTEGER       I, J, M
      REAL*8        APAUX, SHT, VL, VU
c     **
c     ** Array scalars 
c     del      LOGICAL       
c     del      CHARACTER*1   
      INTEGER       IFAIL(50)
      REAL*8        AP(N*(N+1)/2)

c     ** External subroutines 
      EXTERNAL      DSPEVX
c     **
c     ** Intrinsic functions 
c     del      INTRINSIC     
c     .. Start program
      
c     * generate Hamiltonian matrix <<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
c     .. Change between one and two-dimension  
      IF(DIM(1:3).EQ.'.1D')THEN
c     .. Generate half simetrical 1-dimensional hamiltonian matrix 
         CALL AP1D(NP,AP,SHM,VPOT)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ELSEIF(DIM(1:5).EQ.'.2DCT')THEN
c     .. Generate half simetrical 2-dimensional hamiltonian matrix
c     .. with a d2/dxdy kinetic energy cross term 
         CALL AP2DCT(NP,AP,SHM,VPOT)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
c     .. Generate half simetrical 2-dimensional hamiltonian matrix 
         CALL AP2D(NP,AP,SHM,VPOT)
      ELSE
         WRITE(*,1011)
         WRITE(*,*)DIM
      ENDIF
c     .. Diagonalize hamiltonian half simetrical matrix  
      CALL  DSPEVX('V', 'I', 'U', N, AP, VL, VU, IL, IU, ABSTOL,
     &     M, EIGVL, EIGVC, MXDCT, WORK, IWORK, IFAIL, INFO)
c     ..
 1011 FORMAT('<<<>>> Matrix generation error <<<>>>')
 1111 FORMAT(1000(I3,1X))
c     ..
      RETURN
      END
