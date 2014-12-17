      SUBROUTINE MTRXDIAG(DIM, IL, IU, INFO, MXDCT, N, ABSTOL, IWORK, 
     &     NP, EIGVL, SHM, VPOT, WORK, WK, EIGVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) DIM
      INTEGER       IL, IU, INFO, MXDCT, N
      REAL*8        ABSTOL
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
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
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I, J, M
      REAL*8        APAUX, SHT, VL, VU
c     **
c     ** Array scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       IFAIL(50)
      REAL*8        AP(N*(N+1)/2)

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
cdel      REAL*8        
c     **
c     ** External subroutines 
      EXTERNAL      DSPEVX
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program
      
c   * generate Hamiltonian matrix <<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<
c     .. Change between one and two-dimension  
      IF(DIM(1:3).EQ.'.1D')THEN
c     .. Generate half simetrical 1-dimensional hamiltonian matrix 
         DO J=1,N,1
            DO I=1,J,1
               IF(I.EQ.J)THEN
                  APAUX = THIRTY*SHM(1) + VPOT(I)
               ELSEIF(I.EQ.J-1)THEN
                  APAUX = -SIXTEEN*SHM(1)
               ELSEIF(I.EQ.J-2)THEN
                  APAUX = ONE*SHM(1)
               ELSE
                  APAUX = ZERO
               ENDIF
               AP(I+(J-1)*J/2) = APAUX
            ENDDO
         ENDDO
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
c     .. Generate half simetrical 2-dimensional hamiltonian matrix 
         SHT = SHM(1) + SHM(2)
         DO J=1,N,1
            DO I=1,J,1
d               write(*,*)'This option does not work'
d               stop
               IF(I.EQ.J-2*NP(1))THEN
                  APAUX = ONE*SHM(2) 
               ELSEIF(I.EQ.J-NP(1))THEN
                  APAUX = -SIXTEEN*SHM(2) 
               ELSEIF(I.EQ.J-2 .AND. I.GT.2)THEN
                  APAUX = ONE*SHM(1)
               ELSEIF(I.EQ.J-1 .AND. I.GT.1)THEN
                  APAUX = -SIXTEEN*SHM(1) 
               ELSEIF(I.EQ.J)THEN
                  APAUX = THIRTY*SHT + VPOT(I)
               ELSE
                  APAUX = ZERO
               ENDIF

d               IF(I.EQ.J-NP(1)-2)THEN
d                  APAUX = ONE*SHM(2) 
d               ELSEIF(I.EQ.J-NP(1)-1)THEN
d                  APAUX = -SIXTEEN*SHM(2) 
d               ELSEIF(I.EQ.J-2)THEN
d                  APAUX = ONE*SHM(1)
d               ELSEIF(I.EQ.J-1)THEN
d                  APAUX = -SIXTEEN*SHM(1) 
d               ELSEIF(I.EQ.J)THEN
d                  APAUX = THIRTY*SHT + VPOT(I)
d               ELSE
d                  APAUX = ZERO
d               ENDIF
               AP(I+(J-1)*J/2) = APAUX 
            ENDDO
c            write(1,1111)(int(ap(i+(j-1)*j/2)), i=1,j,1)
         ENDDO
      ELSE
         WRITE(*,1011)
         WRITE(*,*)DIM
      ENDIF
c     .. Diagonalize hamiltonian half simetrical matrix  
      CALL  DSPEVX('V', 'I', 'U', N, AP, VL, VU, IL, IU, ABSTOL,
     &     M, EIGVL, EIGVC, MXDCT, WORK, IWORK, IFAIL, INFO)
c     ..
 1011 FORMAT('<<<>>> Matrix generation error <<<>>>')
 1111 FORMAT(100(I3,1X))
c     ..
      RETURN
      END
