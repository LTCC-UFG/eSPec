c     subroutine PSODAUX
c     ==================
      SUBROUTINE PSODAUX(ABSORB, DIM, N, CST, NP, SHM, U1, U2, V1, V2, 
     &     VPOT, VABC, WORK, VAR)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL
      CHARACTER*(*) DIM
      INTEGER       N
      REAL*8        CST
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
      INTEGER       NP(*)
      REAL*8        SHM(*), U1(*), U2(*), V1(*), V2(*), VPOT(*)
      REAL*8        VABC(*), WORK(*), VAR(*)
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
      INTEGER       I
      REAL*8        UAUX, VAUX 

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL      AU
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      CALL AU(DIM, N, NP, SHM, VPOT, V2, WORK, VAR)
c     
      IF(ABSORB)THEN
         DO I=1,N,1  
            UAUX = U2(I)
            U2(I) = U1(I) + CST*WORK(I) - CST*VABC(I)*V1(I) ! WORK = H*U
            U1(I) = UAUX
         ENDDO
      ELSE
         DO I=1,N,1  
            UAUX = U2(I)
            U2(I) = U1(I) + CST*WORK(I) ! WORK = H*U
            U1(I) = UAUX
         ENDDO
      ENDIF
c  
      CALL AU(DIM, N, NP, SHM, VPOT, U1, WORK, VAR)   
c
      DO I=1,N,1  
         VAUX = V2(I)
         V2(I) = V1(I) - CST*WORK(I) + CST*VABC(I)*V1(I) ! WORK = H*U
         V1(I) = VAUX
      ENDDO
c
      RETURN
      END
