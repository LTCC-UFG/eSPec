      SUBROUTINE AU_2D(NP, SHM, VPOT, US, AU)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel     CHARACTER*1      
cdel      INTEGER  
cdel      REAL*8              
c     **
c     ** Array arguments
      INTEGER       NP(*) 
      REAL*8        SHM(*), VPOT(*), US(*), AU(*) 
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
*     (03/02/2003) First version AU_2D written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8          ZERO, SIXTEEN, THIRTY 
      PARAMETER       (ZERO = +0.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = +3.0D+1)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J, K
      INTEGER       KAUX, JAUX1, JAUX2, JAUX3, JAUX4
      INTEGER       NP1AUX, NP2AUX
      REAL*8        CF2, CF4
c     **
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER       I, J, K, KAUX
      REAL*8        CF(9)
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
cdel      INTRINSIC     
c     .. Start program
      CF2 = -SIXTEEN*SHM(2)
      CF4 = -SIXTEEN*SHM(1)
      CF(5) = THIRTY*(SHM(1) + SHM(2))
      NP1AUX = NP(1) - 1
      NP2AUX = NP(2) - 1
c
      DO J=1,NP(2),1
         KAUX = (J - 1)*NP(1)
         JAUX4 = (J - 3)*NP(1)
         JAUX3 = (J - 2)*NP(1)
         JAUX2 = J*NP(1)
         JAUX1 = (J + 1)*NP(1)
c
         CF(1) = SHM(2)
         CF(2) = CF2
         CF(8) = CF(2)
         CF(9) = CF(1)
c
         IF(J.LE.2)CF(1) = ZERO
         IF(J.LE.1)CF(2) = ZERO
         IF(J.GE.NP(2))CF(8) = ZERO
         IF(J.GE.NP2AUX)CF(9) = ZERO
c
         DO I=1,NP(1),1
c
            CF(3) = SHM(1)
            CF(4) = CF4
            CF(6) = CF4
            CF(7) = CF(3)
c
            K = I + KAUX
            IF(I.LE.2)CF(3) = ZERO
            IF(I.LE.1)CF(4) = ZERO
            IF(I.GE.NP(1))CF(6) = ZERO
            IF(I.GE.NP1AUX)CF(7) = ZERO
c
            AU(K) = CF(1)*US(I+JAUX4) 
     &           + CF(2)*US(I+JAUX3) 
     &           + CF(3)*US(K-2) 
     &           + CF(4)*US(K-1)
     &           + (CF(5) + VPOT(K))*US(K)  
     &           + CF(6)*US(K+1) 
     &           + CF(7)*US(K+2)
     &           + CF(8)*US(I+JAUX2) 
     &           + CF(9)*US(I+JAUX1) 
         ENDDO
      ENDDO
c     ..
      RETURN
      END








