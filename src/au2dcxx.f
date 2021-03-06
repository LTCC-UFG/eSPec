      SUBROUTINE AU_2DC(NP, SHM, VPOT, US, AU)
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
      INTEGER       KAUX, JAUX1, JAUX2, JAUX3, JAUX4, IAUX1
      INTEGER       IAUX2, IAUX3, IAUX4, NP1AUX, NP2AUX
      REAL*8        CF1, CF3, CF4, CF5, CF6, CF8, CF9
c     **
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER       I, J, K, KAUX
      REAL*8        CF(17)
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
      CF1 = SHM(3)/4.0D+0
      CF3 = -CF1
      CF4 = -SIXTEEN*CF1
      CF5 = -SIXTEEN*SHM(2)
      CF6 = -CF4
      CF8 = -SIXTEEN*SHM(1)
      CF9 = THIRTY*(SHM(1) + SHM(2))
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
         CF(2) = SHM(2)
         CF(5) = CF5
         CF(13) = CF(5)
         CF(16) = CF(2)
c
         IF(J.LE.2)CF(2) = ZERO
         IF(J.LE.1)CF(5) = ZERO
         IF(J.GE.NP(2))CF(13) = ZERO
         IF(J.GE.NP2AUX)CF(16) = ZERO

                     
         DO I=1,NP(1),1
c
            CF(1) = CF1
            CF(3) = CF3
            CF(4) = CF4
            CF(6) = CF6
            CF(7) = SHM(1)
            CF(8) = CF8
c            CF(9) = CF9
            CF(10) = CF8
            CF(11) = CF(7)
            CF(12) = CF6
            CF(14) = CF4
            CF(15) = CF3
            CF(17) = CF(1)
c
            IAUX4 = I - 2
            IAUX3 = I - 1
            IAUX2 = I + 1
            IAUX1 = I + 2 
c           
            K = I + KAUX
            IF(J.LE.2 .OR. I.LE.2)CF(1) = ZERO
            IF(J.LE.2 .OR. I.GE.NP1AUX)CF(3) = ZERO
            IF(J.LE.1 .OR. I.LE.1)CF(4) = ZERO
            IF(J.LE.1 .OR. I.GE.NP(1))CF(6) = ZERO
            IF(I.LE.2)CF(7) = ZERO
            IF(I.LE.1)CF(8) = ZERO
            IF(I.GE.NP(1))CF(10) = ZERO
            IF(I.GE.NP1AUX)CF(11) = ZERO
            IF(J.GE.NP(2) .OR. I.LE.1)CF(12) = ZERO
            IF(J.GE.NP(2) .OR. I.GE.NP(1))CF(14) = ZERO
            IF(J.GE.NP2AUX .OR. I.LE.2)CF(15) = ZERO
            IF(J.GE.NP2AUX .OR. I.GE.NP1AUX)CF(17) = ZERO
c
            AU(K) = 
     &             CF(1)*US(IAUX4+JAUX4) 
     &           + CF(2)*US(I+JAUX4) 
     &           + CF(3)*US(IAUX1+JAUX4)
c 
     &           + CF(4)*US(IAUX3+JAUX3) 
     &           + CF(5)*US(I+JAUX3) 
     &           + CF(6)*US(IAUX2+JAUX3)
c 
     &           + CF(7)*US(K-2) 
     &           + CF(8)*US(K-1)
     &           + (CF9 + VPOT(K))*US(K)  
     &           + CF(10)*US(K+1) 
     &           + CF(11)*US(K+2)
c
     &           + CF(12)*US(IAUX3+JAUX2) 
     &           + CF(13)*US(I+JAUX2) 
     &           + CF(14)*US(IAUX2+JAUX2)
c 
     &           + CF(15)*US(IAUX4+JAUX1) 
     &           + CF(16)*US(I+JAUX1) 
     &           + CF(17)*US(IAUX1+JAUX1)
         ENDDO
      ENDDO
c     ..
      RETURN
      END








