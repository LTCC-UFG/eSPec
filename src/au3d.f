      SUBROUTINE AU_3D(NP, SHM, VPOT, US, AU)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*1      
cdel      INTEGER  
cdel      REAL*8              
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
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
*     (03/02/2003) First version AU_3D written by Freddy
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0E+0, ONE = +1.0E+0, SIXTEEN = +1.6E+1,
     &     THIRTY = +3.0E+1)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J, K, KKW1, IT01, IT02, IT03, IT04, IT05, IT06
      INTEGER       IT08, IT09, IT10, IT11, IT12, IT13
      REAL*8        SHT

c     **
c     ** External functions 
cdel      LOGICAL
cdel     CHARACTER*
cdel      INTEGER       
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program      
c     .. Initiate and test parameters 
      SHT = SHM(1) + SHM(2) + SHM(3)
c
      DO K=1,NP(3),1
         DO J=1,NP(2),1
            DO I=1,NP(1),1
               KKW1 = I + (J - 1)*NP(1) + (K - 1)*NP(1)*NP(2)
c
               IT01 = ONE
               IT02 = ONE
               IT03 = ONE
               IT04 = ONE
               IT05 = ONE
               IT06 = ONE
               IT08 = ONE
               IT09 = ONE
               IT10 = ONE
               IT11 = ONE
               IT12 = ONE
               IT13 = ONE
c
               IF(K.LE.2 .OR. KKW1-2*NP(1)*NP(2).LE.0) IT01 = ZERO
               IF(K.LE.1 .OR. KKW1-NP(1)*NP(2).LE.0) IT02 = ZERO
               IF(J.LE.2 .OR. KKW1-2*NP(1).LE.0) IT03 = ZERO
               IF(J.LE.1 .OR. KKW1-NP(1).LE.0) IT04 = ZERO
               IF(I.LE.2 .OR. KKW1-2.LE.0) IT05 = ZERO
               IF(I.LE.1 .OR. KKW1-1.LE.0) IT06 = ZERO
               IF(I.GE.NP(1) .OR. KKW1+1.GT.NP(1)*NP(2)*NP(3)) 
     &              IT08 = ZERO
               IF(I.GE.NP(1)-1 .OR. KKW1+2.GT.NP(1)*NP(2)*NP(3)) 
     &              IT09 = ZERO
               IF(J.GE.NP(2) .OR. KKW1+NP(1).GT.NP(1)*NP(2)*NP(3)) 
     &              IT10 = ZERO
               IF(J.GE.NP(2)-1 .OR. KKW1+2*NP(1).GT.NP(1)*NP(2)*NP(3)) 
     &              IT11 = ZERO
               IF(K.GE.NP(3) .OR. KKW1+NP(1)*NP(2).GT.NP(1)*NP(2)*NP(3)) 
     &              IT12 = ZERO
               IF(K.GE.NP(3)-1 .OR. KKW1+2*NP(1)*NP(2).GT.
     &              NP(1)*NP(2)*NP(3)) IT13 = ZERO
c               
               AU(KKW1) =
     &              + IT01*SHM(3)*US(KKW1-2*NP(1)*NP(2))
     &              - IT02*SIXTEEN*SHM(3)*US(KKW1-NP(1)*NP(2))
     &              + IT03*SHM(2)*US(KKW1-2*NP(1))
     &              - IT04*SIXTEEN*SHM(2)*US(KKW1-NP(1))
     &              + IT05*SHM(1)*US(KKW1-2)
     &              - IT06*SIXTEEN*SHM(1)*US(KKW1-1)
     &              + (THIRTY*SHT + VPOT(KKW1))*US(KKW1)
     &              - IT08*SIXTEEN*SHM(1)*US(KKW1+1)
     &              + IT09*SHM(1)*US(KKW1+2)
     &              - IT10*SIXTEEN*SHM(2)*US(KKW1+NP(1))
     &              + IT11*SHM(2)*US(KKW1+2*NP(1))
     &              - IT12*SIXTEEN*SHM(3)*US(KKW1+NP(1)*NP(2))
     &              + IT13*SHM(3)*US(KKW1+2*NP(1)*NP(2))
c
            ENDDO
         ENDDO
      ENDDO
c     
      RETURN
      END




