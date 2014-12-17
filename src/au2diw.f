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
      REAL*8          ONE, SIXTEEN, THIRTY
      PARAMETER       (ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = +3.0D+1) 
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J, K, NP1AUX, JAUX
      REAL*8        CF1, CF2, CF3

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
c     .. ********i=1
      CF1 = THIRTY*(SHM(1) + SHM(2))
      CF2 = -SIXTEEN*SHM(1)
      CF3 = -SIXTEEN*SHM(2)
      NP1AUX = 2*NP(1)
c
c     .. ******** j=1; jaux=0 
c     .. **** i=1   
      K = 1 
      AU(K) = (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c            AU(K) = SHM(2)*US(K-NP1AUX) 
c     &       + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
c     &       + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
c     &       + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=2
      K = 2
      AU(K) = CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=3...np(1)
      DO I=3,NP(1),1
         K = I 
         AU(K) = SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
      ENDDO
c
c     .. ********j=2; jaux=np(1) 
c     .. **** i=1...np(1)
      DO I=1,NP(1),1
         K = I + NP(1)
         AU(K) = CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
      ENDDO
c
c     .. ********j=3...n1-2  "main loop"
      DO J=3,NP(2)-2,1
         JAUX = (J - 1)*NP(1)
c     .. **** i=3...np(1)
         DO I=1,NP(1),1
            K = I + JAUX
            AU(K) = SHM(2)*US(K-NP1AUX) 
     &       + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &       + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &       + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
         ENDDO
      ENDDO
c
c     .. ********j=np(2)-1; jaux=np(2)*np(1) - 2*np(1)
      JAUX = NP(2)*NP(1) - 2*NP(1)
c     .. **** i=1...np(1)   
      DO I=1,NP(1),1
         K = I + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
      ENDDO

c     .. ********j=np(2)
      JAUX = NP(2)*NP(1) - NP(1)
c     .. **** i=1...np(1)-2         
      DO I=1,NP(1)-2,1
         K = I + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) 
      ENDDO
c     ..**** i=np(1)-1  
      K = NP(1) - 1 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
c     ..**** i=np(1)  
      K = NP(1) - 1 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
c     ..
      RETURN
      END




