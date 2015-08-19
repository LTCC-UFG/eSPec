      SUBROUTINE AU_2DCT(NP, SHM, VPOT, US, AU)
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
*     (14/09/2014) AU_2DCT written from AU_2D to include a kinetic cross term d2/dxdy
*     by Vinicius and Freddy
*
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
      INTEGER       I, J, K, NP1AUX
      REAL*8        CF1, CF2, CF3,CF4

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
      I = 1
      J = 1
      K = 1
      CF1 = THIRTY*(SHM(1) + SHM(2))
      CF2 = -SIXTEEN*SHM(1)
      CF3 = -SIXTEEN*SHM(2)
      CF4 = -SIXTEEN*SHM(3)
      NP1AUX = 2*NP(1)
c     
      AU(K) = (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &       + SHM(3)*US(K-NP1AUX-2)
c     &       - SHM(3)*US(K-NP1AUX+2)
c     &       + CF4*US(K-NP(1)-1)
c     &       - CF4*US(K-NP(1)+1)
c     &       - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
     &     + SHM(3)*US(K+NP1AUX+2)
c     
      J = 2
      K = I + (J - 1)*NP(1)
      AU(K) = + CF3*US(K-NP(1))
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &       + SHM(3)*US(K-NP1AUX-2)
c     &       - SHM(3)*US(K-NP1AUX+2)
c     &       + CF4*US(K-NP(1)-1)
     &       - CF4*US(K-NP(1)+1)
c     &       - CF4*US(K+NP(1)-1)
     &       + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
     &       + SHM(3)*US(K+NP1AUX+2)    
     
      DO J=3,NP(2)-2,1
         K = I + (J - 1)*NP(1)
         AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) + SHM(1)*US(K+2)
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &       + SHM(3)*US(K-NP1AUX-2)
     &       - SHM(3)*US(K-NP1AUX+2)
c     &       + CF4*US(K-NP(1)-1)
     &       - CF4*US(K-NP(1)+1)
c     &       - CF4*US(K+NP(1)-1)
     &       + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
     &       + SHM(3)*US(K+NP1AUX+2)        
      ENDDO
      J = NP(2) - 1
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
     &     + CF3*US(K+NP(1))
c     &       + SHM(3)*US(K-NP1AUX-2)
     &       - SHM(3)*US(K-NP1AUX+2)
c     &       + CF4*US(K-NP(1)-1)
     &       - CF4*US(K-NP(1)+1)
c     &       - CF4*US(K+NP(1)-1)
     &       + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
c     &       + SHM(3)*US(K+NP1AUX+2)         
      J = NP(2)
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
c     &       + SHM(3)*US(K-NP1AUX-2)
     &       - SHM(3)*US(K-NP1AUX+2)
c     &       + CF4*US(K-NP(1)-1)
     &       - CF4*US(K-NP(1)+1)
c     &       - CF4*US(K+NP(1)-1)
c     &       + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
c     &       + SHM(3)*US(K+NP1AUX+2)      
     
c     .. ********i=2
      I = 2
      J = 1
      K = I + (J - 1)*NP(1)
      AU(K) = + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
c     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
     &     + SHM(3)*US(K+NP1AUX+2)  
      J = 2
      K = I + (J - 1)*NP(1)
      AU(K) = + CF3*US(K-NP(1))
     &     + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
     &     + SHM(3)*US(K+NP1AUX+2)  
      DO J=3,NP(2)-2,1
         K=I+(J-1)*NP(1)
         AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &        + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) + SHM(1)*US(K+2)
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &     + SHM(3)*US(K-NP1AUX-2)
     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
     &     + SHM(3)*US(K+NP1AUX+2)  
      ENDDO
      J = NP(2) - 1
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) + SHM(1)*US(K+2)
     &     + CF3*US(K+NP(1)) 
c     &     + SHM(3)*US(K-NP1AUX-2)
     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
c
      J = NP(2)
      K = I +(J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1))
     &     + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) 
     &     + SHM(1)*US(K+2)
c     &     + SHM(3)*US(K-NP1AUX-2)
     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
c     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
c
c     .. ********i=3...n1-2      
      DO I=3,NP(1)-2,1
         J = 1
         K = I + (J - 1)*NP(1)
         AU(K) = SHM(1)*US(K-2) + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) + SHM(1)*US(K+2)
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &       + SHM(3)*US(K-NP1AUX-2)
c     &       - SHM(3)*US(K-NP1AUX+2)
c     &       + CF4*US(K-NP(1)-1)
c     &       - CF4*US(K-NP(1)+1)
     &       - CF4*US(K+NP(1)-1)
     &       + CF4*US(K+NP(1)+1)
     &       - SHM(3)*US(K+NP1AUX-2)
     &       + SHM(3)*US(K+NP1AUX+2)
         J = 2
         K = I + (J - 1)*NP(1)
         AU(K) = + CF3*US(K-NP(1))
     &        + SHM(1)*US(K-2) + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) + SHM(1)*US(K+2)
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &       + SHM(3)*US(K-NP1AUX-2)
c     &       - SHM(3)*US(K-NP1AUX+2)
     &       + CF4*US(K-NP(1)-1)
     &       - CF4*US(K-NP(1)+1)
     &       - CF4*US(K+NP(1)-1)
     &       + CF4*US(K+NP(1)+1)
     &       - SHM(3)*US(K+NP1AUX-2)
     &       + SHM(3)*US(K+NP1AUX+2)
d         write(1,'(A6, I6)')'I = ',I
d         write(1,'(A6,A6,A6,A6,A6,A6,A6,A6,A6)')
d     &   'K','-1','1','-16','16','16','-16','1','-1'
         DO J=3,NP(2)-2,1
            K = I + (J - 1)*NP(1)
            AU(K) = SHM(2)*US(K-NP1AUX) 
     &       + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &       + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &       + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
     &       + SHM(3)*US(K-NP1AUX-2)
     &       - SHM(3)*US(K-NP1AUX+2)
     &       + CF4*US(K-NP(1)-1)
     &       - CF4*US(K-NP(1)+1)
     &       - CF4*US(K+NP(1)-1)
     &       + CF4*US(K+NP(1)+1)
     &       - SHM(3)*US(K+NP1AUX-2)
     &       + SHM(3)*US(K+NP1AUX+2)
d            write(1,'(I6,I6, I6, I6, I6, I6, I6, I6, I6)') 
d     &      K,K-NP1AUX-2, K-NP1AUX+2,K-NP(1)-1,K-NP(1)+1,K+NP(1)-1,
d     &      K+NP(1)+1,K+NP1AUX-2,K+NP1AUX+2
         ENDDO
         J = NP(2) - 1
         K = I +(J - 1)*NP(1)
         AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &        + SHM(1)*US(K-2) + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) + SHM(1)*US(K+2)
     &        + CF3*US(K+NP(1))
     &        + SHM(3)*US(K-NP1AUX-2)
     &        - SHM(3)*US(K-NP1AUX+2)
     &        + CF4*US(K-NP(1)-1)
     &        - CF4*US(K-NP(1)+1)
     &        - CF4*US(K+NP(1)-1)
     &        + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
c     &       + SHM(3)*US(K+NP1AUX+2)
         J = NP(2)
         K = I +(J - 1)*NP(1)
         AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &        + SHM(1)*US(K-2) + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) + SHM(1)*US(K+2)
     &        + SHM(3)*US(K-NP1AUX-2)
     &        - SHM(3)*US(K-NP1AUX+2)
     &        + CF4*US(K-NP(1)-1)
     &        - CF4*US(K-NP(1)+1)
c     &        - CF4*US(K+NP(1)-1)
c     &        + CF4*US(K+NP(1)+1)
c     &       - SHM(3)*US(K+NP1AUX-2)
c     &       + SHM(3)*US(K+NP1AUX+2)
         
      ENDDO

cbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
c     .. ********i=n1-1
      I = NP(1) - 1
      J = 1
      K = I +(J - 1)*NP(1)
      AU(K) = SHM(1)*US(K-2) + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) 
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
c     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)
      J = 2
      K = I + (J - 1)*NP(1)
      AU(K) = + CF3*US(K-NP(1))
     &     + SHM(1)*US(K-2) 
     &     + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) 
     &     + CF3*US(K+NP(1)) 
     &     + SHM(2)*US(K+NP1AUX)
c     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)
      DO J=3,NP(2)-2,1      
         K = I + (J - 1)*NP(1)
         AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &        + SHM(1)*US(K-2) + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF2*US(K+1) 
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
     &        + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &        + CF4*US(K-NP(1)-1)
     &        - CF4*US(K-NP(1)+1)
     &        - CF4*US(K+NP(1)-1)
     &        + CF4*US(K+NP(1)+1)
     &        - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)
      ENDDO
      J = NP(2) - 1
      K = I + (J - 1)*NP(1)         
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + SHM(1)*US(K-2) + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1) 
     &     + CF3*US(K+NP(1)) 
     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
      J = NP(2)
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + SHM(1)*US(K-2) + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF2*US(K+1)
     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
     &     - CF4*US(K-NP(1)+1)
c     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + S'HM(3)*US(K+NP1AUX+2)  
c
c     .. ********i=n1
      I = NP(1)
      J = 1
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(1)*US(K-2) + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c      &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
c     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
      J = 2
      K = I + (J - 1)*NP(1)
      AU(K) = + CF3*US(K-NP(1))
     &     + SHM(1)*US(K-2) 
     &     + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF3*US(K+NP(1)) 
     &     + SHM(2)*US(K+NP1AUX)
c      &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
      DO J=3,NP(2)-2,1  
         K=I+(J-1)*NP(1)
         AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &        + SHM(1)*US(K-2) + CF2*US(K-1)
     &        + (CF1 + VPOT(K))*US(K)
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
      ENDDO
      J = NP(2) - 1
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + SHM(1)*US(K-2) + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + CF3*US(K+NP(1))
     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
      J = NP(2)
      K = I + (J - 1)*NP(1)
      AU(K) = SHM(2)*US(K-NP1AUX) + CF3*US(K-NP(1))
     &     + SHM(1)*US(K-2) + CF2*US(K-1)
     &     + (CF1 + VPOT(K))*US(K)
     &     + SHM(3)*US(K-NP1AUX-2)
c     &     - SHM(3)*US(K-NP1AUX+2)
     &     + CF4*US(K-NP(1)-1)
c     &     - CF4*US(K-NP(1)+1)
c     &     - CF4*US(K+NP(1)-1)
c     &     + CF4*US(K+NP(1)+1)
c     &     - SHM(3)*US(K+NP1AUX-2)
c     &     + SHM(3)*US(K+NP1AUX+2)  
c  
d      do i=1,np(1)*np(2),1
d         write(12,*)i,au(i)
d      enddo 
d      write(*,*)'finished!!!stop'
d      stop
c   
      RETURN
      END




