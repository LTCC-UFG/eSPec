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
      REAL*8          ONE, TWO
      PARAMETER       (ONE = +1.0E+0, TWO = +1.0E+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J, K, KAUX, IAUX1, IAUX2, JAUX, JAUX1,JAUX2
cdel      REAL*8    
c     **
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER       I, J, K, KAUX
      REAL*8        CF(4)
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
c      DATA CF / -0.577191652D+0, 0.988205891D+0, -0.897056937D+0, 
c     &     0.918206857D+0 /
c      SHM(I) = ONE/(2.4D+1*SM(I)*SH(I)*SH(I))
      CF(1) = - 6.0D+0*SHM(3)
      CF(2) = - 1.2D+1*SHM(1)
      CF(3) = - 1.2D+1*SHM(2)
      CF(4) = - 2.0D+0*(CF(2) + CF(3)) 
c     .. *****j=1
c     .. ********i=1
      AU(1) =  (CF(4) + VPOT(1))*US(1) + CF(2)*US(2) 
     &     + CF(3)*US(1+NP(1)) + CF(1)*US(2+NP(1)) 
c     .. ********i=2 until n-1
      DO I=2,NP(1)-1,1
         AU(I) = CF(2)*US(I-1) + (CF(4) + VPOT(I))*US(I) 
     &        + CF(2)*US(I+1) - CF(1)*US(I-1+NP(1)) 
     &        + CF(3)*US(I+NP(1)) + CF(1)*US(I+1+NP(1)) 
      ENDDO
c     .. ********i=n
      AU(NP(1)) = CF(2)*US(NP(1)-1) + (CF(4) + VPOT(NP(1)))*US(NP(1)) 
     &     + CF(3)*US(NP(1)+NP(1)) - CF(1)*US(NP(1)-1+NP(1))   
c
c     .. *****j=2 until n-1
      DO J=2,NP(2)-1,1
c     .. ********i=1
         K = ONE + (J - ONE)*NP(1)
         AU(K) = CF(3)*US(1+(J-2)*NP(1)) 
     &        + (CF(4) + VPOT(1+(J-1)*NP(1)))
     &                    *US(1+(J-1)*NP(1)) 
     &        + CF(2)*US(1+1+(J-1)*NP(1)) + CF(3)*US(1+J*NP(1)) 
     &        + CF(1)*US(2+J*NP(1)) - CF(1)*US(2+(J-2)*NP(1))
c     .. ********i=2 until n-1
         KAUX = (J - ONE)*NP(1)
         JAUX = J*NP(1)
         JAUX1 = (J - 1)*NP(1)
         JAUX2 = (J - 2)*NP(1)
         DO I=2,NP(1)-1,1 
            K = I + KAUX
            IAUX1 = I - 1
            IAUX2 = I + 1
            AU(K) = CF(1)*US(IAUX1+JAUX2) 
     &           + CF(2)*US(IAUX1+JAUX1) 
     &           + CF(3)*US(I+JAUX2) 
     &           + (CF(4) + VPOT(K))*US(K) 
     &           + CF(2)*US(IAUX2+JAUX1) 
     &           + CF(3)*US(I+JAUX) 
     &           + CF(1)*US(IAUX2+JAUX) 
     &           - CF(1)*US(IAUX2+JAUX2)
     &           - CF(1)*US(IAUX1+JAUX)  
         ENDDO
c     .. ********i=n 
         K = NP(1) + (J - ONE)*NP(1)
         AU(K) = CF(1)*US(NP(1)-1+(J-2)*NP(1)) 
     &        + CF(2)*US(NP(1)-1+(J-1)*NP(1)) 
     &        + CF(3)*US(NP(1)+(J-2)*NP(1)) 
     &        + (CF(4) + VPOT(NP(1)+(J-1)*NP(1)))
     &                    *US(NP(1)+(J-1)*NP(1)) 
     &        + CF(3)*US(NP(1)+J*NP(1)) 
     &        - CF(1)*US(NP(1)-1+J*NP(1))   
      ENDDO
c     .. *****j=n
c     .. ********i=1
      K = ONE + (NP(2) - ONE)*NP(1)
      AU(K) = CF(3)*US(1+(NP(2)-2)*NP(1)) 
     &     + (CF(4) + VPOT(1+(NP(2)-1)*NP(1)))
     &                 *US(1+(NP(2)-1)*NP(1)) 
     &     + CF(2)*US(2+(NP(2)-1)*NP(1)) - CF(1)*US(2+(NP(2)-2)*NP(1))
c     .. ********i=2 until n-1
      KAUX = (NP(2) - ONE)*NP(1)
      DO I=2,NP(1)-1,1
         K = I + KAUX
         AU(K) = CF(1)*US(I-1+(NP(2)-2)*NP(1)) 
     &        + CF(2)*US(I-1+(NP(2)-1)*NP(1)) 
     &        + CF(3)*US(I+(NP(2)-2)*NP(1)) 
     &        + (CF(4) + VPOT(I+(NP(2)-1)*NP(1)))
     &                    *US(I+(NP(2)-1)*NP(1)) 
     &        + CF(2)*US(I+1+(NP(2)-1)*NP(1)) 
     &        - CF(1)*US(I+1+(NP(2)-2)*NP(1))
      ENDDO
c     .. ********i=n
      K = NP(2)*NP(1) 
      AU(K) = CF(1)*US(NP(1)-1+(NP(2)-2)*NP(1)) 
     &     + CF(2)*US(NP(1)-1+(NP(2)-1)*NP(1)) 
     &     + CF(3)*US(NP(1)+(NP(2)-2)*NP(1)) 
     &     + (CF(4) + VPOT(NP(1)+(NP(2)-1)*NP(1)))
     &                 *US(NP(1)+(NP(2)-1)*NP(1)) 
c     ..
      RETURN
      END
c






