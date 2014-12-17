      SUBROUTINE AU_2DT(NP, VAR, SHM, VPOT, US, AU)
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
      REAL*8        VAR(*), SHM(*), VPOT(*), US(*), AU(*) 
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
      INTEGER       I, J, K, NP1AUX, JAUX, IS
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
      CF1 = THIRTY*(SHM(1) + SHM(2))
      CF2 = -SIXTEEN*SHM(1)
      CF3 = -SIXTEEN*SHM(2)
      NP1AUX = 2*NP(1)
c      write(*,*)var(1)
c      read(*,*)
c     ..
c     .. ******** j=1; jaux=0 
      IS = NINT(1 + VAR(1))
c     .. **** i=1   
      K = 1 
      AU(K) = (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=2
      K = 2
      AU(K) = CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=3...np(1)
      DO I=3,IS-2,1
         K = I 
         AU(K) = SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
      ENDDO
c     .. **** i=np(1) - 1
      K =  IS - 1
      AU(K) = SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=np(1) 
      K =  IS  
      AU(K) = SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     ..
c     .. ********j=2; jaux=np(1)
      IS = NINT(2 + VAR(1))
c     .. **** i=1
      K = 1 + NP(1)
      AU(K) = CF3*US(K-NP(1)) 
     &     + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=2    
      K = 2 + NP(1)
      AU(K) = CF3*US(K-NP(1)) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=3...np(1)-2
      DO I=3,IS-2,1
         K = I + NP(1)
         AU(K) = CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
      ENDDO
c     .. **** i=np(1) - 1
      K = IS - 1 + NP(1)
      AU(K) = CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=np(1) 
      K = IS + NP(1)
CD      AU(K) = CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
CD     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K)  
CD     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
      AU(K) = SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K)  
     &     + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     ..
c     .. ********j=3...n2-v1-2  "1st main loop for np(2)"
      DO J=3,NP(2)-VAR(1)-2,1
         JAUX = (J - 1)*NP(1)
c     .. **** i=1
         K = 1 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) 
     &        + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=2
         K = 2 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     ..
         IS = NINT(J + VAR(1))
c     .. **** i=3...np(1)-2 "1st main loop np(1)"
         DO I=3,IS-2,1
            K = I + JAUX
            AU(K) = SHM(2)*US(K-NP1AUX) 
     &           + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &           + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &           + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
     &           + SHM(2)*US(K+NP1AUX)
         ENDDO
c     .. **** i=is-1
         K = IS - 1 + JAUX
CD         AU(K) = SHM(2)*US(K-NP1AUX) 
CD     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
CD     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
CD     &        + CF3*US(K+NP(1)) 
CD     &        + SHM(2)*US(K+NP1AUX)
         AU(K) = 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + CF3*US(K+NP(1)) 
     &        + SHM(2)*US(K+NP1AUX)
c     .. **** i=is
         K = IS + JAUX
         AU(K) = + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
     &        + CF3*US(K+NP(1)) 
     &        + SHM(2)*US(K+NP1AUX)
      ENDDO
c     ..
c     .. ********j=n2-v1-1
      JAUX = (NP(2) - VAR(1) - 2)*NP(1)
c     .. **** i=1
      K = 1 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) 
     &     + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=2
      K = 2 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=3...np(1)-2 "main loop "
      DO I=3,NP(1)-2,1
         K = I + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
     &        + SHM(2)*US(K+NP1AUX)
      ENDDO
c     .. **** i=np(1)-1
      K = NP(1) - 1 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + CF3*US(K+NP(1)) 
     &     + SHM(2)*US(K+NP1AUX)
c     .. **** i=np(1)
      K = NP(1) + JAUX
CD      AU(K) = SHM(2)*US(K-NP1AUX) 
CD     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
CD     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
CD     &     + CF3*US(K+NP(1)) 
CD     &     + SHM(2)*US(K+NP1AUX)
      AU(K) =  
     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
     &     + CF3*US(K+NP(1)) 
     &     + SHM(2)*US(K+NP1AUX)
c     ..
c     .. ********j=n2-v1...n1-2  "2nd main loop for NP(2)"
      DO J=NP(2)-VAR(1),NP(2)-2,1
         JAUX = (J - 1)*NP(1)
c     .. **** i=1
         K = 1 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) 
     &        + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=2
         K = 2 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) + SHM(2)*US(K+NP1AUX)
c     .. **** i=3...np(1)-2 "main loop "
         DO I=3,NP(1)-2,1
            K = I + JAUX
            AU(K) = SHM(2)*US(K-NP1AUX) 
     &           + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &           + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &           + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
     &           + SHM(2)*US(K+NP1AUX)
         ENDDO
c     .. **** i=np(1)-1
         K = NP(1) - 1 + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + CF3*US(K+NP(1)) 
     &        + SHM(2)*US(K+NP1AUX)
c     .. **** i=np(1)
         K = NP(1) + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
     &        + CF3*US(K+NP(1)) 
     &        + SHM(2)*US(K+NP1AUX)
      ENDDO
c     ..
c     .. ********j=np(2)-1; jaux=np(2)*np(1) - 2*np(1)
      JAUX = NP(2)*NP(1) - 2*NP(1)
c     .. **** i=1
      K = 1 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) 
     &     + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
c     .. **** i=2
      K = 2 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1))  
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
c     .. **** i=3...np(1)-2   
      DO I=3,NP(1)-2,1
         K = I + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) + CF3*US(K+NP(1)) 
      ENDDO
c     .. **** i=np(1)-1
      K = NP(1) - 1 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + CF3*US(K+NP(1)) 
c     .. **** i=np(1)
      K = NP(1) + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
     &     + CF3*US(K+NP(1)) 
c     ..
c     .. ********j=np(2)
      JAUX = NP(2)*NP(1) - NP(1)
c     .. **** i=1
      K = 1 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1))  
     &     + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) 
c     .. **** i=2
      K = 2 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1))  
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &     + SHM(1)*US(K+2) 
c     .. **** i=3...np(1)-2         
      DO I=3,NP(1)-2,1
         K = I + JAUX
         AU(K) = SHM(2)*US(K-NP1AUX) 
     &        + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &        + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
     &        + SHM(1)*US(K+2) 
      ENDDO
c     ..**** i=np(1)-1  
      K = NP(1) - 1 + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) + CF2*US(K+1) 
c     ..**** i=np(1)  
      K = NP(1) + JAUX
      AU(K) = SHM(2)*US(K-NP1AUX) 
     &     + CF3*US(K-NP(1)) + SHM(1)*US(K-2) 
     &     + CF2*US(K-1) + (CF1 + VPOT(K))*US(K) 
c     ..
d      do i=1,np(1)*np(2),1
d         write(13,*)i,au(i)
d      enddo
d      write(*,*)'finished!!!stop'
d      stop
d      read(*,*)
c     ..
      RETURN
      END




