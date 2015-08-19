      SUBROUTINE TRANSDIPOLE(IIL, IIU, IFL, IFU, LDZ, NT, VI, EIGVL, 
     &     EIGVC,DM)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       IIL, IIU, IFL, IFU, LDZ, NT
      REAL*8        VI(*), EIGVL(*), EIGVC(LDZ,*),DM(*) 
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
cdel      REAL*8        
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
*     Vinicius Vaz da Cruz
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     () written based on Freddy's spectrumti.f
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel     LOGICAL       
      CHARACTER*72  TITLE
      INTEGER       I, J, K
      REAL*8        EI, EN, W

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program      

c     
 1044 FORMAT(A72)
 1045 FORMAT(3X,'Energy',I3)
 1046 FORMAT(3X,E20.12)
 1047 FORMAT(3X,'Eigenvector',I3)
 1048 FORMAT(100000000(3X,E20.12)) 
 1049 FORMAT(5X,A8,8X,A6,11X)
 1053 FORMAT(I2,' ->',I2,F12.7)

c           
      WRITE(*,*)
      WRITE(*,*)'Transition Dipole Moments:'
      WRITE(*,*)'========='
      WRITE(*,1049)'*I->F*','*Mij*'   
      DO K=1,IIU
         DO I=K, IIU
            W = ZERO
            DO J=1,NT
               W = W + EIGVC(J, I)*DM(J)*EIGVC(J,K)*0.393456D+0
C 1 DEBYE = 0.393456 A.U.
c               print*, EIGVC(J, I),DM(J),EIGVC(J,K)
            ENDDO
            WRITE(*,1053)K-1,I-1, W
         ENDDO
      ENDDO
c
      RETURN
      END
