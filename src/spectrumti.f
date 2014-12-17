      SUBROUTINE SPECTRUMTI(IIL, IIU, IFL, IFU, LDZ, NT, VI, EIGVL, 
     &     EIGVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       IIL, IIU, IFL, IFU, LDZ, NT
      REAL*8        VI(*), EIGVL(*), EIGVC(LDZ,*) 
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
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     () First version written by 
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
 1049 FORMAT(5X,A8,8X,A6,11X,A4,14X,A6)
 1053 FORMAT(1X,3(F12.7,4X),6X,I2,' ->',I2)
c           
      WRITE(*,*)
      WRITE(*,*)'Spectrum:'
      WRITE(*,*)'========='
      WRITE(*,1049)'*E/a.u.*','*AMPT*','*FC*', '*I->F*'   
      OPEN(1, STATUS='UNKNOWN', FILE='initial_spc.aux')
      READ(1,1044)TITLE
      READ(1,1044)TITLE
      READ(1,1044)TITLE
      DO K=1,IIU-IIL+1,1
c     .. get initials eigenvectors
         READ(1,1044)TITLE
         READ(1,1046)EI
         READ(1,1044)TITLE
         READ(1,1048)(VI(I), I=1,NT,1)
c     .. calculate intensities
         DO J=1,IFU-IFL+1,1 
            EN = EIGVL(J) - EI
            W = ZERO
            DO I=1,NT,1
               W = W + EIGVC(I, J)*VI(I)
            ENDDO
            WRITE(*,1053)EN, W, W*W, IIL+K-2, IFL+J-2
         ENDDO
      ENDDO
c
      RETURN
      END
