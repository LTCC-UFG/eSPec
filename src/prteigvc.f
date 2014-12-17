      SUBROUTINE PRTEIGVC(NAME, MC, MXDCT, ND, NIS, NP, XP, XI, SH, EVT)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) NAME
      INTEGER       MC, MXDCT, ND, NIS
cdel      REAL*8        
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        XP(*) , XI(*) , SH(*), EVT(MXDCT,*) 
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
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I, J, K, L, NC
cdel      REAL*8        

c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       ICHLENGTH
cdel      REAL*8        
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      NC = ICHLENGTH(NAME, 0)
      OPEN(UNIT=9,STATUS='UNKNOWN',FILE=NAME(1:NC))
c      IF(NIS.GE.10)THEN
c         WRITE(MC,1010)('*R',I,'*', I=ND,1,-1),('*C',I-1,'*', 
c     &        I=1,NIS,1)
c      ELSE
         WRITE(MC,1011)('*R',I,'*', I=ND,1,-1),('*C',I-1,'*', 
     &        I=1,NIS,1)
c      ENDIF
c
      DO L=1,NP(3),1
         XP(3) = XI(3) + (L - 1)*SH(3)
         DO K=1,NP(2),1
            XP(2) = XI(2) + (K - 1)*SH(2)
            DO J=1,NP(1),1
               XP(1) = XI(1) + (J - 1)*SH(1)
               NC = J + (K - ONE)*NP(1) + (L - ONE)*NP(1)*NP(2)
               WRITE(MC,1012)(XP(I), I=ND,1,-1),
     &              (EVT(NC, I), I=1,NIS,1)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(9)
c
 1010 FORMAT(5X,100(A2,I2,A1,11X)) 
 1011 FORMAT(5X,100(A2,I1,A1,11X)) 
 1012 FORMAT(1X,100(E12.6,3X))     
c
      RETURN
      END
