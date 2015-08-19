      SUBROUTINE REORT(N, M, LMTREORT, MXDCT, WORK, EIGVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*1      
      INTEGER       N, M, LMTREORT, MXDCT 
cdel      REAL*8      
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)       
cdel      INTEGER  
      REAL*8       WORK(*), EIGVC(MXDCT,*)
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
      REAL*8        ZERO
      PARAMETER     (ZERO = +0.0E+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J
cdel      REAL*8        
c     **
c     ** local arrays 
      REAL*8        RV(LMTREORT) 
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
c     .. init and test parameters 
      DO J=1,M,1 
         RV(J) = ZERO
         DO I=1,N,1
            RV(J) = RV(J) + WORK(I)*EIGVC(I,J)
         ENDDO
c         write(*,*)'[DEBUG: reort] rv(j)',j, rv(j)
      ENDDO
c
      DO J=1,M,1
         DO I=1,N,1
            WORK(I) = WORK(I) - RV(J)*EIGVC(I,J)
         ENDDO
      ENDDO 
c
      RETURN
      END
