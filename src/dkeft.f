      SUBROUTINE DIGT(DIM, NDIM, N, B, B1, B2, NP, SH, SM, TC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL       
      CHARACTER*(*) DIM
      INTEGER       NDIM, N
      REAL*8        B, B1, B2
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        SH(*), SM(*), TC(*)
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
      REAL*8        ZERO, ONE, TWO, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0,
     &     TWOPI = +6.283185307D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I, J, K
      REAL*8        P1, P2, P3
c     **
c     ** Local vectors
cdel      LOGICAL       
cdel      CHARACTER*1   
cdel      INTEGER       
      REAL*8        DP(3)
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
      IF(DIM(1:3).EQ.'.1D')THEN
         NDIM = 1
         B = ONE/NP(1)
         B1 = SH(1)
         B2 = ONE/(NP(1)*SH(1))
         DP(1) = TWOPI/(NP(1)*SH(1))
         P1 = ZERO
         DO I=1,N/2,1
            TC(I) = P1**2/(TWO*SM(1))
            P1 = P1 + DP(1)
         ENDDO
         P1 = -P1 
         DO I=N/2+1,N,1
            TC(I) = P1**2/(TWO*SM(1))
            P1 = P1 + DP(1)
         ENDDO
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
         NDIM = 2
         B = ONE/(NP(1)*NP(2))
         DP(1) = TWOPI/(NP(1)*SH(1)) 
         DP(2) = TWOPI/(NP(2)*SH(2)) 
         P1 = ZERO
         P2 = ZERO
         DO J=1,NP(2)/2,1
            DO I=1,NP(1)/2,1
               TC(I) = P1**2/(TWO*SM(1)) + P2**2/(TWO*SM(2))
               P1 = P1 + DP(1)
            ENDDO
            P2 = P2 + DP(2)
         ENDDO
         P1 = -(P1 + DP(1))
         P2 = -(P2 + DP(2))
         DO J=NP(2)/2+1,NP(2),1
            DO I=NP(1)/2+1,NP(1),1
               TC(I) = P1**2/(TWO*SM(1)) + P2**2/(TWO*SM(2))
               P1 = P1 + DP(1)
            ENDDO
            P2 = P2 + DP(2)
         ENDDO

      ELSEIF(DIM(1:3).EQ.'.3D')THEN
         NDIM = 3
         B = ONE/(NP(1)*NP(2)*NP(3))
         DP(1) = TWOPI/(NP(1)*SH(1)) 
         DP(2) = TWOPI/(NP(2)*SH(2)) 
         DP(3) = TWOPI/(NP(3)*SH(3)) 
         P1 = ZERO
         P2 = ZERO
         P3 = ZERO
         DO K=1,NP(3)/2,1
            DO J=1,NP(2)/2,1
               DO I=1,NP(1)/2,1
                  TC(I) = P1**2/(TWO*SM(1)) + P2**2/(TWO*SM(2)) 
     &                 + P3**2/(TWO*SM(3))
                  P1 = P1 + DP(1)
               ENDDO
               P2 = P2 + DP(2)
            ENDDO
            P3 = P3 + DP(3)
         ENDDO                
         P1 = -(P1 + DP(1))
         P2 = -(P2 + DP(2))
         P3 = -(P3 + DP(3))
         DO K=NP(3)/2+1,NP(3),1
            DO J=NP(2)/2,NP(2),1
               DO I=NP(1)/2,NP(1),1
                  TC(I) = P1**2/(TWO*SM(1)) + P2**2/(TWO*SM(2)) 
     &                 + P3**2/(TWO*SM(3))
                  P1 = P1 + DP(1)
               ENDDO
               P2 = P2 + DP(2)
            ENDDO
            P3 = P3 + DP(3)
         ENDDO       
      ELSE
         WRITE(*,*)'<<<>>> Dimension error <<<>>>' 
         STOP
      ENDIF
c
      RETURN
      END
