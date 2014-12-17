      SUBROUTINE DFDXI(DR, N, ND, NP, CS, SH, PL, WORK, F) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
c      LOGICAL
      CHARACTER*1    DR      
      INTEGER        N, ND
      REAL*8         CS     
c     **
c     ** Array arguments
c      LOGICAL
      INTEGER        NP(*)  
      REAL*8         SH(*), PL(*), WORK(*), F(*)       
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
*     (19/03/2003) First version DFDXI written by Freddy [1D, 2D] 
*     
*
*
c     **
c     ** Parameters 
      REAL*8         ZERO, ONE, TWO
      PARAMETER      (ZERO = +0.0E+0, ONE = +1.0E+0, TWO = +2.0E+0)
c     **
c     ** Local scalars 
c      LOGICAL
c      CHARACTER
      INTEGER         I, II, J, K
c      REAL*8          DELX
c     **
c     ** External functions 
c      LOGICAL
c      CHARACTER
c      INTEGER   
      REAL*8          FCORREL
c     **
c     ** External subroutines 
c      EXTERNAL      
c     **
c     ** Intrinsic functions 
      INTRINSIC       SQRT     
c     .. Start program
      IF(ND.EQ.1)THEN
         WORK(1) = ZERO
         DO I=2,NP(1)-1,1
            WORK(I) = PL(1)*(F(I-1) - F(I+1))/(TWO*SH(1))
         ENDDO
         WORK(NP(1)) = ZERO
      ELSEIF(ND.EQ.2)THEN
         IF(DR.EQ.'i')THEN
            II = ZERO
            DO J=1,NP(2),1
               II = II + ONE
               WORK(II) = ZERO
               DO I=2,NP(1)-1,1
                  II = II + ONE
                  WORK(II) = PL(1)*(F(II-1) - F(II+1))/(TWO*SH(1))
               ENDDO
               II = II +1
               WORK(II) = ZERO
            ENDDO
         ELSEIF(DR.EQ.'j')THEN
            DO J=1,NP(1),1
               II = J 
               WORK(II) = ZERO
               DO I=2,NP(2)-1,1
                  II = J + (I - 1)*NP(1)
                  WORK(II) = PL(1)*(F(II-NP(1)) - F(II+NP(1)))
     &                 /(TWO*SH(2))
               ENDDO
               II = J + (NP(2) - 1)*NP(1)
               WORK(II) = ZERO
            ENDDO
         ELSE
            WRITE(*,1011)
            STOP
         ENDIF
c
      ELSEIF(ND.EQ.3)THEN
         IF(DR.EQ.'i')THEN
            II = ZERO
            DO K=1,NP(3),1
               II = II + ONE
               WORK(II) = ZERO
               DO J=2,NP(2)-1,1
                  II = II + ONE
                  WORK(II) = ZERO
                  DO I=2,NP(1)-1,1
                     II = II + ONE
                     WORK(II) = PL(1)*(F(II-1) - F(II+1))/(TWO*SH(1))
                  ENDDO
                  II = II +1
                  WORK(II) = ZERO
               ENDDO
               II = II +1
               WORK(II) = ZERO
            ENDDO
         ELSEIF(DR.EQ.'j')THEN
            II = ZERO
            DO K=1,NP(1),1
               II = K 
               WORK(II) = ZERO
               DO J=2,NP(3)-1,1
                  II = K + (J - 1)*NP(2)
                  WORK(II) = ZERO
                  DO I=2,NP(2)-1,1
                     II = K + (J - 1)*NP(2) + (I - 1)*NP(1)*NP(2)
                     WORK(II) = PL(1)*(F(II-NP(1)) - F(II+NP(1)))
     &                    /(TWO*SH(1))
                  ENDDO
                  II = K + (J - 1)*NP(2) + (NP(2) - 1)*NP(1)*NP(2)
                  WORK(II) = ZERO
               ENDDO
               II = K + (J - 1)*NP(2) + (NP(2) - 1)*NP(1)*NP(2)
               WORK(II) = ZERO
            ENDDO

         ELSEIF(DR.EQ.'k')THEN
         
         ELSE
            WRITE(*,1011)
            STOP
         ENDIF   
c
      ELSE
         WRITE(*,*)
         STOP
      ENDIF
c
      CS = ONE/SQRT(FCORREL(N, WORK, WORK))
      DO I=1,N,1
         F(I) = CS*WORK(I)
      ENDDO
 1011 FORMAT('<<<>>> Direction error <<<>>>')
c     
      RETURN
      END
