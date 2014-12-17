      SUBROUTINE RIF(DR, N, ND, NP, CS, XI, XF, SH, PL, WORK, F) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
c      LOGICAL
c      CHARACTER      
      INTEGER        N, ND
      REAL*8         CS     
c     **
c     ** Array arguments
c      LOGICAL
      CHARACTER*1    DR      
      INTEGER        NP(*)  
      REAL*8         XI(*), XF(*), SH(*), PL(*), WORK(*), F(*)       
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
*     (19/03/2003) First version RIF written by Freddy [1D, 2D, 3D] 
*     
**
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
      REAL*8          AUX1
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
c      INTRINSIC     
c     .. Start program
c
c      write(*,*)'oi'
      IF(ND.EQ.1)THEN
         DO I=1,NP(1),1
            WORK(I) = (XI(1) + (I - 1)*SH(1))*F(I)
         ENDDO
      ELSEIF(ND.EQ.2)THEN
         IF(DR.EQ.'i')THEN
            II = ZERO
            DO J=1,NP(2),1
               DO I=1,NP(1),1
                  II = II + ONE
                  WORK(II) = (XI(1) + (I - 1)*SH(1))*F(II)
               ENDDO
            ENDDO
         ELSEIF(DR.EQ.'j')THEN
            DO J=1,NP(1),1
               AUX1 = (XI(2) + (J - 1)*SH(2))
               DO I=1,NP(2),1
                  II = J + (I - 1)*NP(1)
                  WORK(II) = AUX1*F(II)/(TWO*SH(2))
               ENDDO
            ENDDO
         ELSE
            WRITE(*,1011)
            STOP
         ENDIF
ccccccccccccccccccc
ccccccccccccccccccccccccccccccc
      ELSEIF(ND.EQ.3)THEN
         IF(DR.EQ.'i')THEN
            II = ZERO
            DO K=1,NP(3),1
               DO J=1,NP(2),1
                  DO I=1,NP(1),1
                     II = II + ONE
                     WORK(II) = (XI(1) + (I - 1)*SH(1))*F(II)
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF(DR.EQ.'j')THEN
            DO K=1,NP(1),1 
               DO J=1,NP(3),1
                  DO I=1,NP(2),1
                     II = I + (J - 1)*NP(1) 
                     WORK(II) = (XI(2) + (I - 1)*SH(2))*F(II)
     &                    /(TWO*SH(2))
                  ENDDO
               ENDDO  
            ENDDO
         ELSEIF(DR.EQ.'k')THEN
            DO K=1,NP(2),1 
               DO J=1,NP(1),1
                  DO I=1,NP(3),1
                     II = K + (J - 1)*NP(1) + (I - 1)*NP(1)*NP(2) 
                     WORK(II) = (XI(3) + (I - 1)*SH(3))*F(II)
     &                    /(TWO*SH(2))
                  ENDDO
               ENDDO  
            ENDDO
            
         ELSE
            WRITE(*,1011)
            STOP
         ENDIF   
c
      ELSE
         WRITE(*,*)'<<<>>> Dimesion error <<<>>>'
         STOP
      ENDIF
c
      CS = ONE/SQRT(FCORREL(N, WORK, WORK))
c
      DO I=1,N,1
         F(I) = CS*WORK(I)
      ENDDO
c 
 1011 FORMAT('<<<>>> Direction error <<<>>>')
c     
      RETURN
      END
