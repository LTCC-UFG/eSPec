      SUBROUTINE GVPOT(SHIFTP, POTCH, DIM, NP, XI, XF, XMIN, SM, SH, 
     &     CST, VPOT)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL           SHIFTP
      CHARACTER*(*)     POTCH, DIM  
cdel      INTEGER  
      REAL*8            XMIN
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        XI(*), XF(*), SM(*), SH(*), CST(*), VPOT(*)
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
*     (13/03/2003) First version GVPOT written by Freddy. 
*
c     **
c     ** Parameters 
cdel      REAL*8        ZERO, ONE
cdel      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, J, K, L 
      REAL*8        RX, RY, RZ, VPOTAUX
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER       
      REAL*8        POT1D, POT2D, POT3D
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program
      VPOTAUX = 1.0D+99
c
      IF(DIM(1:3).EQ.'.1D')THEN
         DO I=1,NP(1),1
            RX = XI(1) + (I - 1)*SH(1)
            VPOT(I) = POT1D(POTCH, CST, RX)
            IF(SHIFTP .AND. VPOT(I).LT.VPOTAUX)THEN
               VPOTAUX = VPOT(I)
               XMIN = RX
            ENDIF
         ENDDO
         IF(SHIFTP .AND. VPOTAUX.GT.1.0D-6)THEN
            DO I=1,NP(1),1
               VPOT(I) = VPOT(I) - VPOTAUX
            ENDDO
         ENDIF       
      ELSEIF(DIM(1:3).EQ.'.2D')THEN
         DO J=1,NP(2),1
c            write(*,*)POTCH
            RY = XI(2) + (J - 1)*SH(2)
            DO I=1,NP(1),1
               RX = XI(1) + (I - 1)*SH(1)
               L = I + (J - 1)*NP(1)
               VPOT(L) = POT2D(POTCH, SM, CST, RX, RY)
               IF(SHIFTP .AND. VPOT(L).LT.VPOTAUX)VPOTAUX = VPOT(L)
            ENDDO
         ENDDO
         IF(SHIFTP .AND. ABS(VPOTAUX).GT.1.0D-6)THEN
            DO I=1,NP(1)*NP(2),1
               VPOT(I) = VPOT(I) - VPOTAUX
            ENDDO
         ENDIF
      ELSEIF(DIM(1:3).EQ.'.3D')THEN
         DO K=1,NP(3),1
            RZ = XI(3) + (K - 1)*SH(3)
            DO J=1,NP(2),1
               RY = XI(2) + (J - 1)*SH(2)
               DO I=1,NP(1),1
                  RX = XI(1) + (K - 1)*SH(1)
                  L = I + (J - 1)*NP(1) + (K - 1)*NP(1)*NP(2)
                  VPOT(L) = POT3D(POTCH, CST, RX, RY, RZ)
                  IF(SHIFTP .AND. VPOT(L).LT.VPOTAUX)VPOTAUX = VPOT(L)
               ENDDO
            ENDDO
         ENDDO 
         IF(SHIFTP .AND. VPOTAUX.GT.1.0D-6)THEN
            DO I=1,NP(1)*NP(2)*NP(3),1
               VPOT(I) = VPOT(I) - VPOTAUX
            ENDDO
         ENDIF
      ENDIF
c
      RETURN
      END


