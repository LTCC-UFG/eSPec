      SUBROUTINE SPECTRUM(TPWIND, M, DT, TF, WP, U, V) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) TPWIND
      INTEGER       M
      REAL*8        DT, TF, WP
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)      
cdel      INTEGER   
      REAL*8        U(*), V(*)       
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
*     (21/03/2003) First version SPECTRUM written by Freddy. 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO 
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, ITYPE, N
      REAL*8        CSTL, EN, ENC, T, TAUX, TMAX, UAUX, VAUX

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
      REAL*8        FCORREL
c     **
c     ** External subroutines 
      EXTERNAL      FFTNT
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT
c     .. Start program      
      ITYPE = ONE
      N = TWO**M
      ENC = ONE/(N*DT)
      EN = ZERO
c
      IF(TPWIND(1:3).EQ.'.SG')THEN
         TMAX = N*DT
         TAUX = TMAX**2/(LOG(ONE/WP)*2.302585093)
         DO I=1,N,1
            T = -TF + (I - ONE)*DT
            CSTL = DEXP(-T**2/TAUX)
            U(I) = CSTL*U(I)
            V(I) = CSTL*V(I)
         ENDDO
      ELSEIF(TPWIND(1:5).EQ.'.NONE')THEN
         CONTINUE
      ELSE
         WRITE(*,1111)
         STOP
      ENDIF
c
      DO I=1,N/2,1
         UAUX = U(I)
         U(I) = U(N/2+I)
         U(N/2+I) = UAUX
         VAUX = V(I)
         V(I) = V(N/2+I)
         V(N/2+I) = VAUX
      ENDDO
c
      CALL FFTNT(U, V, N, M, ITYPE)
c
c      DO I=1,N/2,1
c         VAUX = U(I)
c         IF(VAUX.GT.UAUX)UAUX = VAUX
c         UAUX = U(I)
c         VAUX = V(I)
c         U(I) = DT*U(I+N/2)
c         V(I) = DT*V(I+N/2)
c         U(I+N/2) = DT*UAUX
c         V(I+N/2) = DT*VAUX
c      ENDDO
      UAUX = SQRT(FCORREL(N, U, U))
c
 1011 FORMAT(6X,A7,10X,A3,11X,A6,10X,A6)
 1021 FORMAT(1X,4(F12.7,4X))
      WRITE(*,*)
      WRITE(*,*)'Spectrum:'
      WRITE(*,*)'========='
      WRITE(*,1011)'*E/fHz*','*W*','*FFTr*','*FFTi*'
      DO I=1,N/2,1 
         WRITE(*,1021)EN,U(I)/UAUX, U(I), V(I)
         EN = EN + ENC
      ENDDO
c
 1111 FORMAT('<<<>>> Windowing function error <<<>>>')
c
      RETURN
      END
