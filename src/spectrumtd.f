      SUBROUTINE SPECTRUMTD(TPWIND, M, DT, TF, WP, DE, U, V) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) TPWIND
      INTEGER       M
      REAL*8        DT, TF, WP, DE
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
      REAL*8        ZERO, ONE, TWO, FATE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0,
     &     FATE = 27.2113834D0)
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
      ITYPE = - ONE
      N = TWO**M
      ENC = ONE/(N*DT)*1.519829846003D-1
      ENC = ENC*FATE
      EN = -N/TWO*ENC
      IF(DE.LT.ZERO)DE = ZERO
c
      IF(TPWIND(1:3).EQ.'.SG')THEN
         TMAX = N*DT
         TAUX = TMAX**2/(LOG(ONE/WP))**2
c         write(1,*)'#', wp, taux
         DO I=1,N,1
            T = -TF + (I - ONE)*DT
            CSTL = DEXP(-T**2/TAUX)
c            write(1,*)t, CSTL
            U(I) = CSTL*U(I)
            V(I) = CSTL*V(I)
c            write(1,*)t, u(i), v(i)
         ENDDO
      ELSEIF(TPWIND(1:5).EQ.'.NONE' .OR. 
     &        TPWIND(1:5).EQ.'.NULL')THEN
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
      UAUX = SQRT(FCORREL(N, U, U))
c
      WRITE(*,*)
      WRITE(*,*)'Spectrum:'
      WRITE(*,*)'========='
      WRITE(*,1011)
      DO I=N/2+1,N,1
         WRITE(*,1021)EN + DE, U(I), U(I), V(I)
         EN = EN + ENC
      ENDDO
      DO I=1,N/2,1
         WRITE(*,1021)EN + DE, U(I), U(I), V(I)
         EN = EN + ENC
      ENDDO
c
 1011 FORMAT(5X,'*E (eV)*',11X,'*W*',11X,'*FFTr*',10X,'*FFTi*')
 1021 FORMAT(1X,4(F12.7,4X))
 1111 FORMAT('<<<>>> Windowing function error <<<>>>')
c
      RETURN
      END
