c     subroutine PSPOFFTAUX
c     =====================
      SUBROUTINE PSPOFFTAUX(NDIM, N, B, B1, B2, CST1, CST2, NP,
     &     U, V, TC, VPOT)
      IMPLICIT NONE
c     **
c     ** Scalar arguments
cdel      LOGICAL
cdel      CHARACTER*(*) DIM
      INTEGER       NDIM, N
      REAL*8        B, B1, B2, CST1, CST2
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*)
      REAL*8        U(*), V(*), TC(*), VPOT(*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     |p^(2) = exp(-i*dt/hbar*V/2)*exp(-i*dt/hbar*T)*exp(-i*dt/hbar*V/2)|p^(1).
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
*     (02/2004) First version PSPOAUX written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, NF, NH
      REAL*8        UAUX, VAUX, TH1, TH2, CTH1, CTH2, STH1, STH2 

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL      FFT
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      do i=1,n,1
c         write(11,*)i, tc(I)!U(I), V(I)
c      enddo
c      write(11,*)
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
      NH = N/TWO
      DO I=1,NH,1
         TH1 = CST2*VPOT(I)
         CTH1 = DCOS(TH1)
         STH1 = DSIN(TH1)
         UAUX = CTH1*U(I) + STH1*V(I)
         VAUX = CTH1*V(I) - STH1*U(I)
         TH2 = CST2*VPOT(I+NH)
         CTH2 = DCOS(TH2)
         STH2 = DSIN(TH2)
         U(I) = CTH2*U(I+NH) + STH2*V(I+NH)
         V(I) = CTH2*V(I+NH) - STH2*U(I+NH)
         U(I+NH) = UAUX
         V(I+NH) = VAUX
      ENDDO
c
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      do i=1,n,1
c         write(12,*)i,U(I), V(I)
c      enddo
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      write(*,*)ndim,n,np(1),np(2),np(3)
c      write(*,*)DCOS(CST2*VPOT(I)), DSIN(CST2*VPOT(I+NH)),
c     &     DCOS(CST1*TC(I)), DSIN(CST1*TC(I))
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(U, V, N, NP(I), NF, 1)
      ENDDO
c
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      do i=1,n,1
c         write(13,*)i,U(I), V(I)
c      enddo
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
      DO I=1,N,1
         TH1 = CST1*TC(I)
         CTH1 = DCOS(TH1)
         STH1 = DSIN(TH1)
         U(I) = CTH1*U(I) + STH1*V(I)
         V(I) = CTH1*V(I) - STH1*U(I)
      ENDDO
c
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      do i=1,n,1
c         write(14,*)i,U(I), V(I)
c      enddo
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
      NF = ONE
      DO I=1,NDIM,1
         NF = NF*NP(I)
         CALL FFT(U, V, N, NP(I), NF, -1)
      ENDDO
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      do i=1,n,1
c         write(15,*)i,U(I), V(I)
c      enddo
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c
      DO I=1,NH,1
         TH1 = CST2*VPOT(I+NH)
         CTH1 = DCOS(TH1)
         STH1 = DSIN(TH1)
         UAUX = CTH1*B*U(I) + STH1*B*V(I)
         VAUX = CTH1*B*V(I) - STH1*B*U(I)
         TH2 = CST2*VPOT(I)
         CTH2 = DCOS(TH2)
         STH2 = DSIN(TH2)
         U(I) = CTH2*B*U(I+NH) + STH2*B*V(I+NH)
         V(I) = CTH2*B*V(I+NH) - STH2*B*U(I+NH)
         U(I+NH) = UAUX
         V(I+NH) = VAUX
      ENDDO
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c      do i=1,n,1
c         write(16,*)i,U(I), V(I)
c      enddo
cccccccccccccccccccccccccccc
cccccccccccccccccccccccccccc
c     ..
      RETURN
      END
