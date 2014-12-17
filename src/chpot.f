      SUBROUTINE CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, DM, VP, 
     &     VEFF)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL      
      CHARACTER*(*) EFC
      INTEGER       N, KL
      REAL*8        E0, T0, TD, TP, T, OMG, SNI
c     **
c     ** Array arguments
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        DM(*), VP(*), VEFF(*)
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
*     (24/06/2003) First version CHPOT written by Freddy. 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL       
cdel      CHARACTER*1   
      INTEGER       I
      REAL*8        EFAUX
c     **
c     ** External functions 
cdel      LOGICAL       
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        EF
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
      EFAUX = EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
c      if(t.eq.t0)write(*,*)'DM(I)*EFAUX'
      DO I=1,N,1
         VEFF(I) = VP(I) - DM(I)*EFAUX
      ENDDO
c
      RETURN
      END

