      SUBROUTINE PSPOFFT(DIM, PRTCRL, INFO, LC, MP, MPI, MPR, MXDCT, N, 
     &     DT, TI, TF, TOL, NP, SHM, SH, SM, U1, V1, XI, XF, 
     &     TC, VPOT, LNZVC, WRK1, WRK2, VAR) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) DIM, PRTCRL
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N
      REAL*8        DT, TI, TF, TOL
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        SHM(*), SH(*), SM(*), U1(*), V1(*), VAR(*)
      REAL*8        XI(*), XF(*), TC(*), VPOT(*), WRK1(*), WRK2(*)
      REAL*8        LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     |p^(2) = exp(-i*dt*V/(2*hbar))*exp(-i*dt*T/hbar)exp(-i*dt*V/(2*hbar)).
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
*     (02/2004) First version PSPO written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0,
     &     TWOPI = +6.283185307D0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, MC, NC, NDIM
      REAL*8        CR, CI, CT, CST1, CST2, E0, ER, EI, ET, F0, FNR 
      REAL*8        FNI, FNT, SHFS, T, B, B1, B2
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER       I
      REAL*8        U0(N), V0(N)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
      REAL*8        FCORREL, EIGENERG, EFFT4  
c     **
c     ** External subroutines 
      EXTERNAL      PSPOFFTAUX
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT 
c     .. Start program      
c     ..
c     .. Starting values
      INFO = ZERO
      MC = ZERO
      NC = ONE
      T = TI
      SHFS = 2.4188843243D-2 !fs/a.u.(t)
      CST1 = DT/SHFS 
      CST2 = DT/(TWO*SHFS) 
      DO I=1,N,1
         U0(I) = U1(I)
         V0(I) = V1(I)
      ENDDO
c
      CR = FCORREL(N, U0, U1) + FCORREL(N, V0, V1)
      CI = FCORREL(N, U0, V1) - FCORREL(N, V0, U1)
      CT = SQRT(CR**2 + CI**2)
      LNZVC(NC,MPR) = CR 
      LNZVC(NC,MPI) = CI
      ER = EIGENERG(DIM, N, NP, SHM, U1, VPOT, WRK1, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V1, VPOT, WRK1, VAR)
      ET  = (ER + EI)
      E0 = ET
      FNR = FCORREL(N, U1, U1) + FCORREL(N, V1, V1)
      FNI = FCORREL(N, U1, V1) - FCORREL(N, V1, U1)
      FNT = SQRT(FNR**2 + FNI**2)
      F0 = FNT
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         WRITE(*,*)
         WRITE(*,*)'Partial auto-correlation function:'
         WRITE(*,*)'=================================='
         WRITE(*,1011)'*t/fs*', '*E/a.u.*', '*Cr(t)*', '*Ci(t)*',
     &        '*|C(t)|*', '*|C|*'
         WRITE(*,1021)T, ET, CR, CI, CT, FNT
      ENDIF
c
      CALL DIGT(DIM, NDIM, N, B, B1, B2, NP, SH, SM, TC)
c
      DO T=TI+DT,TF+DT,DT
c     
         CALL PSPOFFTAUX(NDIM, N, B, B1, B2, CST1, CST2, NP, U1, V1, 
     &        TC, VPOT)
c
         MC = MC + ONE
         IF(MC.EQ.LC)THEN            
            NC = NC + ONE
            MC = ZERO
            CR = FCORREL(N, U0, U1) + FCORREL(N, V0, V1)
            CI = FCORREL(N, U0, V1) - FCORREL(N, V0, U1)
            LNZVC(NC,MPR) = CR 
            LNZVC(NC,MPI) = CI
         ENDIF
c     
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
            CR = FCORREL(N, U0, U1) + FCORREL(N, V0, V1)
            CI = FCORREL(N, U0, V1) - FCORREL(N, V0, U1)
            CT = SQRT(CR**2 + CI**2)
            ET  = EFFT4(NDIM, N, NP, B, U1, V1, WRK1, WRK2, TC, VPOT)
            FNR = FCORREL(N, U1, U1) + FCORREL(N, V1, V1)
            FNI = FCORREL(N, U1, V1) - FCORREL(N, V1, U1)
            FNT = SQRT(FNR**2 + FNI**2)
            IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL .OR. 
     &           ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
            WRITE(*,1021)T, ET, CR, CI, CT, FNT
         ENDIF     
      ENDDO
c     ..
 1011 FORMAT(/,7X,A6,7X,A8,8X,A7,7X,A7,6X,A8,6X,A5)
 1021 FORMAT(1X,F14.8,3X,E12.6,2(4X,F10.7),4X,F10.8,4X,F8.6)
      RETURN
      END
