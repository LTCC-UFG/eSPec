      SUBROUTINE PSODFFT(ABSORB, DIM, PRTCRL, TPABSOR, INFO, LC, MP, 
     &     MPI, MPR, MXDCT, N, DT, TI, TF, TOL, NP, SHM, SH, SM, U1, 
     &     U2, V1, V2, XI, XF, TC, VPOT, VABC, WRKR, WRKI, VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB
      CHARACTER*(*) DIM, PRTCRL, TPABSOR
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N
      REAL*8        DT, TI, TF, TOL
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*) 
      INTEGER       NP(*)
      REAL*8        SHM(*), SH(*), SM(*), U1(*), U2(*), V1(*), V2(*)
      REAL*8        XI(*), XF(*), TC(*), VPOT(*), VABC(*), WRKR(*)
      REAL*8        WRKI(*), VAR(*), LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
*     Wave function propagation by SOD approach:
*     |p_(j+1) = |p^(j-1) + 2*i*dt*|H* |p_(j-1)/hbar (j = 0, 1, 2, ...). 
*
*     This subroutine does the propagation of a real or an imaginary 
*     eigenvector using the SOD method.
*
*     ..
*     Arguments
*     =========
*
*
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (02/2004) First version PSODFFT written by Freddy based in PSOD subrout. 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0,
     &     TWOPI = ++6.283185307179586476925287D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER       I, MC, NC, NDIM
      REAL*8        CR, CI, CT, CST, E0, ET, F0, FNR, FNI, FNT
      REAL*8        SHFS, T, B, B1, B2!, et2, et4, EI, ER
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
      REAL*8        EFFT4, FCORREL!, EIGENERG EFFT2   
c     **
c     ** External subroutines 
      EXTERNAL      PSODFFTAUX
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT 
c     .. Start program      
c     ..
c     .. Starting values
      INFO = ZERO
      NC = ONE
      T = TI
      SHFS = 2.4188843243D-2 !fs/a.u.(t)
      CST = TWO*DT/SHFS 
      DO I=1,N,1
         U0(I) = U1(I)
         U2(I) = U0(I)
         V0(I) = V1(I)
         V2(I) = V0(I)
      ENDDO
c
      CALL DIGT(DIM, NDIM, N, B, B1, B2, NP, SH, SM, TC)
c
      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
      CT = SQRT(CR**2 + CI**2)
      LNZVC(NC,MPR) = CR 
      LNZVC(NC,MPI) = CI
D      ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WRKR)
D      EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WRKR)
D      ET  = (ER + EI)
      ET  = EFFT4(NDIM, N, NP, B, U2, V2, WRKR, WRKI, TC, VPOT)
      E0 = ET
      FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
      FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
      FNT = SQRT(FNR**2 + FNI**2)
      F0 = FNT
c
      CALL DIGT(DIM, NDIM, N, B, B1, B2, NP, SH, SM, TC)
c      DO I=1,N,1
c         TC(I) = TWOPI*TC(I)
c      ENDDO
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
      CALL PSODFFTAUX(ABSORB, TPABSOR, NDIM, N, B, CST/2D0, NP, U1, U2, 
     &     V1, V2, TC, VPOT, VABC, WRKR, WRKI)
      CALL PSODFFTAUX(ABSORB, TPABSOR, NDIM, N, B, CST/2D0, NP, U1, U2, 
     &     V1, V2, TC, VPOT, VABC, WRKR, WRKI)
c     
      MC = ONE
      T = TI + DT
      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
      CT = SQRT(CR**2 + CI**2)
      IF(MC.EQ.LC)THEN
         NC = NC + ONE
         MC = ZERO
         LNZVC(NC,MPR) = CR 
         LNZVC(NC,MPI) = CI
      ENDIF
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &     PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
D         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WRKR)
D         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WRKR)
D         ET  = (ER + EI)
         ET  = EFFT4(NDIM, N, NP, B, U2, V2, WRKR, WRKI, TC, VPOT)
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2) 
         FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
         FNT = SQRT(FNR**2 + FNI**2)
         WRITE(*,1021)T, ET, CR, CI, CT, FNT
      ENDIF
c
      DO I=1,N,1
         U1(I) = U0(I)
         V1(I) = V0(I)
      ENDDO
c
      DO T=TI+2*DT,TF+DT,DT
c     
         CALL PSODFFTAUX(ABSORB, TPABSOR, NDIM, N, B, CST, NP, U1, U2, 
     &        V1, V2, TC, VPOT, VABC, WRKR, WRKI) 
c
         MC = MC + ONE
         IF(MC.EQ.LC)THEN 
            CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
            CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
            NC = NC + ONE
            MC = ZERO
            LNZVC(NC,MPR) = CR 
            LNZVC(NC,MPI) = CI
         ENDIF
c     
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
            CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
            CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
            CT = SQRT(CR**2 + CI**2)
            ET  = EFFT4(NDIM, N, NP, B, U2, V2, WRKR, WRKI, TC, VPOT)
            FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
            FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
            FNT = SQRT(FNR**2 + FNI**2)
            IF(ABSORB)THEN
               IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL)INFO = ONE
            ELSE
               IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL 
     &              .OR. ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
            ENDIF
            WRITE(*,1021)T, ET, CR, CI, CT, FNT
ccccccccccccccccccccccccccc
D            ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WRKR)
D            EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WRKR)
D            ET4 = (ER + EI)
D            ET2 = EFFT2(NDIM, N, NP, B, B1, U2, V2, TC, VPOT)
D            write(*,'(E12.6,3X,E12.6)')et2, et4
D            write(*,*)'ET4,ET',ET4,ET,ER,EI
D            read(*,*)
cccccccccccccccccccccccccc
         ENDIF     
      ENDDO
c
      ET  = EFFT4(NDIM, N, NP, B, U2, V2, WRKR, WRKI, TC, VPOT)
      FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
      FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
      FNT = SQRT(FNR**2 + FNI**2)
      IF(ABSORB)THEN
         IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL)INFO = ONE
      ELSE
         IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL 
     &        .OR. ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
      ENDIF
c     ..
      RETURN
 1011 FORMAT(/,7X,A6,7X,A8,8X,A7,7X,A7,6X,A8,6X,A5)
 1021 FORMAT(1X,F14.8,3X,E12.6,2(4X,F10.7),4X,F10.8,4X,F8.6)
      END
