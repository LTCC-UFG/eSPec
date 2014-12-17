      SUBROUTINE P2ECSOD(ABSORB, DIM, PRTCRL, TPABSOR, INFO, LC, MP, MPI, 
     &     MPR, MXDCT, N, DT, TI, TF, TOL, NP, SHM, U1, U2, V1, V2, XI, 
     &     XF, VPOT, VABC, WORK, VAR, LNZVC) 
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
      REAL*8        SHM(*), XI(*), XF(*), VPOT(*), VABC(*)
      REAL*8        U1(*), U2(*), V1(*), V2(*), WORK(*), VAR(*)
      REAL*8        LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
*     Wave function propagation by SOD aproach:
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
*     (21/03/2003) First version PSOD written by Freddy.  
*     (29/07/2003) Changed in order to include propagation of an initial 
*     imaginary eigenvector.  
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
      INTEGER       I, MC, NC
      REAL*8        CR, CI, CT, CST, E0, ER, EI, ET, F0, FNR, FNI, FNT
      REAL*8        SHFS, T
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
      REAL*8        FCORREL, EIGENERG    
c     **
c     ** External subroutines 
      EXTERNAL      PSODAUX
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
         U2(I) = U1(I)
         V0(I) = V1(I)
         V2(I) = V1(I)
      ENDDO
c
      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
      CT = SQRT(CR**2 + CI**2)
      LNZVC(NC,MPR) = CR 
      LNZVC(NC,MPI) = CI
      ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
      ET  = (ER + EI)
      E0 = ET
      FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
      FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
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
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/2, NP, SHM, U1, U2, V1, 
     &        V2, VPOT, VABC, WORK, VAR)
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/2, NP, SHM, U1, U2, V1, 
     &        V2, VPOT, VABC, WORK, VAR)
c
      MC = ONE
      T = TI + DT
      IF(MC.EQ.LC)THEN
         NC = NC + ONE
         MC = ZERO  
         CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
         CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
         LNZVC(NC,MPR) = CR 
         LNZVC(NC,MPI) = CI
      ENDIF
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &     PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
         CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
         CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
         CT = SQRT(CR**2 + CI**2)
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
         ET  = (ER + EI)
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2) 
         FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
         FNT = SQRT(FNR**2 + FNI**2)
         WRITE(*,1021)T, ET, CR, CI, CT, FNT
      ENDIF
c
      DO I=1,N,2
         U1(I) = U0(I)
         U1(I+1) = U0(I+1)
         V1(I) = V0(I)
         V1(I+1) = V0(I)
      ENDDO
c
      DO T=TI+2*DT,TF+DT,DT
c     
         CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST, NP, SHM, U1, U2, V1, 
     &        V2, VPOT, VABC, WORK, VAR) 
c     
         MC = MC + ONE
         IF(MC.EQ.LC)THEN
            NC = NC + ONE
            MC = ZERO
            CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
            CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
            LNZVC(NC,MPR) = CR 
            LNZVC(NC,MPI) = CI
         ENDIF
c     
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
            CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
            CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
            CT = SQRT(CR**2 + CI**2)
            ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
            EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
            ET  = (ER + EI)
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
         ENDIF
      ENDDO
c
      ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
      ET  = (ER + EI)
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
