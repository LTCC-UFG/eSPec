      SUBROUTINE PSODWLT(DIM, PRTCRL, INFO, LC, MP, MPI, MPR, MXDCT, N, 
     &     DT, TI, TF, TOL, NP, SHM, SH, SM, U1, U2, V1, V2, XI, XF, 
     &     TC, VPOT, WRKR, WRKI, VAR, LNZVC) 
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
      REAL*8        SHM(*), SH(*), SM(*), U1(*), U2(*), V1(*), V2(*)
      REAL*8        XI(*), XF(*), TC(*), VPOT(*), WRKR(*), WRKI(*)
      REAL*8        VAR(*), LNZVC(MXDCT,*)
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
      INTEGER       I, MC, NC, NDIM
      REAL*8        CR, CI, CT, CST, E0, ER, EI, ET, F0, FNR, FNI, FNT
      REAL*8        SHFS, T, B
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
      EXTERNAL      PSODWLTAUX
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
      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
      CT = SQRT(CR**2 + CI**2)
      LNZVC(NC,MPR) = CR 
      LNZVC(NC,MPI) = CI
c
      ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WRKR, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WRKR, VAR)
      ET  = (ER + EI)
      E0 = ET
      FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
      FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
      FNT = SQRT(FNR**2 + FNI**2)
      F0 = FNT
c
      CALL DIGT(DIM, NDIM, N, B, NP, SH, SM, TC)
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
 1011    FORMAT(/,7X,A6,7X,A8,8X,A7,7X,A7,6X,A8,6X,A5)
 1021    FORMAT(1X,F14.8,3X,E12.6,2(4X,F10.7),4X,F10.8,4X,F8.6)
         WRITE(*,*)
         WRITE(*,*)'Partial auto-correlation function:'
         WRITE(*,*)'=================================='
         WRITE(*,1011)'*t/fs*', '*E/a.u.*', '*Cr(t)*', '*Ci(t)*',
     &        '*|C(t)|*', '*|C|*'
         WRITE(*,1021)T, ET, CR, CI, CT, FNT
      ENDIF
c

      CALL PSODWLTAUX(NDIM, N, B, CST/2D0, NP, U1, U2, V1, V2, 
     &     TC, VPOT, WRKR, WRKI)
      CALL PSODWLTAUX(NDIM, N, B, CST/2D0, NP, U1, U2, V1, V2, 
     &     TC, VPOT, WRKR, WRKI)
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
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WRKR, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WRKR, VAR)
         ET  = (ER + EI)
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
         CALL PSODWLTAUX(NDIM, N, B, CST, NP, U1, U2, V1, V2, 
     &        TC, VPOT, WRKR, WRKI) 
c     
         CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
         CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
         CT = SQRT(CR**2 + CI**2)
c
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WRKR, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WRKR, VAR)
         ET  = (ER + EI)
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
         FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
         FNT = SQRT(FNR**2 + FNI**2)
c
         IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL 
     &        .OR. ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
c     
         MC = MC + ONE
         IF(MC.EQ.LC)THEN            
            NC = NC + ONE
            MC = ZERO
            LNZVC(NC,MPR) = CR 
            LNZVC(NC,MPI) = CI
         ENDIF
c     
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
            WRITE(*,1021)T, ET, CR, CI, CT, FNT
         ENDIF
c     
      ENDDO
c     ..
      RETURN
      END
