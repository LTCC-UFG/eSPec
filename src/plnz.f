      SUBROUTINE PLNZ(ABSORB, CHANGE, DIM, PRTCRL, TPABSOR, INFO, LC, 
     &     LMTREORT, MP, MPI, MPR, MXDCT, N, ABSTOL, DT, TI, TF, TOL, 
     &     NP, IWORK, SHM, VPOT, VABC, U, V, XI, XF, WU, WV, WP, WORK, 
     &     VAR, LNZVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, CHANGE
      CHARACTER*(*) DIM, PRTCRL, TPABSOR          
      INTEGER       INFO, LC, LMTREORT, MP, MPI, MPR, MXDCT, N
      REAL*8        ABSTOL, DT, TI, TF, TOL
c     **
c     ** Array arguments
cdel      LOGICAL
      INTEGER       NP(*), IWORK(*)
      REAL*8        SHM(*), VPOT(*), VABC(*), U(*), V(*), XI(*), XF(*) 
      REAL*8        WU(*), WV(*), WP(*), WORK(*), VAR(*), LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
*     Wave function propagation by LANCZOS aproach:
*     |p_(j+1) =  |L'*|Q'*e^(2*i*|D*dt/hbar)*|Q*|L*|p_(j). 
*     ..
*     Arguments
*     =========
*     MP number of lanczos vector to be used.
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (19/03/2003) First version LNZ written by Freddy. 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, 
     &     TWOPI = +6.283185307179586476925287D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER        I, MC, NC
      REAL*8         CR, CI, CT, E0, ER, EI, ET, F0, FNR, FNI, FNT
      REAL*8         SHFS, CST, T
c     **
c     ** Local arrays
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
      REAL*8         V0(N), U0(N)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
      REAL*8         FCORREL, EIGENERG!, ECNORM
c     **
c     ** External subroutines 
      EXTERNAL       PLNZAUX  
c     **
c     ** Intrinsic functions 
      INTRINSIC      SQRT     
c     .. Start program
c     ..
c     .. Initiate valors
      NC = ONE
      SHFS = 2.4188843243D-2 !fs/a.u.(t)
      CST = DT/SHFS
      T = TI
      DO I=1,N,1
         U0(I) = U(I)
         V0(I) = V(I)
         WU(I) = ZERO
         WV(I) = ZERO
         WP(I) = ZERO
         WORK(I) = ZERO
      ENDDO
c
      CR = FCORREL(N, U0, U) + FCORREL(N, V0, V)
      CI = FCORREL(N, U0, V) - FCORREL(N, V0, U)
      CT = CR**2 + CI**2
      LNZVC(NC,MPR) = CR 
      LNZVC(NC,MPI) = CI
      ER = EIGENERG(DIM, N, NP, SHM, U, VPOT, WORK, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V, VPOT, WORK, VAR)
      ET  = (ER + EI)
      E0 = ET
      FNR = FCORREL(N, U, U) + FCORREL(N, V, V)
c     FNI must be igual to zero
c     FNI = FCORREL(N, U, V) - FCORREL(N, V, U)
      FNT = FNR
      F0 = FNT
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         WRITE(*,*)
         WRITE(*,*)'Partial auto-correlation function:'
         WRITE(*,*)'=================================='
         WRITE(*,1011)'*t*', '*Energy*', '*Cr(t)*', '*Ci(t)*',
     &        '*|C(t)|^2*','*|C|*'
         WRITE(*,1021)T, ET, CR, CI, CT, FNT 
      ENDIF
c
      MC = ZERO
      DO T=TI+DT,TF+DT,DT
c         
         CALL PLNZAUX(ABSORB, CHANGE, DIM, TPABSOR, INFO, LMTREORT, 
     &        MP, MXDCT, N, ABSTOL, CST, NP, IWORK, SHM, VPOT, VABC, 
     &        U, V, WU, WV, WP, WORK, VAR, LNZVC)
c
         MC = MC + ONE
         IF(MC.EQ.LC)THEN
            NC = NC + ONE
            MC = ZERO
            CR = FCORREL(N, U0, U) + FCORREL(N, V0, V)
            CI = FCORREL(N, U0, V) - FCORREL(N, V0, U)
            LNZVC(NC,MPR) = CR 
            LNZVC(NC,MPI) = CI
         ENDIF
c
         IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL' 
     &        .AND. MC.EQ.ZERO)THEN              
            CR = FCORREL(N, U0, U) + FCORREL(N, V0, V)
            CI = FCORREL(N, U0, V) - FCORREL(N, V0, U)
            CT = CR**2 + CI**2
            ER = EIGENERG(DIM, N, NP, SHM, U, VPOT, WORK, VAR)
            EI = EIGENERG(DIM, N, NP, SHM, V, VPOT, WORK, VAR)
            ET  = (ER + EI)
            FNR = FCORREL(N, U, U) + FCORREL(N, V, V)
c     FNI must be igual to zero
c            FNI = FCORREL(N, U, V) - FCORREL(N, V, U)
            FNT = FNR
            IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL 
     &           .OR. ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
            WRITE(*,1021)T, ET, CR, CI, CT, FNT
         ENDIF
      ENDDO
c
      ER = EIGENERG(DIM, N, NP, SHM, U, VPOT, WORK, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V, VPOT, WORK, VAR)
      ET  = (ER + EI)
      FNR = FCORREL(N, U, U) + FCORREL(N, V, V)
c     FNI must be igual to zero
c     FNI = FCORREL(N, U, V) - FCORREL(N, V, U)
      FNT = FNR
      IF(ABS((ET-E0)/E0)*1.0D+2.GT.TOL 
     &     .OR. ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
c
      RETURN
 1011 FORMAT(/,7X,A6,7X,A8,8X,A7,7X,A7,4X,A10,6X,A5)
 1021 FORMAT(1X,F14.8,3X,E12.6,2(4X,F10.7),4X,F10.8,4X,F8.6)
      END
