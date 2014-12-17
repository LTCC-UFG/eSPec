      SUBROUTINE PPSOD(PRTVEFF, PRTEIGVC2, PRTPULS, DIM, EFC, PRTCRL, 
     &     INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR, E0, T0, TD, TP,
     &     OMG, SNI, DT, TI, TF, TOL, NP, ND, SH, SHM, U1, U2, V1, V2, 
     &     XI, XP, XF, DM, VPOT, WORK, VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR
      REAL*8        DT, TI, TF, TOL, E0, T0, TD, TP, OMG, SNI
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*), ND(*)
      REAL*8        SH(*), SHM(*), U1(*), U2(*), V1(*), V2(*), XI(*)
      REAL*8        XP(*), DM(*), XF(*), VPOT(*), WORK(*), VAR(*)
      REAL*8        LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
*     Wave function propagation by SOD aproach:
*     |p_(j+1) = |p^(j-1) + 2*i*dt*|H* |p_(j-1)/hbar (j = 0, 1, 2, ...). 
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
*     (21/03/2003) First version PSOD written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0,
     &     TWOPI = +6.283185307D0)
c     **
c     ** Local scalars 
cdel      LOGICAL
      CHARACTER*3   CHNUM   
      CHARACTER*14  NEWNAM2, NEWNAM3
      CHARACTER*15  NEWNAM1
      INTEGER       I, MC, NC, MC1, MC2
      REAL*8        CR, CI, CST, EINI, ER, EI, ET, ER0, EI0, ET0, F0
      REAL*8        FNR, FNI, FNT, SHFS, T, EFAUX
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*3
cdel      INTEGER       
      REAL*8        U0(N), V0(N), V(N), W1(N), W2(N), W3(N), W4(N) 
      REAL*8        W5(N), W6(N), CT(7)  
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
      REAL*8        FCORREL, EIGENERG, EF    
c     **
c     ** External subroutines 
      EXTERNAL      PSODAUX, PRTPT
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
      MC1 = ZERO
      DO I=1,N,1
         V(I) = VPOT(I)
         U0(I) = U1(I)
         U2(I) = U0(I) 
         V0(I) = V1(I)
         V2(I) = V0(I) 
         IF(NPR-2.GE.ZERO)W1(I) = LNZVC(I,2)
         IF(NPR-3.GE.ZERO)W2(I) = LNZVC(I,3)
         IF(NPR-4.GE.ZERO)W3(I) = LNZVC(I,4)
         IF(NPR-5.GE.ZERO)W4(I) = LNZVC(I,5)
         IF(NPR-6.GE.ZERO)W5(I) = LNZVC(I,6)
         IF(NPR-7.GE.ZERO)W6(I) = LNZVC(I,7)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         write(*,*)'u(i),v(i)',u2(i),v2(i)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ENDDO
c
      IF(NPR-2.GE.ZERO)THEN
         CR = FCORREL(N, W1, U2)
         CI = FCORREL(N, W1, V2)
         CT(2) = CR**2 + CI**2
      ENDIF
      IF(NPR-3.GE.ZERO)THEN
         CR = FCORREL(N, W2, U2)
         CI = FCORREL(N, W2, V2)
         CT(3) = CR**2 + CI**2
      ENDIF
      IF(NPR-4.GE.ZERO)THEN
         CR = FCORREL(N, W3, U2)
         CI = FCORREL(N, W3, V2)
         CT(4) = CR**2 + CI**2
      ENDIF
      IF(NPR-5.GE.ZERO)THEN
         CR = FCORREL(N, W4, U2)
         CI = FCORREL(N, W4, V2)
         CT(5) = CR**2 + CI**2
      ENDIF
      IF(NPR-6.GE.ZERO)THEN
         CR = FCORREL(N, W5, U2)
         CI = FCORREL(N, W5, V2)
         CT(6) = CR**2 + CI**2
      ENDIF
      IF(NPR-7.GE.ZERO)THEN
         CR = FCORREL(N, W6, U2)
         CI = FCORREL(N, W6, V2)
         CT(7) = CR**2 + CI**2
      ENDIF
      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
      CT(1) = CR**2 + CI**2
      LNZVC(NC,MPR) = CR
      LNZVC(NC,MPI) = CI
c
      ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
      EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
      ET  = (ER + EI)
      ER0 = EIGENERG(DIM, N, NP, SHM, U2, V, WORK, VAR)
      EI0 = EIGENERG(DIM, N, NP, SHM, V2, V, WORK, VAR)
      ET0  = (ER0 + EI0)
c        
      EINI = ET
      FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
      FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
      FNT = SQRT(FNR**2 + FNI**2)
      F0 = FNT
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
 1011    FORMAT(/,8X,'*t/fs*',5X,'*ET/a.u.*',6X,'*E0/a.u.*',6X,'*Pr(t)*'
     &        ,4X,'*Pi(t)*',3X,7(A3,I1,A4,3X))
 1021    FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),10(F8.6,3X))
         WRITE(*,*)
         WRITE(*,*)'Partial auto-correlation function:'
         WRITE(*,*)'=================================='
         WRITE(*,1011)('*|P',I-ONE,'|^2*', I=1,NPR,1),'*|c',0,'|^2*'         
         WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR,1), FNT
         MC1 = ONE
         MC2 = ONE
         WRITE(CHNUM,'(I3)')MC2
         NEWNAM2='veff_00'//CHNUM(3:3)//'.dat'
         NEWNAM1='eigvc_00'//CHNUM(3:3)//'.dat'
         IF(PRTVEFF)THEN
            CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT)
         ENDIF
         IF(PRTEIGVC2)THEN
            DO I=1,N,1
               WORK(I) = (U2(I)*U2(I)+V2(I)*V2(I))**2
     &              + (U2(I)*V2(I)+V2(I)*U2(I))**2
     &              + ET
            ENDDO
            CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
         ENDIF
         IF(PRTPULS)THEN
            OPEN(UNIT=8,STATUS='UNKNOWN',FILE='pulse.dat')
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, T)
            WRITE(8,*)T, EFAUX/E0/+2.29371276D+17
         ENDIF 
      ENDIF
c
      CALL CHPOT(EFC, N, E0, T0, TD, TP, (T+DT)/TWO, OMG, SNI, DM, V, 
     &     VPOT)
      CALL PSODAUX(DIM, N, CST/2, NP, SHM, U1, U2, V1, V2, VPOT, 
     &     WORK, VAR)
      CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, SNI, DM, V, 
     &     VPOT)
      CALL PSODAUX(DIM, N, CST/2, NP, SHM, U1, U2, V1, V2, VPOT, 
     &     WORK, VAR)
c
      MC = ONE
      T = TI + DT
      IF(NPR-2.GE.ZERO)THEN
         CR = FCORREL(N, W1, U2)
         CI = FCORREL(N, W1, V2)
         CT(2) = CR**2 + CI**2
      ENDIF
      IF(NPR-3.GE.ZERO)THEN
         CR = FCORREL(N, W2, U2)
         CI = FCORREL(N, W2, V2)
         CT(3) = CR**2 + CI**2
      ENDIF
      IF(NPR-4.GE.ZERO)THEN
         CR = FCORREL(N, W3, U2)
         CI = FCORREL(N, W3, V2)
         CT(4) = CR**2 + CI**2
      ENDIF
      IF(NPR-5.GE.ZERO)THEN
         CR = FCORREL(N, W4, U2)
         CI = FCORREL(N, W4, V2)
         CT(5) = CR**2 + CI**2
      ENDIF
      IF(NPR-6.GE.ZERO)THEN
         CR = FCORREL(N, W5, U2)
         CI = FCORREL(N, W5, V2)
         CT(6) = CR**2 + CI**2
      ENDIF
      IF(NPR-7.GE.ZERO)THEN
         CR = FCORREL(N, W6, U2)
         CI = FCORREL(N, W6, V2)
         CT(7) = CR**2 + CI**2
      ENDIF
      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
      CT(1) = CR**2 + CI**2
      IF(MC.EQ.LC)THEN
         NC = NC + ONE
         MC = ZERO
         LNZVC(NC,MPR) = CR 
         LNZVC(NC,MPI) = CI
      ENDIF
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &     PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
         ET  = (ER + EI)
         ER0 = EIGENERG(DIM, N, NP, SHM, U2, V, WORK, VAR)
         EI0 = EIGENERG(DIM, N, NP, SHM, V2, V, WORK, VAR)
         ET0  = (ER0 + EI0)
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
         FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
         FNT = SQRT(FNR**2 + FNI**2)
         WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR,1), FNT
         IF(MC.EQ.ZERO) MC1 = MC1 + ONE
         IF(MC1.EQ.NSHOT)THEN
            MC1 = ZERO
            MC2 = MC2 + ONE
            WRITE(CHNUM,'(I3)')MC2
            NEWNAM2 = 'veff_00'//CHNUM(3:3)//'.dat'
            NEWNAM1 = 'eigvc_00'//CHNUM(3:3)//'.dat'
            IF(PRTVEFF)THEN
               CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT)
            ENDIF
            IF(PRTEIGVC2)THEN
               DO I=1,N,1
                  WORK(I) = U2(I)*U2(I)+V2(I)*V2(I) + ET
               ENDDO
               CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
            ENDIF
         ENDIF
         IF(PRTPULS)THEN
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, T)
            WRITE(8,*)T,EFAUX/E0/+2.29371276D+17
         ENDIF 
      ENDIF
c
      DO I=1,N,1
         U1(I) = U0(I)
         V1(I) = V0(I)
      ENDDO
c
      DO T=TI+2*DT,TF+DT,DT
c     
         CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, DM, V, 
     &     VPOT)
         CALL PSODAUX(DIM, N, CST, NP, SHM, U1, U2, V1, V2, VPOT, 
     &        WORK, VAR) 
c    
         IF(NPR-2.GE.ZERO)THEN
            CR = FCORREL(N, W1, U2)
            CI = FCORREL(N, W1, V2)
            CT(2) = CR**2 + CI**2
         ENDIF
         IF(NPR-3.GE.ZERO)THEN
            CR = FCORREL(N, W2, U2)
            CI = FCORREL(N, W2, V2)
            CT(3) = CR**2 + CI**2
         ENDIF
         IF(NPR-4.GE.ZERO)THEN
            CR = FCORREL(N, W3, U2)
            CI = FCORREL(N, W3, V2)
            CT(4) = CR**2 + CI**2
         ENDIF
         IF(NPR-5.GE.ZERO)THEN
            CR = FCORREL(N, W4, U2)
            CI = FCORREL(N, W4, V2)
            CT(5) = CR**2 + CI**2
         ENDIF
         IF(NPR-6.GE.ZERO)THEN
            CR = FCORREL(N, W5, U2)
            CI = FCORREL(N, W5, V2)
            CT(6) = CR**2 + CI**2
         ENDIF
         IF(NPR-7.GE.ZERO)THEN
            CR = FCORREL(N, W6, U2)
            CI = FCORREL(N, W6, V2)
            CT(7) = CR**2 + CI**2
         ENDIF
         CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
         CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
         CT(1) = CR**2 + CI**2
c     
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
         ET  = (ER + EI)
         ER0 = EIGENERG(DIM, N, NP, SHM, U2, V, WORK, VAR)
         EI0 = EIGENERG(DIM, N, NP, SHM, V2, V, WORK, VAR)
         ET0  = (ER0 + EI0)
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
         FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
         FNT = SQRT(FNR**2 + FNI**2)
c     
         IF(ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
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
            WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR,1), FNT
            IF(MC.EQ.ZERO) MC1 = MC1 + ONE
            IF(MC1.EQ.NSHOT)THEN
               MC1 = ZERO
               MC2 = MC2 + ONE
               WRITE(CHNUM,'(I3)')MC2
               IF(MC2.LT.10)THEN
                  NEWNAM3 = 'ReIm_00'//CHNUM(3:3)//'.dat'
                  NEWNAM2 = 'veff_00'//CHNUM(3:3)//'.dat'
                  NEWNAM1 = 'eigvc_00'//CHNUM(3:3)//'.dat'
               ELSEIF(MC2.LT.100)THEN
                  NEWNAM3 = 'ReIm_0'//CHNUM(2:3)//'.dat'
                  NEWNAM2 = 'veff_0'//CHNUM(2:3)//'.dat'
                  NEWNAM1 = 'eigvc_0'//CHNUM(2:3)//'.dat'
               ELSE
                  NEWNAM3 = 'ReIm_'//CHNUM(1:3)//'.dat'
                  NEWNAM2 = 'veff_'//CHNUM(1:3)//'.dat'
                  NEWNAM1 = 'eigvc_'//CHNUM(1:3)//'.dat'
               ENDIF
               IF(PRTVEFF)THEN
                  CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT)
               ENDIF
               IF(PRTEIGVC2)THEN
                  DO I=1,N,1
                     WORK(I) = U2(I)*U2(I)+V2(I)*V2(I) + ET
                  ENDDO
                  CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
                  CALL PRPT2(NEWNAM3, 9, ND, T, NP, XP, XI, SH, U2, V2)
               ENDIF
            ENDIF
            IF(PRTPULS)THEN
               EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, T)
               WRITE(8,*)T,EFAUX/E0/+2.29371276D+17
            ENDIF 
         ENDIF
c     
      ENDDO
c     ..
      RETURN
      END
