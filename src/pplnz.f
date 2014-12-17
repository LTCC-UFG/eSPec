      SUBROUTINE PPLNZ(ABSORB, CHANGE, PRTVEFF, PRTEIGVC2, PRTPULS, 
     &     DIM, EFC, PRTCRL, TPABSOR, INFO, LC, LMTREORT, MP, MPI, MPR, 
     &     MXDCT, N, NSHOT, KL, NPR, KP, ABSTOL, DT, TI, TF, TOL, E0, 
     &     T0, TD, TP, OMG, SNI, NP, ND, IWORK, SH, SHM, VPOT, VABC, U, 
     &     V, XI, XF, XP, DM, WU, WV, WP, WORK, VAR, LNZVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, CHANGE, PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL,  TPABSOR          
      INTEGER       INFO, LC, LMTREORT, MP, MPI, MPR, MXDCT, N, NSHOT
      INTEGER       KL, NPR, KP, K, NF
      REAL*8        ABSTOL, DT, TI, TF, TOL, E0, T0, TD, TP, OMG 
      REAL*8        SNI
c     **
c     ** Array arguments
cdel      LOGICAL
      INTEGER       NP(*), ND(*), IWORK(*)
      REAL*8        SH(*), SHM(*), VPOT(*), U(*), V(*), XI(*), XF(*)
      REAL*8        XP(*), DM(*), WU(*), WV(*), WP(*), WORK(*), VABC(*)
      REAL*8        VAR(*), LNZVC(MXDCT,*)
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
      CHARACTER*4   CHNUM
      CHARACTER*14  NEWNAM2, NEWNAM3    
      CHARACTER*15  NEWNAM1
      INTEGER       I, MC, NC, MC1, MC2
      REAL*8        CR, CI, EINI, ER, EI, ET, F0, FNR, FNI, FNT
      REAL*8        SHFS, T, EFAUX, ER0, EI0, ET0  
c     **
c     ** Local arrays
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
      REAL*8        U0(N), V0(N), VP(N), W1(N), W2(N), W3(N), W4(N) 
      REAL*8        W5(N), W6(N), CT(7) 
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
      REAL*8        FCORREL, EIGENERG, EF, CST
c     **
c     ** External subroutines 
      EXTERNAL      PLNZAUX, PRTPT
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT     
c     .. Start program
c     ..
c     .. Initiate valors
      INFO = ZERO
      NC = ONE
      T = TI
      SHFS = 2.4188843243D-2 !fs/a.u.(t) [hbar/E_{h}]
      CST = DT/SHFS
      MC1 = ZERO
      MC = ZERO
      DO I=1,N,1
         VP(I) = VPOT(I)
         U0(I) = U(I)
         V0(I) = V(I)
         WU(I) = ZERO
         WV(I) = ZERO
         WP(I) = ZERO
         WORK(I) = ZERO
         IF(NPR-1.GE.ZERO)W1(I) = LNZVC(I,KP+1)
         IF(NPR-2.GE.ZERO)W2(I) = LNZVC(I,KP+2)
         IF(NPR-3.GE.ZERO)W3(I) = LNZVC(I,KP+3)
         IF(NPR-4.GE.ZERO)W4(I) = LNZVC(I,KP+4)
         IF(NPR-5.GE.ZERO)W5(I) = LNZVC(I,KP+5)
         IF(NPR-6.GE.ZERO)W6(I) = LNZVC(I,KP+6)
      ENDDO
c 
c      CR = FCORREL(N, U0, U) + FCORREL(N, V0, V)
c      CI = FCORREL(N, U0, V) - FCORREL(N, V0, U)
c      LNZVC(NC,MPR) = CR
c      LNZVC(NC,MPI) = CI
c      FNR = FCORREL(N, U, U) + FCORREL(N, V, V)
c      FNI = FCORREL(N, U, V) - FCORREL(N, V, U)
c      FNT = SQRT(FNR**2 + FNI**2)
c      F0 = FNT
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         IF(NPR-2.GE.ZERO)THEN
            CR = FCORREL(N, W1, U)
            CI = FCORREL(N, W1, V)
            CT(2) = CR**2 + CI**2
         ENDIF
         IF(NPR-3.GE.ZERO)THEN
            CR = FCORREL(N, W2, U)
            CI = FCORREL(N, W2, V)
            CT(3) = CR**2 + CI**2
         ENDIF
         IF(NPR-4.GE.ZERO)THEN
            CR = FCORREL(N, W3, U)
            CI = FCORREL(N, W3, V)
            CT(4) = CR**2 + CI**2
         ENDIF
         IF(NPR-5.GE.ZERO)THEN
            CR = FCORREL(N, W4, U)
            CI = FCORREL(N, W4, V)
            CT(5) = CR**2 + CI**2
         ENDIF
         IF(NPR-6.GE.ZERO)THEN
            CR = FCORREL(N, W5, U)
            CI = FCORREL(N, W5, V)
            CT(6) = CR**2 + CI**2
         ENDIF
         IF(NPR-7.GE.ZERO)THEN
            CR = FCORREL(N, W6, U)
            CI = FCORREL(N, W6, V)
            CT(7) = CR**2 + CI**2
         ENDIF
         CR = FCORREL(N, U0, U) + FCORREL(N, V0, V)
         CI = FCORREL(N, U0, V) - FCORREL(N, V0, U)
         LNZVC(NC,MPR) = CR
         LNZVC(NC,MPI) = CI
         CT(1) = CR**2 + CI**2 
c    
         FNR = FCORREL(N, U, U) + FCORREL(N, V, V)
         FNT = FNR**2 
         F0 = FNT
c
         ER = EIGENERG(DIM, N, NP, SHM, U, VPOT, WORK, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V, VPOT, WORK, VAR)
         ET  = (ER + EI)
         ER0 = EIGENERG(DIM, N, NP, SHM, U, VP, WORK, VAR)
         EI0 = EIGENERG(DIM, N, NP, SHM, V, VP, WORK, VAR)
         ET0  = (ER0 + EI0)
         EINI = ET
c
         WRITE(*,*)
         IF(DIM.EQ.'.1D')THEN
            WRITE(*,*)'Auto-correlation function:'
         ELSE
            WRITE(*,*)'Partial auto-correlation function:'
         ENDIF
         WRITE(*,*)'=================================='
         WRITE(*,1011)('*|P',INT(I-ONE),'|^2*', I=1,NPR,1),
     &        '*|c',0,'|^2*'         
         WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR,1), FNT
         MC1 = ZERO
         MC2 = ONE
         WRITE(CHNUM,'(I4)')MC2
         NEWNAM3 = 'ReIm_000'//CHNUM(4:4)//'.dat'
         NEWNAM2='veff_000'//CHNUM(4:4)//'.dat'
         NEWNAM1='eigvc_000'//CHNUM(4:4)//'.dat'
         IF(PRTVEFF .OR. PRTEIGVC2)THEN
            CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT)
         ENDIF
         IF(PRTEIGVC2)THEN
            DO I=1,N,1
               WORK(I) = (U(I)*U(I)+V(I)*V(I))**2
     &              + ET0
            ENDDO
            OPEN(UNIT=22,STATUS='UNKNOWN',FILE='movie.gplt')
            WRITE(22,*)'set xrange[*:*]'
            WRITE(22,*)'set yranginput.spce[*:*]'
            WRITE(22,*)'set xlabel "r/a.u."'
            WRITE(22,*)'set ylabel "Energy/a.u"'
            WRITE(22,*)'#set output "',NEWNAM1,'"'
            WRITE(22,1031)NEWNAM1,t
            CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
            CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, U, V)
         ENDIF
         IF(PRTPULS)THEN
            OPEN(UNIT=21,STATUS='UNKNOWN',FILE='pulse.dat')
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T, EFAUX/SQRT(E0)/2.1132D-13
         ENDIF 
      ENDIF
c
      DO T=TI+DT,TF+DT/2,DT
c      T = TI 
c      NF = NINT((TF+DT/2 - T)/DT)
c      DO K=1,NF,1
c         T = T + DT
c 
         CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, DM, VP, 
     &        VPOT)
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
            CT(1) = CR**2 + CI**2     
            LNZVC(NC,MPR) = CR
            LNZVC(NC,MPI) = CI
         ENDIF
c
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
c    
            IF(NPR-2.GE.ZERO)THEN
               CR = FCORREL(N, W1, U)
               CI = FCORREL(N, W1, V)
               CT(2) = CR**2 + CI**2
            ENDIF
            IF(NPR-3.GE.ZERO)THEN
               CR = FCORREL(N, W2, U)
               CI = FCORREL(N, W2, V)
               CT(3) = CR**2 + CI**2
            ENDIF
            IF(NPR-4.GE.ZERO)THEN
               CR = FCORREL(N, W3, U)
               CI = FCORREL(N, W3, V)
               CT(4) = CR**2 + CI**2
            ENDIF
            IF(NPR-5.GE.ZERO)THEN
               CR = FCORREL(N, W4, U)
               CI = FCORREL(N, W4, V)
               CT(5) = CR**2 + CI**2
            ENDIF
            IF(NPR-6.GE.ZERO)THEN
               CR = FCORREL(N, W5, U)
               CI = FCORREL(N, W5, V)
               CT(6) = CR**2 + CI**2
            ENDIF
            IF(NPR-7.GE.ZERO)THEN
               CR = FCORREL(N, W6, U)
               CI = FCORREL(N, W6, V)
               CT(7) = CR**2 + CI**2
            ENDIF
            CR = FCORREL(N, U0, U) + FCORREL(N, V0, V)
            CI = FCORREL(N, U0, V) - FCORREL(N, V0, U)
            CT(1) = CR**2 + CI**2
            ER = EIGENERG(DIM, N, NP, SHM, U, VPOT, WORK, VAR)
            EI = EIGENERG(DIM, N, NP, SHM, V, VPOT, WORK, VAR)
            ET  = (ER + EI)
            ER0 = EIGENERG(DIM, N, NP, SHM, U, VP, WORK, VAR)
            EI0 = EIGENERG(DIM, N, NP, SHM, V, VP, WORK, VAR)
            ET0  = (ER0 + EI0)
            FNR = FCORREL(N, U, U) + FCORREL(N, V, V)
c     FNI must be igual to zero
c            FNI = FCORREL(N, U, V) - FCORREL(N, V, U)
            FNT = FNR
c               
            IF(ABSORB)THEN
               CONTINUE
            ELSE
               IF(ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
            ENDIF
c
            WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR,1), FNT
c
            IF(MC.EQ.ZERO) MC1 = MC1 + ONE
            IF(MC1.EQ.NSHOT)THEN
               MC1 = ZERO
               MC2 = MC2 + ONE
               WRITE(CHNUM,'(I4)')MC2
               IF(MC2.LT.10)THEN
                  NEWNAM3 = 'ReIm_000'//CHNUM(4:4)//'.dat'
                  NEWNAM2 = 'veff_000'//CHNUM(4:4)//'.dat'
                  NEWNAM1 = 'eigvc_000'//CHNUM(4:4)//'.dat'
               ELSEIF(MC2.LT.100)THEN
                  NEWNAM3 = 'ReIm_00'//CHNUM(3:4)//'.dat'
                  NEWNAM2 = 'veff_00'//CHNUM(3:4)//'.dat'
                  NEWNAM1 = 'eigvc_00'//CHNUM(3:4)//'.dat'
               ELSEIF(MC2.LT.1000)THEN
                  NEWNAM3 = 'ReIm_0'//CHNUM(2:4)//'.dat'
                  NEWNAM2 = 'veff_0'//CHNUM(2:4)//'.dat'
                  NEWNAM1 = 'eigvc_0'//CHNUM(2:4)//'.dat'
               ELSEIF(MC2.LT.10000)THEN
                  NEWNAM3 = 'ReIm_'//CHNUM(1:4)//'.dat'
                  NEWNAM2 = 'veff_'//CHNUM(1:4)//'.dat'
                  NEWNAM1 = 'eigvc_'//CHNUM(1:4)//'.dat'
               ELSE
                  WRITE(*,*)'Too many files to be pinted! Stopping'
                  STOP
               ENDIF
               IF(PRTVEFF)THEN
                  CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT)
               ENDIF
               IF(PRTEIGVC2)THEN
                  DO I=1,N,1
                     WORK(I) = U(I)*U(I)+V(I)*V(I) + ET0
                  ENDDO
                  WRITE(22,*)'#set output "',NEWNAM1,'"'
                  WRITE(22,1031)NEWNAM1,t
                  CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
                  CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, U, V)
               ENDIF
            ENDIF
            IF(PRTPULS)THEN
               EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
               WRITE(21,*)T,EFAUX  !/SQRT(E0)/2.1132D-13
            ENDIF 
         ENDIF
c     
      ENDDO
c                        
      IF(PRTCRL(1:7).EQ.'.LASTWP')THEN
         NEWNAM1 = 'eigvc_fwp.dat'
         NEWNAM3 = 'ReIm_fwp.dat'
         CALL PRTPT(NEWNAM1, 9, ND, T-DT, NP, XP, XI, SH, WORK)
         CALL PRPT2(NEWNAM3, 8, ND, T-DT, NP, XP, XI, SH, U, V)
      ENDIF
c     ..
 1011 FORMAT(/,8X,'*t/fs*',5X,'*ET/a.u.*',6X,'*E0/a.u.*',6X,'*Pr(t)*'
     &     ,4X,'*Pi(t)*',3X,7(A3,I1,A4,3X))
 1021 FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),10(F8.6,3X))
 1031 FORMAT('plot "',A14,'" title `',F9.4,1X,
     &     'fs` with lines lw 2.0, ',
     &     '"veff_0001.dat" notitle with lines lt 3 lw 2.1')
      RETURN
      END


