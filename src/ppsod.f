      SUBROUTINE PPSOD(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS,  DIM, EFC, 
     &     PRTCRL, TPABSOR, INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, 
     &     NPR, KP, NISPG, E0, T0, TD, TP, OMG, SNI, KL, DT, TI, TF, 
     &     TOL, NP, ND, SH, SHM, U1, U2, V1, V2, XI, XP, XF, DM, VPOT, 
     &     VABC, WORK, VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL, TPABSOR
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR, KP
      INTEGER       NISPG, KL, NF
      REAL*8        DT, TI, TF, TOL, E0, T0, TD, TP, OMG, SNI
c     **
c     ** Array arguments
cdel      LOGICAL 
cdel      CHARACTER*(*)
      INTEGER       NP(*), ND(*)
      REAL*8        SH(*), SHM(*), U1(*), U2(*), V1(*), V2(*), XI(*)
      REAL*8        XP(*), DM(*), XF(*), VPOT(*), VABC(*), WORK(*)
      REAL*8        VAR(*), LNZVC(MXDCT,*)
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
     &     TWOPI = +6.283185307179586476925287D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
      CHARACTER*4   CHNUM   
      CHARACTER*14  NEWNAM2, NEWNAM3
      CHARACTER*15  NEWNAM1
      INTEGER       I, MC, NC, MC1, MC2, K, KST
      REAL*8        CR, CI, CST, EINI, ER, EI, ET, ER0, EI0, ET0, F0
      REAL*8        FNR, FNT, SHFS, T, EFAUX!, FNI
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*3
cdel      INTEGER       
      REAL*8        U0(N), V0(N), VP(N), W1(N), W2(N), W3(N), W4(N) 
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
      SHFS = 2.4188843243D-2 !fs/a.u.(t) [hbar/E_{h}]
      CST = TWO*DT/SHFS
      MC1 = ZERO
      DO I=1,N,1
         VP(I) = VPOT(I)
         U0(I) = U1(I)
         U2(I) = U0(I) 
         V0(I) = V1(I)
         V2(I) = V0(I) 
         IF(NPR-1.GE.ZERO)W1(I) = LNZVC(I,KP+1)
         IF(NPR-2.GE.ZERO)W2(I) = LNZVC(I,KP+2)
         IF(NPR-3.GE.ZERO)W3(I) = LNZVC(I,KP+3)
         IF(NPR-4.GE.ZERO)W4(I) = LNZVC(I,KP+4)
         IF(NPR-5.GE.ZERO)W5(I) = LNZVC(I,KP+5)
         IF(NPR-6.GE.ZERO)W6(I) = LNZVC(I,KP+6)
      ENDDO
c 
D      CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
D      CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
D      LNZVC(NC,MPR) = CR
D      LNZVC(NC,MPI) = CI
D      FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
c     FNI must be igual to zero
D      FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
D      FNT = FNR
D      F0 = FNT
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         IF(NPR-1.GE.ZERO)THEN
            CR = FCORREL(N, W1, U2)
            CI = FCORREL(N, W1, V2)
            CT(2) = CR**2 + CI**2
         ENDIF
         IF(NPR-2.GE.ZERO)THEN
            CR = FCORREL(N, W2, U2)
            CI = FCORREL(N, W2, V2)
            CT(3) = CR**2 + CI**2
         ENDIF
         IF(NPR-3.GE.ZERO)THEN
            CR = FCORREL(N, W3, U2)
            CI = FCORREL(N, W3, V2)
            CT(4) = CR**2 + CI**2
         ENDIF
         IF(NPR-4.GE.ZERO)THEN
            CR = FCORREL(N, W4, U2)
            CI = FCORREL(N, W4, V2)
            CT(5) = CR**2 + CI**2
         ENDIF
         IF(NPR-5.GE.ZERO)THEN
            CR = FCORREL(N, W5, U2)
            CI = FCORREL(N, W5, V2)
            CT(6) = CR**2 + CI**2
         ENDIF
         IF(NPR-6.GE.ZERO)THEN
            CR = FCORREL(N, W6, U2)
            CI = FCORREL(N, W6, V2)
            CT(7) = CR**2 + CI**2
         ENDIF
         CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
         CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
         LNZVC(NC,MPR) = CR
         LNZVC(NC,MPI) = CI
         CT(1) = CR**2 + CI**2
c   
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
         FNT = FNR
         F0 = FNT
c
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
D         print *,ER
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
D         print *,EI
         ET  = (ER + EI)
         ER0 = EIGENERG(DIM, N, NP, SHM, U2, VP, WORK, VAR)
         EI0 = EIGENERG(DIM, N, NP, SHM, V2, VP, WORK, VAR)
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
         IF(NPR+KP-1.GT.10)THEN
            WRITE(*,1012)'*|P',NISPG-1,'|^2*', 
     &           ('*|P',INT(I-1+KP),'|^2*', I=1,NPR,1),
     &           '*|c',NISPG-1,'|^2*' 
         ELSE
            WRITE(*,1011)'*|P',NISPG-1,'|^2*', 
     &           ('*|P',INT(I-1+KP),'|^2*', I=1,NPR,1),
     &           ' |c',NISPG-1,'|^2 '
         ENDIF
         WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR+1,1), FNT
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
               WORK(I) = U2(I)*U2(I) + V2(I)*V2(I)
!     &              + ET0
            ENDDO
            OPEN(UNIT=22,STATUS='UNKNOWN',FILE='movie.gplt')
            WRITE(22,*)'set xrange[*:*]'
            WRITE(22,*)'set yrange[*:*]'
            WRITE(22,*)'set xlabel "P/a.u."'
            WRITE(22,*)'set ylabel "S/a.u."'
            WRITE(22,*)'#set output "',NEWNAM1,'"'
            WRITE(22,1031)NEWNAM1,t,NEWNAM3
            CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
            CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, U2, V2)
         ENDIF
         IF(PRTPULS)THEN
            OPEN(UNIT=21,STATUS='UNKNOWN',FILE='pulse.dat')
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T, EFAUX        !/SQRT(E0)/2.1132D-9
         ENDIF 
      ENDIF
c     .. First half step in the propagation (Taylor)
      CALL CHPOT(EFC, N, E0, T0, TD, TP, (T+DT)/TWO, OMG, SNI, KL, DM, 
     &     VP, VPOT)
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/4D0, NP, SHM, U1, U2, 
     &     V1, V2, VPOT, VABC, WORK, VAR)
c     .. Second half step in the propagation towards the point T+DT (SOD) 
      CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, SNI, KL, DM, VP, 
     &     VPOT)
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, NP, SHM, U1, U2, 
     &     V1, V2, VPOT, VABC, WORK, VAR)
c
      MC = ONE
      T = TI + DT
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &     PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.LC)THEN
         IF(NPR-1.GE.ZERO)THEN
            CR = FCORREL(N, W1, U2)
            CI = FCORREL(N, W1, V2)
c            write(11,*)T, CR, CI
            CT(2) = CR**2 + CI**2
         ENDIF
         IF(NPR-2.GE.ZERO)THEN
            CR = FCORREL(N, W2, U2)
            CI = FCORREL(N, W2, V2)
c            write(12,*)T, CR, CI
            CT(3) = CR**2 + CI**2
         ENDIF
         IF(NPR-3.GE.ZERO)THEN
            CR = FCORREL(N, W3, U2)
            CI = FCORREL(N, W3, V2)
            CT(4) = CR**2 + CI**2
         ENDIF
         IF(NPR-4.GE.ZERO)THEN
            CR = FCORREL(N, W4, U2)
            CI = FCORREL(N, W4, V2)
            CT(5) = CR**2 + CI**2
         ENDIF
         IF(NPR-5.GE.ZERO)THEN
            CR = FCORREL(N, W5, U2)
            CI = FCORREL(N, W5, V2)
            CT(6) = CR**2 + CI**2
         ENDIF
         IF(NPR-6.GE.ZERO)THEN
            CR = FCORREL(N, W6, U2)
            CI = FCORREL(N, W6, V2)
            CT(7) = CR**2 + CI**2
         ENDIF
         CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
         CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
c         write(10,*)T, CR, CI
         CT(1) = CR**2 + CI**2
         IF(MC.EQ.LC)THEN
            NC = NC + ONE
            MC = ZERO
            LNZVC(NC,MPR) = CR 
            LNZVC(NC,MPI) = CI
         ENDIF
         ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
         EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
         ET  = (ER + EI)
         ER0 = EIGENERG(DIM, N, NP, SHM, U2, VP, WORK, VAR)
         EI0 = EIGENERG(DIM, N, NP, SHM, V2, VP, WORK, VAR)
         ET0  = (ER0 + EI0)
         FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
c     FNI must be igual to zero
D         FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
         FNT = FNR
c
         WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR+1,1), FNT
c
c         IF(MC.EQ.ZERO) MC1 = MC1 + ONE
         IF(MC1.EQ.NSHOT)THEN
            MC1 = ZERO
            MC2 = MC2 + ONE
            WRITE(CHNUM,'(I4)')MC2
            NEWNAM3 = 'ReIm_000'//CHNUM(4:4)//'.dat'
            NEWNAM2 = 'veff_000'//CHNUM(4:4)//'.dat'
            NEWNAM1 = 'eigvc_000'//CHNUM(4:4)//'.dat'
            IF(PRTVEFF)THEN
               CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT)
            ENDIF
            IF(PRTEIGVC2)THEN
               DO I=1,N,1
                  WORK(I) = U2(I)*U2(I)+V2(I)*V2(I) !+ ET0
               ENDDO
               WRITE(22,*)'#set output "',NEWNAM1,'"'
               WRITE(22,1031)NEWNAM1,t,NEWNAM3
               CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
               CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, U2, V2)
            ENDIF
         ENDIF
         IF(PRTPULS)THEN
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T,EFAUX       !/SQRT(E0)/2.1132D-9
         ENDIF 
      ENDIF
c     
      DO I=1,N,1
         U1(I) = U0(I)
         V1(I) = V0(I)
      ENDDO
c
      DO T=TI+2*DT,TF+DT/2,DT
c      T = TI + DT
c      NF = NINT((TF+DT/2 - T)/DT)
c      DO K=1,NF,1
c         T = T + DT
c     .. Steps in the propagation towards t (SOD).
         CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, DM, VP, 
     &        VPOT)
         CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST, NP, SHM, U1, U2, V1, 
     &        V2, VPOT, VABC, WORK, VAR)
c     
         MC = MC + ONE
         IF(MC.EQ.LC)THEN
            NC = NC + ONE
            MC = ZERO
            CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
            CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
            CT(1) = CR**2 + CI**2     
            LNZVC(NC,MPR) = CR
            LNZVC(NC,MPI) = CI
         ENDIF
c
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
c    
            IF(NPR-1.GE.ZERO)THEN
               CR = FCORREL(N, W1, U2)
               CI = FCORREL(N, W1, V2)
c               write(11,*)T, CR, CI
               CT(2) = CR**2 + CI**2
            ENDIF
            IF(NPR-2.GE.ZERO)THEN
               CR = FCORREL(N, W2, U2)
               CI = FCORREL(N, W2, V2)
c               write(12,*)T, CR, CI
               CT(3) = CR**2 + CI**2
            ENDIF
            IF(NPR-3.GE.ZERO)THEN
               CR = FCORREL(N, W3, U2)
               CI = FCORREL(N, W3, V2)
               CT(4) = CR**2 + CI**2
            ENDIF
            IF(NPR-4.GE.ZERO)THEN
               CR = FCORREL(N, W4, U2)
               CI = FCORREL(N, W4, V2)
               CT(5) = CR**2 + CI**2
            ENDIF
            IF(NPR-5.GE.ZERO)THEN
               CR = FCORREL(N, W5, U2)
               CI = FCORREL(N, W5, V2)
               CT(6) = CR**2 + CI**2
            ENDIF
            IF(NPR-6.GE.ZERO)THEN
               CR = FCORREL(N, W6, U2)
               CI = FCORREL(N, W6, V2)
               CT(7) = CR**2 + CI**2
            ENDIF
            CR = FCORREL(N, U0, U2) + FCORREL(N, V0, V2)
            CI = FCORREL(N, U0, V2) - FCORREL(N, V0, U2)
c            write(10,*)T, CR, CI
            CT(1) = CR**2 + CI**2
            ER = EIGENERG(DIM, N, NP, SHM, U2, VPOT, WORK, VAR)
            EI = EIGENERG(DIM, N, NP, SHM, V2, VPOT, WORK, VAR)
            ET  = (ER + EI)
            ER0 = EIGENERG(DIM, N, NP, SHM, U2, VP, WORK, VAR)
            EI0 = EIGENERG(DIM, N, NP, SHM, V2, VP, WORK, VAR)
            ET0  = (ER0 + EI0)
            FNR = FCORREL(N, U2, U2) + FCORREL(N, V2, V2)
c     FNI must be igual to zero
D            FNI = FCORREL(N, U2, V2) - FCORREL(N, V2, U2)
            FNT = FNR
c
            IF(ABSORB)THEN
               CONTINUE
            ELSE
               IF(ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
            ENDIF
c
            WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR+1,1), FNT
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
                     WORK(I) = U2(I)*U2(I) + V2(I)*V2(I) !+ ET0
                  ENDDO
                  WRITE(22,*)'#set output "',NEWNAM1,'"'
                  WRITE(22,1031)NEWNAM1,t,NEWNAM3
                  CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
                  CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, U2, V2)
               ENDIF
            ENDIF
            IF(PRTPULS)THEN
               EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
               WRITE(21,*)T,EFAUX         !/SQRT(E0)/2.1132D-9
            ENDIF 
         ENDIF
c     
      ENDDO
c                        
      IF(PRTCRL(1:7).EQ.'.LASTWP')THEN
         NEWNAM1 = 'eigvc_fwp.dat'
         NEWNAM3 = 'ReIm_fwp.dat'
         DO I=1,N,1
            WORK(I) = U2(I)*U2(I) + V2(I)*V2(I) !+ ET0
         ENDDO
         CALL PRTPT(NEWNAM1, 9, ND, T-DT, NP, XP, XI, SH, WORK)
         CALL PRPT2(NEWNAM3, 8, ND, T-DT, NP, XP, XI, SH, U2, V2)
c      ELSEIF(PRTCRL(1:7).EQ.'.LASTWP')THEN
c         
      ENDIF
c     ..
 1011 FORMAT(/,8X,'*t/fs*',5X,'*ET/a.u.*',6X,'*E0/a.u.*',6X,'*Pr(t)*'
     &     ,4X,'*Pi(t)*',3X,8(A3,I1,A4,3X)) 
 1012 FORMAT(/,8X,'*t/fs*',5X,'*ET/a.u.*',6X,'*E0/a.u.*',6X,'*Pr(t)*'
     &     ,4X,'*Pi(t)*',3X,8(A3,I2,A4,3X))
 1021 FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),10(F8.6,3X))
 1031 FORMAT('plot "',A14,'" title `',F9.4,1X,
     &     'fs` with lines lw 2.0, ','"',A13,
     &     '" notitle with lines lw 2.1, ',
     &     '"veff_0001.dat" notitle with lines lt 3 lw 2.2')
      RETURN
      END
