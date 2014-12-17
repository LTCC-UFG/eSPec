      SUBROUTINE S2PPABMN(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, EFC, 
     &     PRTCRL, TPABSOR, INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, 
     &     NPR, KP, NISPG, KL, NPES, E0, T0, TD, TP, OMG, SNI, DT, 
     &     TI, TF, TOL, DELQ, QC, AA, NP, ND, GAMMA, SH, SHM, 
     &     U1, U2, V1, V2, XI, XP, XF, DM, VPOT, VABC, WORK, VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL, TPABSOR
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR, KP
      INTEGER       NISPG, KL, NPES
      REAL*8        DT, TI, TF, TOL, E0, T0, TD, TP, OMG, SNI
      REAL*8        DELQ, QC, AA
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*), ND(*)
      REAL*8        GAMMA(*), SH(*), SHM(*), U1(*), U2(*), V1(*), V2(*)
      REAL*8        XI(*), XP(*), DM(*), XF(*), VPOT(*), VABC(*)
      REAL*8        WORK(*), VAR(*), LNZVC(MXDCT,*)
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
*     (20/06/2005) First version written by Freddy based on PSOD
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, THREE, FOUR, TWOPI, SHBAR
      PARAMETER(
     &     ZERO  = +0.0D+0, 
     &     ONE   = +1.0D+0, 
     &     TWO   = +2.0D+0,
     &     THREE = +3.0D+0, 
     &     FOUR  = +4.0D+0, 
     &     TWOPI = +6.283185307179586476925287D+0,
     &     SHBAR = +6.58211889D-1 !ev*fs
     &     )
c     **
c     ** Local scalars 
      LOGICAL       DERIV, FLG
      CHARACTER*4   CHNUM, CHNUM2   
      INTEGER       I, J, MC, MC1, MC2, NP1, NN, NM, NNM1, TWON
      REAL*8        ET, ET0, FNT, FNT0, SHFS, T, EFAUX, CT, CST
      REAL*8        CSQC, XD, OMGT, THREECST, DTD3, CRB2A, CIB2A
      REAL*8        CRB3A, CIB3A, COMGT, SOMGT, EFXDM
c     ** Local arrays 
cdel      LOGICAL
      CHARACTER*15  NEWNAM2(NPES), NEWNAM3(NPES)
      CHARACTER*16  NEWNAM1(NPES)
cdel      INTEGER    
      REAL*8        FNR(NPES), FNTP(NPES), F0(NPES)!, FNI(NPES)
      REAL*8        CSGAMMA1(NPES), CSGAMMA2(NPES), SAUXR(NPES)
      REAL*8        CRB1(NPES-1), CRB2(NPES-1), CRB3(NPES-1)
      REAL*8        CIB1(NPES-1), CIB2(NPES-1), CIB3(NPES-1)
      REAL*8        SUMR(NPES-1), SUMI(NPES-1), SAUXI(NPES-1)
      REAL*8        CR(NPES), CI(NPES), CTP(NPES), ER(NPES), EI(NPES)
      REAL*8        U0(NPES*N), V0(NPES*N), V(NPES*N), ETP(NPES)
      REAL*8        HU1(NPES*N), HV1(NPES*N)
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
c      DATA ET0 / NPES* 0 /
c
      DERIV = .FALSE.
D      DERIV = .TRUE.
      NPES = 2
c
      INFO  = ZERO
      T     = TI
      SHFS  = 2.4188843243D-2   !fs/a.u.(t) [hbar/E_{h}]
      CST   = TWO*DT/SHFS
      DTD3  = DT/THREE
      MC1   = ZERO
      NP1   = N + ONE
      TWON  = TWO*N
      THREECST = THREE*CST
c
      DO I=1,NPES,1
         CSGAMMA2(I) = CST*GAMMA(I)
         CSGAMMA1(I) = THREE*CSGAMMA2(I)
      ENDDO
c
      IF(DERIV)THEN
         CSQC = CST*AA/(6.0D0*SH(1))
      ELSE
         CSQC = CST*AA
      ENDIF
c
      DO I=1,NPES*N,1
         V(I)  = VPOT(I)
         U0(I) = U1(I)
         U2(I) = U0(I) 
         V0(I) = V1(I)
         V2(I) = V0(I) 
      ENDDO
c
      IF(PRTEIGVC2)THEN
         OPEN(UNIT=22,STATUS='UNKNOWN',FILE='movie.gplt')
         WRITE(22,*)'set xrange[*:*]'
         WRITE(22,*)'set yrange[*:*]'
         WRITE(22,*)'set xlabel "P/a.u."'
         WRITE(22,*)'set ylabel "S/a.u."' 
      ENDIF
c 
      FNT0 = ZERO
      ET0 = ZERO
      DO I=1,NPES,1
         NN = (I - ONE)*N + ONE
         FNR(I) = FCORREL(N, U2(NN), U2(NN)) 
     &        + FCORREL(N, V2(NN), V2(NN))
c     FNI must be igual to zero
D         FNI(I) = FCORREL(N, U2(NN), V2(NN)) 
D     &        - FCORREL(N, V2(NN), U2(NN))
         FNTP(I) = FNR(I)**2
         FNT0 = FNT0 + FNTP(I)
c
         ER(I) = EIGENERG(DIM, N, NP, SHM, U2(NN), V(NN), WORK(NN), VAR)
         EI(I) = EIGENERG(DIM, N, NP, SHM, V2(NN), V(NN), WORK(NN), VAR)
         ET0  = ET0 + ER(I) + EI(I) 
c     ..   variables used in the integration (1/3 Simpson first step)
         IF(I.GT.ONE)THEN
            CRB1A = FCORREL(N, U2(1), U2(NN)) 
     &           + FCORREL(N, V2(1), V2(NN))
            CIB1A = FCORREL(N, U2(1), V2(NN)) 
     &           - FCORREL(N, V2(1), U2(NN))
            CRB1(I-1) = EFXDM*(CRB1A*COMGT - CIB1A*SOMGT)
            CIB1(I-1) = EFXDM*(CRB1A*SOMGT + CIB1A*COMGT)
         ENDIF
      ENDDO
c
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         CT = ZERO
         ET = ZERO
         FNT = ZERO
         DO I=1,NPES,1
            NN = (I - ONE)*N + ONE
            CR(I) = FCORREL(N, U0(NN), U2(NN)) 
     &           + FCORREL(N, V0(NN), V2(NN))
            CI(I) = FCORREL(N, U0(NN), V2(NN)) 
     &           - FCORREL(N, V0(NN), U2(NN))
            CTP(I) = CR(I)**2 + CI(I)**2
            CT = CT + CTP(I)
c
            FNR(I) = FCORREL(N, U2(NN), U2(NN)) 
     &           + FCORREL(N, V2(NN), V2(NN))
c     FNI must be igual to zero
D            FNI(I) = FCORREL(N, U2(NN), V2(NN)) 
D     &           - FCORREL(N, V2(NN), U2(NN))
            FNTP(I) = FNR(I)
            FNT = FNT + FNTP(I)
c
            ER(I) = EIGENERG(DIM, N, NP, SHM, U2(NN), VPOT(NN), 
     &           WORK(NN), VAR)
            EI(I) = EIGENERG(DIM, N, NP, SHM, V2(NN), VPOT(NN), 
     &           WORK(NN), VAR)
            ETP(I)  = ER(I) + EI(I)
            ET = ET + ETP(I)
         ENDDO
c
         WRITE(*,*)
         IF(DIM.EQ.'.1D')THEN
            WRITE(*,*)'Auto-correlation function:'
         ELSE
            WRITE(*,*)'Partial auto-correlation function:'
         ENDIF
         WRITE(*,*)'=================================='
         WRITE(*,1011)
         WRITE(*,1021)T, ET, ET0, (CTP(I), I=1,NPES,1), 
     &        (FNTP(I), I=1,NPES,1), FNT 
         MC1 = ZERO
         MC2 = ONE
         WRITE(CHNUM,'(I4)')MC2
         DO I=1,NPES,1
            WRITE(CHNUM2,'(I1)')I
            NNM1 = (I - ONE)*N
            NN = NNM1 + ONE
            NEWNAM3(I) = 'ReIm'//CHNUM2//'_000'//CHNUM(4:4)//'.dat'
            NEWNAM2(I) = 'veff'//CHNUM2//'_000'//CHNUM(4:4)//'.dat'
            NEWNAM1(I) = 'eigvc'//CHNUM2//'_000'//CHNUM(4:4)//'.dat'
c  .. PRTEIGVC2 was introduced here to print the stationary PES
            IF(PRTVEFF .OR. PRTEIGVC2)THEN
               CALL PRTPT(NEWNAM2(I), 9, ND, T, NP, XP, XI, SH, 
     &              VPOT(NN))
            ENDIF
            IF(PRTEIGVC2)THEN
               WRITE(22,*)'#set output "',NEWNAM1(I),'"'
               WRITE(22,1031)NEWNAM1(I),T,NEWNAM3(I)
               DO J=1,N,1
                  WORK(J) = ETP(I) + U2(J+NNM1)*U2(J+NNM1) 
     &                 + V2(J+NNM1)*V2(J+NNM1)
               ENDDO
               CALL PRTPT(NEWNAM1(I), 9, ND, T, NP, XP, XI, SH, WORK)
               CALL PRPT2(NEWNAM3(I), 8, ND, T, NP, XP, XI, SH, 
     &              U2(NN), V2(NN))
            ENDIF
         ENDDO
         IF(PRTPULS)THEN
            OPEN(UNIT=21,STATUS='UNKNOWN',FILE='pulses.dat')
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T, EFAUX 
         ENDIF 
      ENDIF
c     .. computation of the 
      DO I=1,NPES,1
         DO J=1,NPES,1
            LS = (I - ONE)*NPES + J
            IF(I.NE.J)THEN
               DO K=1,N,1  
                  L1 = (I - ONE)*N + K
                  L2 = (J - ONE)*N + K 
                  UU(L1) = UU(L1) + TDCS1(LS)*TRAS(J,L1)*V1(L2)
                  VV(L1) = VV(L1) - TDCS1(LS)*TRAS(J,L1)*U1(L2)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c
      DO I=1,NPES,1
         NN = (I - ONE)*N + ONE
         IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, (T+DT)/TWO, 
     &        OMG, SNI, KL, DM(NN), V(NN), VPOT(NN))
         CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, THREECST/TWO, 
     &        CSGAMMA1(I), CSGAMMA2(I), NP, SHM, U1(NN), U2(NN), 
     &        V1(NN), V2(NN), HU1(NN), HV1(NN), VPOT(NN), VABC(NN), 
     &        WORK, VAR)
      ENDDO
c
      IF(DERIV)THEN
         DO I=1,NPES,1
            DO J=1,NPES,1    
               LS = (I - ONE)*NPES + J
               IF(I.NE.J)THEN
c     K=1
                  L1 = (I - ONE)*N + 1
                  L2 = (J - ONE)*N + 1
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                 + V1(L2+1)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2)
     &                 + U1(L2+1)) + VV(L1)
c     K=2
                  L1 = (I - ONE)*N + K
                  L2 = (J - ONE)*N + K 
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                 + V1(L2+1) - V1(L2-1)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                 + U1(L2+1) - U1(L2-1)) + VV(L1)
c     K=3:N-2
                  DO K=3,N-2,1  
                     L1 = (I - ONE)*N + K
                     L2 = (J - ONE)*N + K 
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                    + V1(L2+1) - V1(L2-1) - V1(L2-2)) + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                    + U1(L2+1) - U1(L2-1) - U1(L2-2)) + VV(L1)
                  ENDDO
c     K=N-1
                  L1 = (I - ONE)*N + N - 1
                  L2 = (J - ONE)*N + N - 1
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                 + V1(L2+1) - V1(L2-1) - V1(L2-2)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                 + U1(L2+1) - U1(L2-1) - U1(L2-2)) + VV(L1)
c     K=N
                  L1 = (I - ONE)*N + N
                  L2 = (J - ONE)*N + N 
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(- V1(L2-1) 
     &                 - V1(L2-2)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(- U1(L2-1)
     &                 - U1(L2-2)) + VV(L1)
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO I=1,NPES,1
            DO J=1,NPES,1
               LS = (I - ONE)*NPES + J
               IF(I.NE.J)THEN
                  DO K=1,N,1  
                     L1 = (I - ONE)*N + K
                     L2 = (J - ONE)*N + K 
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*V1(L2)
     &                    + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*U1(L2)
     &                    + VV(L1)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF
c
      DO I=1,NPES,1
         DO J=1,NPES,1
            LS = (I - ONE)*NPES + J
            IF(I.NE.J)THEN
               DO K=1,N,1  
                  L1 = (I - ONE)*N + K
                  L2 = (J - ONE)*N + K 
                  UU(L1) = UU(L1) + TDCS1(LS)*TRAS(J,L1)*V1(L2)
                  VV(L1) = VV(L1) - TDCS1(LS)*TRAS(J,L1)*U1(L2)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
c
      DO I=1,NPES,1
         NN = (I - ONE)*N + ONE
         IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, 
     &        SNI, KL, DM(NN), V(NN), VPOT(NN))
         CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, THREECST/TWO, 
     &        CSGAMMA1(I), CSGAMMA2(I), NP, SHM, U1(NN), U2(NN), 
     &        V1(NN), V2(NN), HU1(NN), HV1(NN), VPOT(NN), VABC(NN), 
     &        WORK, VAR)
      ENDDO
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(DERIV)THEN
         DO I=1,NPES,1
            DO J=1,NPES,1    
               LS = (I - ONE)*NPES + J
               IF(I.NE.J)THEN
c     K=1
                  L1 = (I - ONE)*N + 1
                  L2 = (J - ONE)*N + 1
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                 + V1(L2+1)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2)
     &                 + U1(L2+1)) + VV(L1)
c     K=2
                  L1 = (I - ONE)*N + K
                  L2 = (J - ONE)*N + K 
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                 + V1(L2+1) - V1(L2-1)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                 + U1(L2+1) - U1(L2-1)) + VV(L1)
c     K=3:N-2
                  DO K=3,N-2,1  
                     L1 = (I - ONE)*N + K
                     L2 = (J - ONE)*N + K 
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                    + V1(L2+1) - V1(L2-1) - V1(L2-2)) + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                    + U1(L2+1) - U1(L2-1) - U1(L2-2)) + VV(L1)
                  ENDDO
c     K=N-1
                  L1 = (I - ONE)*N + N - 1
                  L2 = (J - ONE)*N + N - 1
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                 + V1(L2+1) - V1(L2-1) - V1(L2-2)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                 + U1(L2+1) - U1(L2-1) - U1(L2-2)) + VV(L1)
c     K=N
                  L1 = (I - ONE)*N + N
                  L2 = (J - ONE)*N + N 
                  U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(- V1(L2-1) 
     &                 - V1(L2-2)) + UU(L1)
                  V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(- U1(L2-1)
     &                 - U1(L2-2)) + VV(L1)
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO I=1,NPES,1
            DO J=1,NPES,1
               LS = (I - ONE)*NPES + J
               IF(I.NE.J)THEN
                  DO K=1,N,1  
                     L1 = (I - ONE)*N + K
                     L2 = (J - ONE)*N + K 
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*V1(L2)
     &                    + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*U1(L2)
     &                    + VV(L1)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF
c     ..   variables used in the integration (1/3 Simpson second step)
      DO I=1,NPES-1,1
         NM = I*NPES + ONE
         CRB2A = FCORREL(N, U2(1), U2(NM)) 
     &        + FCORREL(N, V2(1), V2(NM))
         CIB2A = FCORREL(N, U2(1), V2(NM)) 
     &        - FCORREL(N, V2(1), U2(NM))
         CRB2(I) = EFXDM*(CRB2A*COMGT - CIB2A*SOMGT)
         CIB2(I) = EFXDM*(CRB2A*SOMGT + CIB2A*COMGT)
      ENDDO
c     
      MC = ONE
      T = TI + DT
      IF(MC.EQ.LC)MC = ZERO
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &     PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.LC)THEN
         CT = ZERO
         ET = ZERO
         FNT = ZERO
         DO I=1,NPES,1
            NN = (I - ONE)*N + ONE
            CR(I) = FCORREL(N, U0(NN), U2(NN)) 
     &           + FCORREL(N, V0(NN), V2(NN))
            CI(I) = FCORREL(N, U0(NN), V2(NN)) 
     &           - FCORREL(N, V0(NN), U2(NN))
            CTP(I) = CR(I)**2 + CI(I)**2
            CT = CT + CTP(I)
c     
            FNR(I) = FCORREL(N, U2(NN), U2(NN)) 
     &           + FCORREL(N, V2(NN), V2(NN))
c     FNI must be igual to zero
D     FNI(I) = FCORREL(N, U2(NN), V2(NN)) 
D     &           - FCORREL(N, V2(NN), U2(NN))
            FNTP(I) = FNR(I)
            FNT = FNT + FNTP(I)
c     
            ER(I) = EIGENERG(DIM, N, NP, SHM, U2(NN), VPOT(NN), 
     &           WORK(NN), VAR)
            EI(I) = EIGENERG(DIM, N, NP, SHM, V2(NN), VPOT(NN), 
     &           WORK(NN), VAR)
            ETP(I)  = (ER(I) + EI(I))
            ET = ET + ETP(I)
         ENDDO
c     
         WRITE(*,1021)T, ET, ET0, (CTP(I), I=1,NPES,1), 
     &        (FNTP(I), I=1,NPES,1), FNT 
c     
         MC1 = ZERO
         MC2 = ONE
         WRITE(CHNUM,'(I4)')MC2
         DO I=1,NPES,1
            NNM1 =  (I - ONE)*N
            NN =  NNM1 + ONE
            WRITE(CHNUM2,'(I1)')NN
            NEWNAM3(I) = 'ReIm'//CHNUM2//'_000'//CHNUM(4:4)//'.dat'
            NEWNAM2(I) = 'veff'//CHNUM2//'_000'//CHNUM(4:4)//'.dat'
            NEWNAM1(I) = 'eigvc'//CHNUM2//'_000'//CHNUM(4:4)//'.dat'
            IF(PRTVEFF)THEN
               NN = (I - ONE)*N + ONE
               CALL PRTPT(NEWNAM2(I), 9, ND, T, NP, XP, XI, SH, 
     &              VPOT(NN))
            ENDIF    
            IF(PRTEIGVC2)THEN
               WRITE(22,*)'#set output "',NEWNAM1(I),'"'
               WRITE(22,1031)NEWNAM1(I),t,NEWNAM3(I)
               DO J=1,N,1
                  WORK(J) = ETP(I) + U2(J+NNM1)*U2(J+NNM1) 
     &                 + V2(J+NNM1)*V2(J+NNM1)
               ENDDO
               CALL PRTPT(NEWNAM1(I), 9, ND, T, NP, XP, XI, SH, WORK)
               CALL PRPT2(NEWNAM3(I), 8, ND, T, NP, XP, XI, SH, 
     &              U2(NN), V2(NN))
            ENDIF
         ENDDO
      ENDIF
      IF(PRTPULS)THEN
         EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
         WRITE(21,*)T, EFAUX 
      ENDIF 
c      
      DO I=1,NPES*N,1
         U1(I) = U0(I)
         V1(I) = V0(I)
      ENDDO
      DO I=1,NPES,1
         NN = (I - ONE)*N + ONE
         CALL AU(DIM, N, NP, SHM, VPOT(NN), V0(NN), HU1(NN), VAR)
         CALL AU(DIM, N, NP, SHM, VPOT(NN), U0(NN), HV1(NN), VAR)
         SUMR(I) = ZERO
         SUMI(I) = ZERO
      ENDDO
      FLG = .FALSE.
c
      DO T=TI+2*DT,TF+DT/2,DT
c
         EFX     = EXP(-((T - T0X)/TPA)**KL2)
         EFXDM   = EFX*DMX 
         CSEFXDM = CST*EFXDM
         OMGT    = CSOMEGA*T
         COMGT   = CSEFXDM*DCOS(OMGT)
         SOMGT   = CSEFXDM*DSIN(OMGT)
         COMGT1  = CSEFXDM*DCOS(OMGT)/TWO
         SOMGT1  = CSEFXDM*DSIN(OMGT)/TWO
         COMGT3  = THREE/TWO*CSEFXDM*DCOS(OMGT)
         SOMGT3  = THREE/TWO*CSEFXDM*DSIN(OMGT)
c     
         DO I=1,NPES*N,1
            UU(I) = ZERO
            VV(I) = ZERO
         ENDDO
         DO I=1,NPES,1
            DO J=1,NPES,1
               LS = (I - ONE)*NPES + J
               IF(I.NE.J)THEN
                  DO K=1,N,1  
                     L1 = (I - ONE)*N + K
                     L2 = (J - ONE)*N + K 
                     UU(L1) = UU(L1) + TDCS1(LS)*TRAS(J,L1)*V1(L2)
                     VV(L1) = VV(L1) - TDCS1(LS)*TRAS(J,L1)*U1(L2)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
c
         DO I=1,NPES,1
            NN = (I - ONE)*N + ONE
            IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, 
     &           SNI, KL, DM(NN), V(NN), VPOT(NN))
            CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST, THREECST, 
     &           CSGAMMA1(I), CSGAMMA2(I), NP, SHM, U1(NN), U2(NN), 
     &           V1(NN), V2(NN), HU1(NN), HV1(NN), VPOT(NN), VABC(NN), 
     &           WORK, VAR)
         ENDDO
c
         IF(DERIV)THEN
            DO I=1,NPES,1
               DO J=1,NPES,1    
                  LS = (I - ONE)*NPES + J
                  IF(I.NE.J)THEN
c     K=1
                     L1 = (I - ONE)*N + 1
                     L2 = (J - ONE)*N + 1
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                    + V1(L2+1)) + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2)
     &                    + U1(L2+1)) + VV(L1)
c     K=2
                     L1 = (I - ONE)*N + K
                     L2 = (J - ONE)*N + K 
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                    + V1(L2+1) - V1(L2-1)) + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                    + U1(L2+1) - U1(L2-1)) + VV(L1)
c     K=3:N-2
                     DO K=3,N-2,1  
                        L1 = (I - ONE)*N + K
                        L2 = (J - ONE)*N + K 
                        U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                       + V1(L2+1) - V1(L2-1) - V1(L2-2)) + UU(L1)
                        V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                       + U1(L2+1) - U1(L2-1) - U1(L2-2)) + VV(L1)
                     ENDDO
c     K=N-1
                     L1 = (I - ONE)*N + N - 1
                     L2 = (J - ONE)*N + N - 1
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(V1(L2+2)
     &                    + V1(L2+1) - V1(L2-1) - V1(L2-2)) + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(U1(L2+2) 
     &                    + U1(L2+1) - U1(L2-1) - U1(L2-2)) + VV(L1)
c     K=N
                     L1 = (I - ONE)*N + N
                     L2 = (J - ONE)*N + N 
                     U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*(- V1(L2-1) 
     &                    - V1(L2-2)) + UU(L1)
                     V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*(- U1(L2-1)
     &                    - U1(L2-2)) + VV(L1)
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO I=1,NPES,1
               DO J=1,NPES,1
                  LS = (I - ONE)*NPES + J
                  IF(I.NE.J)THEN
                     DO K=1,N,1  
                        L1 = (I - ONE)*N + K
                        L2 = (J - ONE)*N + K 
                        U2(L1) = U2(L1) + TDCS3(LS)*TRAS(LS,L1)*V1(L2)
     &                       + UU(L1)
                        V2(L1) = V2(L1) - TDCS3(LS)*TRAS(LS,L1)*U1(L2)
     &                       + VV(L1)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
c     ..   1/3 Simpson integration
         IF(FLG)THEN
            DO I=1,NPES-1,1
               NM = I*NPES + ONE
               CRB2A = FCORREL(N, U2(1), U2(NM)) 
     &              + FCORREL(N, V2(1), V2(NM))
               CIB2A = FCORREL(N, U2(1), V2(NM)) 
     &              - FCORREL(N, V2(1), U2(NM))
               CRB2(I) = EFXDM*(CRB2A*COMGT - CIB2A*SOMGT)
               CIB2(I) = EFXDM*(CRB2A*SOMGT + CIB2A*COMGT)
            ENDDO
            FLG = .FALSE.
         ELSE   
            DO I=1,NPES-1,1
               NM = I*NPES + ONE
               CRB3A = FCORREL(N, U2(1), U2(NM)) 
     &              + FCORREL(N, V2(1), V2(NM))
               CIB3A = FCORREL(N, U2(1), V2(NM)) 
     &              - FCORREL(N, V2(1), U2(NM))
               CRB3(I) = EFXDM*(CRB3A*COMGT - CIB3A*SOMGT)
               CIB3(I) = EFXDM*(CRB3A*SOMGT + CIB3A*COMGT)
c
               SAUXR(I) = CRB1(I) + FOUR*CRB2(I) + CRB3(I)
               SUMR(I) = SUMR(I) + DTD3*SAUXR(I)
               CRB1(I) = CRB3(I)
c
               SAUXI(I) = CIB1(I) + FOUR*CIB2(I) + CIB3(I)
               SUMI(I) = SUMI(I) + DTD3*SAUXI(I)
               CIB1(I) = CIB3(I)
            ENDDO
            FLG = .TRUE.
         ENDIF
c     
         MC = MC + ONE
         IF(MC.EQ.LC)MC = ZERO
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
            CT = ZERO
            ET = ZERO
            FNT = ZERO
            DO I=1,NPES,1
               NN = (I - ONE)*N + ONE
               CR(I) = FCORREL(N, U0(NN), U2(NN)) 
     &              + FCORREL(N, V0(NN), V2(NN))
               CI(I) = FCORREL(N, U0(NN), V2(NN)) 
     &              - FCORREL(N, V0(NN), U2(NN))
               CTP(I) = CR(I)**2 + CI(I)**2
               CT = CT + CTP(I)
c     
               FNR(I) = FCORREL(N, U2(NN), U2(NN)) 
     &              + FCORREL(N, V2(NN), V2(NN))
c     FNI must be igual to zero
D               FNI(I) = FCORREL(N, U2(NN), V2(NN)) 
D     &              - FCORREL(N, V2(NN), U2(NN))
               FNTP(I) = FNR(I)
               FNT = FNT + FNTP(I)
c     
               ER(I) = EIGENERG(DIM, N, NP, SHM, U2(NN), VPOT(NN), 
     &              WORK(NN), VAR)
               EI(I) = EIGENERG(DIM, N, NP, SHM, V2(NN), VPOT(NN), 
     &              WORK(NN), VAR)
               ETP(I)  = (ER(I) + EI(I))
               ET = ET + ETP(I)
            ENDDO
c     
            WRITE(*,1021)T, ET, ET0, (CTP(I), I=1,NPES,1), 
     &           (FNTP(I), I=1,NPES,1), FNT 
c     
            IF(ABSORB)THEN
               CONTINUE
            ELSE
               IF(ABS((FNT-FNT0)/FNT0)*1.0D+2.GT.TOL)INFO = ONE
            ENDIF
c
            IF(MC.EQ.ZERO) MC1 = MC1 + ONE
            IF(MC1.EQ.NSHOT)THEN
               MC1 = ZERO
               MC2 = MC2 + ONE
               WRITE(CHNUM,'(I4)')MC2
               DO I=1,NPES,1 
                  WRITE(CHNUM2,'(I1)')I
                  NNM1 = (I - ONE)*N
                  NN = NNM1 + ONE
                  IF(MC2.LT.10)THEN
                     NEWNAM3(I) = 'ReIm'//CHNUM2//'_000'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM2(I) = 'veff'//CHNUM2//'_000'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM1(I) = 'eigvc'//CHNUM2//'_000'//CHNUM(4:4)//
     &                    '.dat'
                  ELSEIF(MC2.LT.100)THEN
                     NEWNAM3(I) = 'ReIm'//CHNUM2//'_00'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM2(I) = 'veff'//CHNUM2//'_00'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM1(I) = 'eigvc'//CHNUM2//'_00'//CHNUM(4:4)//
     &                    '.dat'
                  ELSEIF(MC2.LT.1000)THEN
                     NEWNAM3(I) = 'ReIm'//CHNUM2//'_0'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM2(I) = 'veff'//CHNUM2//'_0'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM1(I) = 'eigvc'//CHNUM2//'_0'//CHNUM(4:4)//
     &                    '.dat'
                  ELSEIF(MC2.LT.10000)THEN
                     NEWNAM3(I) = 'ReIm'//CHNUM2//'_'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM2(I) = 'veff'//CHNUM2//'_'//CHNUM(4:4)//
     &                    '.dat'
                     NEWNAM1(I) = 'eigvc'//CHNUM2//'_'//CHNUM(4:4)//
     &                    '.dat'
                  ELSE
                     WRITE(*,*)'Too many files to be printed! Stopping'
                     STOP
                  ENDIF
                  IF(PRTVEFF)THEN
                     CALL PRTPT(NEWNAM2(I), 9, ND, T, NP, XP, XI, SH, 
     &                    VPOT(NN))
                  ENDIF
                  IF(PRTEIGVC2)THEN
                     WRITE(22,*)'#set output "',NEWNAM1(I),'"'
                     WRITE(22,1031)NEWNAM1(I), T, NEWNAM3(I)
                     DO J=1,N,1
                        WORK(J) = ETP(I) + U2(J+NNM1)*U2(J+NNM1) 
     &                       + V2(J+NNM1)*V2(J+NNM1)
                     ENDDO
                     CALL PRTPT(NEWNAM1(I), 9, ND, T, NP, XP, XI, SH, 
     &                    WORK)
                     CALL PRPT2(NEWNAM3(I), 8, ND, T, NP, XP, XI, SH, 
     &                    U2(NN), V2(NN))
                  ENDIF
               ENDDO
               IF(PRTPULS)THEN
                  EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
                  WRITE(21,*)T, EFAUX 
               ENDIF 
            ENDIF
         ENDIF
      ENDDO
c                        
      IF(PRTCRL(1:7).EQ.'.LASTWP')THEN
         DO I=1,NPES,1
            WRITE(CHNUM2,'(I1)')I
            NEWNAM1(I) = 'eigvc'//CHNUM2//'_fwp.dat'
            NEWNAM3(I) = 'ReIm'//CHNUM2//'_fwp.dat'
            NNM1 = (I - ONE)*N
            NN = NNM1 + ONE
            DO J=1,N,1
               WORK(J) = U2(J+NNM1)*U2(J+NNM1) + V2(J+NNM1)*V2(J+NNM1) 
            ENDDO
            CALL PRTPT(NEWNAM1(I), 9, ND, T-DT, NP, XP, XI, SH, WORK)
            CALL PRPT2(NEWNAM3(I), 8, ND, T-DT, NP, XP, XI, SH, U2(NN),
     &        V2(NN))
         ENDDO
      ENDIF
c
c     ..
 1011 FORMAT(/,8X,'*t/fs*',6X,'*ET/a.u.*',6X,'*E0/a.u.*',6X,'*CIA*',
     &     6X,'*CRA*',6X,'*CTA*',8X,'*CIB*',6X,'*CRB*',6X,'*CTB*',
     &     8X,'*FTA*',8X,'*FTB*',8X,'*FT*')
 1021 FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),E10.4,3X,2(F8.5,3X),
     &     E10.4,3X,2(E10.4,3X),F8.6)
 1031 FORMAT('plot "',A15,'" title `',F9.4,1X,
     &     'fs` with lines lw 2.0, ','"',A14,
     &     '" notitle with lines lw 2.1, ',
     &     '"veffA_0001.dat" notitle with lines lt 3 lw 2.2')
      RETURN
      END
