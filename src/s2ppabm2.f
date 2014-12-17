      SUBROUTINE S2PPABM2(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, EFC, 
     &     PRTCRL, TPABSOR, INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, 
     &     NPR, KP, NISPG, KL, KLX, E0, T0, TD, TP, OMG, SNI, DT, 
     &     TI, TF, TOL, GAMMA, OMEGA, DMX, T0X, TPX, AA, QC, DELQ, NP, 
     &     ND, SH, SHM, U1, U2, V1, V2, XI, XP, XF, DM, VPOT, VABC, 
     &     WORK, SM, VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL, TPABSOR
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR, KP
      INTEGER       NISPG, KL, KLX
      REAL*8        DT, TI, TF, TOL, E0, T0, TD, TP, OMG, SNI
      REAL*8        GAMMA, OMEGA, DMX, T0X, TPX
      REAL*8        AA, QC, DELQ         
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*), ND(*)
      REAL*8        SH(*), SHM(*), U1(*), U2(*), V1(*), V2(*), XI(*)
      REAL*8        XP(*), DM(*), XF(*), VPOT(*), VABC(*), WORK(*)
      REAL*8        LNZVC(MXDCT,*), SM(*), VAR(*)
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
*     (22/03/2005) First version written by Freddy based on ppsod.f
*
c     ** 
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, THREE, FOUR, TWOPI, SHBAR, EHT
      PARAMETER(
     &     ZERO  = +0.0D+0, 
     &     ONE   = +1.0D+0, 
     &     TWO   = +2.0D+0,
     &     THREE = +3.0D+0, 
     &     FOUR  = +4.0D+0, 
     &     EHT   = +8.0D+0,
     &     TWOPI = +6.283185307179586476925287D+0,
     &     SHBAR = +6.58211889D-1 !ev*fs
     &     )
c     **
c     ** Local scalars 
      LOGICAL       FLG
      CHARACTER*4   CHNUM   
      CHARACTER*15  NEWNAM2, NEWNAM3, NEWNAM5, NEWNAM6
      CHARACTER*16  NEWNAM1, NEWNAM4
      INTEGER       I, J, MC, MC1, MC2, NP1, KL2, TWON, NN
      REAL*8        CST, EINI, ET, ET0, FNT, SHFS, T, EFAUX
      REAL*8        CRA, CIA, CTA, ERA, EIA, ETA, ER0A, EI0A, ET0A
      REAL*8        FNRA, FNTA, F0A
      REAL*8        CRB, CIB, CTB, ERB, EIB, ETB, ER0B, EI0B, ET0B
      REAL*8        FNRB, FNTB!, FNIA, FNIB 
      REAL*8        CSTD2, DTD3, EFX, EFXDM, SUMR, SUMI, SAUXR, SAUXI
      REAL*8        CRB1, CRB2, CRB3, CIB1, CIB2, CIB3, TPA, A1, OMGT
      REAL*8        COMGT, SOMGT, CSEFXDM, CSGAMMA1, CSGAMMA2, CSOMEGA
      REAL*8        THREECST, CRB2A, CIB2A, CRB3A, CIB3A
      REAL*8        COMGT1, SOMGT1, COMGT3, SOMGT3
      REAL*8        CSQC, CSQC1, CSQC2, XD, EP, CEP2, ANT, UAUX, VAUX
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*3
cdel      INTEGER       
      REAL*8        U0(2*N), V0(2*N), V(2*N), HU1(2*N), HV1(2*N), 
     &     HV2(2*N)
c      REAL*8        UU(2*N), VV(2*N)
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
      INFO  = ZERO
      T     = TI
      SHFS  = 2.4188843243D-2 !fs/a.u.(t) [hbar/E_{h}]
      CST   = DT/(TWO*SHFS)
      THREECST = THREE*CST
      CSGAMMA2 = CST*GAMMA!*0
      CSGAMMA1 = THREE*CSGAMMA2!*0
      CSQC  =  AA/(TWO*SH(1))/SM(1)
      CSQC  =  AA/(1.2D+1*SH(1))/SM(1)
      CSQC1 =  AA/SM(1)
      CSQC2 =  AA*AA/(2.0D0*SM(1))
      CSOMEGA = OMEGA
      MC1 = ZERO
d      write(*,*)CSQC,CSQC1,CSQC2
d      read(*,*)
c
      CRB1 = ZERO
      CRB2 = ZERO
      CRB3 = ZERO
      CIB1 = ZERO
      CIB2 = ZERO
      CIB3 = ZERO
c
      NP1   = N + ONE 
      TWON  = TWO*N
      CSTD2 = CST/TWO           ! DT/SHFS
      DTD3  = DT/THREE 
      KL2   = TWO*KLX
      A1    = ONE/KL2
      TPA   = TPX/(0.693147181**A1)
c
      DO I=1,2*N,1
         V(I) = VPOT(I)
         U0(I) = U1(I)
         U2(I) = U0(I) 
         V0(I) = V1(I)
         V2(I) = V0(I)
         HU1(I) = ZERO
         HV1(I) = ZERO 
         HV2(I) = ZERO
      ENDDO
c 
      FNRA = FCORREL(N, U2(1), U2(1)) 
     &     + FCORREL(N, V2(1), V2(1))
c     FNI must be igual to zero
D      FNIA = FCORREL(N, U2(1), V2(1)) 
D     &     - FCORREL(N, V2(1), U2(1))
      FNTA = FNRA !+ FNIA
      F0A = FNTA
c
      FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &     + FCORREL(N, V2(NP1), V2(NP1))
c     FNI must be igual to zero
D      FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &     - FCORREL(N, V2(NP1), U2(NP1))
      FNTB = FNRB !+ FNIB
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         CRA = FCORREL(N, U0(1), U2(1)) 
     &        + FCORREL(N, V0(1), V2(1))
         CIA = FCORREL(N, U0(1), V2(1)) 
     &        - FCORREL(N, V0(1), U2(1))
         CTA = CRA**2 + CIA**2
c
         ERA = EIGENERG(DIM, N, NP, SHM, U2(1), VPOT(1), WORK(1), VAR)
         EIA = EIGENERG(DIM, N, NP, SHM, V2(1), VPOT(1), WORK(1), VAR)
         ETA  = (ERA + EIA)
         ER0A = EIGENERG(DIM, N, NP, SHM, U2(1), V(1), WORK(1), VAR)
         EI0A = EIGENERG(DIM, N, NP, SHM, V2(1), V(1), WORK(1), VAR)
         ET0A  = (ER0A + EI0A)  
c
         CRB = FCORREL(N, U0(1), U2(NP1)) 
     &        + FCORREL(N, V0(1), V2(NP1))
         CIB = FCORREL(N, U0(1), V2(NP1)) 
     &        - FCORREL(N, V0(1), U2(NP1))  
         CTB = CRB**2 + CIB**2
c    
         ERB = EIGENERG(DIM, N, NP, SHM, U2(NP1), VPOT(NP1),
     &        WORK(NP1), VAR)
         EIB = EIGENERG(DIM, N, NP, SHM, V2(NP1), VPOT(NP1), 
     &        WORK(NP1), VAR)
         ETB  = (ERB + EIB)
         ER0B = EIGENERG(DIM, N, NP, SHM, U2(NP1), V(NP1), 
     &        WORK(NP1), VAR)
         EI0B = EIGENERG(DIM, N, NP, SHM, V2(NP1), V(NP1), 
     &        WORK(NP1), VAR)
         ET0B  = (ER0B + EI0B)
c
         ET = ETA + ETB
         ET0 = ET0A + ET0B
         FNT = FNTA + FNTB
         EINI = ET
c
         WRITE(*,*)
         IF(DIM.EQ.'.1D')THEN
            WRITE(*,*)'Auto-correlation function:'
            WRITE(*,*)'=========================='
         ELSE
            WRITE(*,*)'Partial auto-correlation function:'
            WRITE(*,*)'=================================='
         ENDIF
         WRITE(*,1011)
         WRITE(*,1021)T, ET, ET0, CRA, CIA, CTA, CTB, FNTA, FNTB, FNT 
         MC1 = ZERO
         MC2 = ONE
         WRITE(CHNUM,'(I4)')MC2
         NEWNAM6 = 'ReImB_000'//CHNUM(4:4)//'.dat'
         NEWNAM5 ='veffB_000'//CHNUM(4:4)//'.dat'
         NEWNAM4 ='eigvcB_000'//CHNUM(4:4)//'.dat'
         NEWNAM3 = 'ReImA_000'//CHNUM(4:4)//'.dat'
         NEWNAM2 ='veffA_000'//CHNUM(4:4)//'.dat'
         NEWNAM1 ='eigvcA_000'//CHNUM(4:4)//'.dat'
         IF(PRTVEFF .OR. PRTEIGVC2)THEN
            CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, VPOT(1))
            CALL PRTPT(NEWNAM5, 9, ND, T, NP, XP, XI, SH, VPOT(NP1))
         ENDIF
         IF(PRTEIGVC2)THEN
            OPEN(UNIT=22,STATUS='UNKNOWN',FILE='movie.gplt')
            WRITE(22,*)'set xrange[*:*]'
            WRITE(22,*)'set yrange[*:*]'
            WRITE(22,*)'set xlabel "P/a.u."'
            WRITE(22,*)'set ylabel "S/a.u."'
            WRITE(22,*)'#set output "',NEWNAM1,'"'
            WRITE(22,1031)NEWNAM1,t,NEWNAM3
c
            DO I=1,N,1
               WORK(I) = U2(I)*U2(I) + V2(I)*V2(I)
     &              + ET0A
            ENDDO
            CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
            CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, 
     &           U2(1), V2(1))
c            
            DO I=NP1,TWON,1
               WORK(I) = U2(I)*U2(I) + V2(I)*V2(I)
     &              + ET0B
            ENDDO
            CALL PRTPT(NEWNAM4, 9, ND, T, NP, XP, XI, SH, WORK)
            CALL PRPT2(NEWNAM6, 8, ND, T, NP, XP, XI, SH, 
     &           U2(NP1), V2(NP1))
         ENDIF
         IF(PRTPULS)THEN
            OPEN(UNIT=21,STATUS='UNKNOWN',FILE='pulses.dat')
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T, EFAUX 
         ENDIF 
      ENDIF
c     .. First half step in the propagation (Taylor)
      DO I=1,2,1
         NN = (I - ONE)*N + ONE
         IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT/TWO, 
     &        OMG, SNI, KL, DM(1), V(1), VPOT(1))
         CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST, NP, SHM, 
     &        U1(NN), U2(NN), V1(NN), V2(NN), VPOT(NN), VABC(NN), 
     &        WORK, VAR)
D         CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, THREECST/TWO, 
D     &        CSGAMMA1/TWO, CSGAMMA2/TWO, NP, SHM, U1(NN), U2(NN), 
D     &        V1(NN), V2(NN), HU1(NN), HV1(NN), VPOT(NN), VABC(NN), 
D     &        WORK)
      ENDDO
c
      DO I=2,N-1,1
         J = I + N
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = CSQC1*(XD - QC)/DELQ**2*EP
         U2(I) = U2(I) + CEP2*V1(I)
     &        + OMGT*(V1(J+1) - V1(J-1)) + ANT*V1(J)
         V2(I) = V2(I) - CEP2*U1(I)
     &        - OMGT*(U1(J+1) - U1(J-1)) - ANT*U1(J)
         U2(J) = U2(J) + CEP2*V1(J)
     &        - OMGT*(V1(I+1) - V1(I-1)) - ANT*V1(I)
         V2(J) = V2(J) - CEP2*U1(J)
     &        + OMGT*(U1(I+1) - U1(I-1)) + ANT*U1(I)
      ENDDO
c     .. Second half step in the propagation ABM      DO I=1,2,1
      DO I=1,2,1
         NN = (I - ONE)*N + ONE
         CALL AU(DIM, N, NP, SHM, VPOT(NN), V0(NN), HU1(NN), VAR)
         CALL AU(DIM, N, NP, SHM, VPOT(NN), U0(NN), HV1(NN), VAR)
      ENDDO
      DO I=1,N,1
         J = I + N
         HU1(I) = HU1(I) - GAMMA*U0(I)
     &        + CEP2*V0(I) + ANT*V0(J)
     &        + OMGT*(V0(J+1) - V0(J-1)) 
         HU1(J) = HU1(J) - GAMMA*U0(J)
     &        + CEP2*V0(J) - ANT*V0(I)
     &        - OMGT*(V0(I+1) - V0(I-1))
         HV1(I) = HV1(I) + GAMMA*V0(I) 
     &        - CEP2*U0(I) - ANT*U0(J)
     &        - OMGT*(U0(J+1) - U0(J-1))
         HV1(J) = HV1(J) + GAMMA*V0(J) 
     &        - CEP2*U0(J) + ANT*U0(I)
     &        + OMGT*(U0(I+1) - U0(I-1)) 
      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO I=1,2,1
         NN = (I - ONE)*N + ONE
         CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, SNI, KL, 
     &        DM(NN), V(NN), VPOT(NN))
         CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, THREECST/TWO, 
     &        CSGAMMA1/TWO, CSGAMMA2/TWO, NP, SHM, U1(NN), U2(NN), 
     &        V1(NN), V2(NN), HU1(NN), HV1(NN), VPOT(NN), VABC(NN), 
     &        WORK, VAR)
      ENDDO
c     
c      DO I=1,2,1
c         NN = (I - ONE)*N + ONE
c         CALL AU(DIM, N, NP, SHM, VPOT(NN), V2(NN), WORK(NN))
c      ENDDO
c      XD = XI(1)   
c      DO I=2,N-1,1
c         XD = XD + SH(1)
c         EP = EXP(-((XD - QC)/DELQ)**2)
c         CEP2 = CSQC2/TWO*EP*EP
c         OMGT = CSQC/TWO*EP
c         ANT  = CSQC1/TWO*(XD - QC)/DELQ**2*EP
c         J = I + N
c         UAUX  = U2(I)
c         U2(I) = UAUX + THREECST/TWO*WORK(I) - CSGAMMA1/TWO*UAUX
c     &        + THREECST/TWO*CEP2*V2(I) + THREECST/TWO*ANT*V2(J)
c     &        + THREECST/TWO*OMGT*(V2(J+1) - V2(J-1)) 
c     &        - CST/TWO*HU1(I) 
c         HU1(I) = WORK(I) - GAMMA*U1(I)
c     &        + CEP2*V1(I) + ANT*V1(J)
c     &        + OMGT*(V1(J+1) - V1(J-1)) 
c         HV2(I) = - CEP2*U1(I) - ANT*U1(J)
c     &        - OMGT*(U1(J+1) - U1(J-1)) 
c         U1(I) = UAUX
c     
c         UAUX  = U2(J)
c         U2(J) = UAUX + THREECST/TWO*WORK(J) - CSGAMMA1/TWO*UAUX
c     &        + THREECST/TWO*CEP2*V2(J) - THREECST/TWO*ANT*V2(I)
c     &        - THREECST/TWO*OMGT*(V2(I+1) - V2(I-1)) 
c     &        - CST/TWO*HU1(J) 
c         HU1(J) = WORK(J) - GAMMA*U1(J)
c     &        + CEP2*V1(J) - ANT*V1(I)
c     &        - OMGT*(V1(I+1) - V1(I-1))
c         HV2(J) = - CEP2*U1(J) + ANT*U1(I)
c     &        + OMGT*(U1(I+1) - U1(I-1)) 
c         U1(J) = UAUX
c      ENDDO
c     
c      DO I=1,2,1
c         NN = (I - ONE)*N + ONE
c         CALL AU(DIM, 2*N, NP, SHM, VPOT(NN), U1(NN), WORK(NN)) 
c      ENDDO
c      XD = XI(1)
c      DO I=2,N-1,1
c         XD = XD + SH(1)
c         EP = EXP(-((XD - QC)/DELQ)**2)
c         CEP2 = CSQC2*EP*EP
c         OMGT = CSQC*EP
c         ANT  = CSQC1*(XD - QC)/DELQ**2*EP
c         J = I + N
c         VAUX = V2(I)
c         V2(I) = VAUX - THREECST/TWO*WORK(I) - CSGAMMA1/TWO*VAUX
c     &        - THREECST/TWO*CEP2*U1(I) - THREECST/TWO*ANT*U1(J)
c     &        - THREECST/TWO*OMGT*(U1(J+1) - U1(J-1)) 
c     &        + CST/TWO*HV1(I)
c         HV1(I) = WORK(I) + GAMMA*V1(I) + HV2(I)
c         V1(I) = VAUX
c         VAUX = V2(J)
c         V2(J) = VAUX - THREECST/TWO*WORK(J) - CSGAMMA1/TWO*VAUX 
c     &        - THREECST/TWO*CEP2*U1(J) + THREECST/TWO*ANT*U1(I)
c     &        + THREECST/TWO*OMGT*(U1(I+1) - U1(I-1)) 
c     &        + CST/TWO*HV1(J) 
c         HV1(J) = WORK(J) + GAMMA*V1(J) + HV2(I)
c         V1(J) = VAUX
c      ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c
      MC = ONE
      T = TI + DT
      IF(MC.EQ.LC)MC = ZERO
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &     PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.LC)THEN
c         
         CRA = FCORREL(N, U0(1), U2(1)) 
     &        + FCORREL(N, V0(1), V2(1))
         CIA = FCORREL(N, U0(1), V2(1)) 
     &        - FCORREL(N, V0(1), U2(1))
         FNRA = FCORREL(N, U2(1), U2(1)) 
     &        + FCORREL(N, V2(1), V2(1))
c     FNI must be igual to zero
D         FNIA = FCORREL(N, U2(1), V2(1)) 
D     &        - FCORREL(N, V2(1), U2(1))
         FNTA = FNRA !+ FNIA
         CTA = CRA**2 + CIA**2
c     
         ERA = EIGENERG(DIM, N, NP, SHM, U2(1), VPOT(1), WORK(1), VAR)
         EIA = EIGENERG(DIM, N, NP, SHM, V2(1), VPOT(1), WORK(1), VAR)
         ETA  = (ERA + EIA)
         ER0A = EIGENERG(DIM, N, NP, SHM, U2(1), V(1), WORK(1), VAR)
         EI0A = EIGENERG(DIM, N, NP, SHM, V2(1), V(1), WORK(1), VAR)
         ET0A  = (ER0A + EI0A)  
c
c
         CRB = FCORREL(N, U0(1), U2(NP1)) 
     &        + FCORREL(N, V0(1), V2(NP1))
         CIB = FCORREL(N, U0(1), V2(NP1)) 
     &        - FCORREL(N, V0(1), U2(NP1))
         FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &        + FCORREL(N, V2(NP1), V2(NP1))
c     FNI must be igual to zero
D         FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &        - FCORREL(N, V2(NP1), U2(NP1))
         FNTB = FNRB !+ FNIB
         CTB = CRB**2 + CIB**2
c    
         ERB = EIGENERG(DIM, N, NP, SHM, U2(NP1), VPOT(NP1),
     &        WORK(NP1), VAR)
         EIB = EIGENERG(DIM, N, NP, SHM, V2(NP1), VPOT(NP1), 
     &        WORK(NP1), VAR)
         ETB  = (ERB + EIB)
         ER0B = EIGENERG(DIM, N, NP, SHM, U2(NP1), V(NP1), 
     &        WORK(NP1), VAR)
         EI0B = EIGENERG(DIM, N, NP, SHM, V2(NP1), V(NP1), 
     &        WORK(NP1), VAR)
         ET0B  = (ER0B + EI0B)
c     
         ET = ETA + ETB
         ET0 = ET0A + ET0B
         FNT = FNTA + FNTB
c
         WRITE(*,1021)T, ET, ET0, CRA, CIA, CTA, CTB, FNTA, FNTB, FNT 
c

         IF(MC1.EQ.NSHOT)THEN
            MC1 = ZERO
            MC2 = MC2 + ONE
            WRITE(CHNUM,'(I4)')MC2
            NEWNAM6 = 'ReImB_000'//CHNUM(4:4)//'.dat'
            NEWNAM5 = 'veffB_000'//CHNUM(4:4)//'.dat'
            NEWNAM4 = 'eigvcB_000'//CHNUM(4:4)//'.dat'
            NEWNAM3 = 'ReImA_000'//CHNUM(4:4)//'.dat'
            NEWNAM2 = 'veffA_000'//CHNUM(4:4)//'.dat'
            NEWNAM1 = 'eigvcA_000'//CHNUM(4:4)//'.dat'
            IF(PRTVEFF)THEN
               CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, 
     &              VPOT(1))
               CALL PRTPT(NEWNAM5, 9, ND, T, NP, XP, XI, SH, 
     &              VPOT(NP1))
            ENDIF
            IF(PRTEIGVC2)THEN
               WRITE(22,*)'#set output "',NEWNAM1,'"'
               WRITE(22,1031)NEWNAM1,t,NEWNAM3
               DO I=1,N,1
                  WORK(I) = U2(I)*U2(I) + V2(I)*V2(I) + ET0
               ENDDO
               CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
               CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, 
     &              U2(1), V2(1))
               DO I=NP1,TWON,1
                  WORK(I) = U2(I)*U2(I)+V2(I)*V2(I) + ET0
               ENDDO
               CALL PRTPT(NEWNAM4, 9, ND, T, NP, XP, XI, SH, WORK)
               CALL PRPT2(NEWNAM6, 8, ND, T, NP, XP, XI, SH, 
     &              U2(NP1), V2(NP1))
            ENDIF
         ENDIF
         IF(PRTPULS)THEN
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T,EFAUX       !/SQRT(E0)/2.1132D-9
         ENDIF 
      ENDIF
c      
      DO I=1,2*N,1
         U1(I) = U0(I)
         V1(I) = V0(I)
      ENDDO 
      DO I=1,2,1
         NN = (I - ONE)*N + ONE
         CALL AU(DIM, N, NP, SHM, VPOT(NN), V0(NN), HU1(NN), VAR)
         CALL AU(DIM, N, NP, SHM, VPOT(NN), U0(NN), HV1(NN), VAR)
      ENDDO
      DO I=1,N,1
         J = I + N
         HU1(I) = HU1(I) - GAMMA*U0(I)
     &        + CEP2*V0(I) + ANT*V0(J)
     &        + OMGT*(V0(J+1) - V0(J-1)) 
         HU1(J) = HU1(J) - GAMMA*U0(J)
     &        + CEP2*V0(J) - ANT*V0(I)
     &        - OMGT*(V0(I+1) - V0(I-1))
         HV1(I) = HV1(I) + GAMMA*V0(I) 
     &        - CEP2*U0(I) - ANT*U0(J)
     &        - OMGT*(U0(J+1) - U0(J-1))
         HV1(J) = HV1(J) + GAMMA*V0(J) 
     &        - CEP2*U0(J) + ANT*U0(I)
     &        + OMGT*(U0(I+1) - U0(I-1)) 
      ENDDO

      SUMR = ZERO
      SUMI = ZERO
      FLG = .FALSE.
c      
      DO T=TI+2*DT,TF+DT/2,DT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         DO I=1,2,1
            NN = (I - ONE)*N + ONE
            IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, 
     &           SNI, KL, DM(NN), V(NN), VPOT(NN))
D            CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST, THREECST, 
D     &        CSGAMMA1, CSGAMMA2, NP, SHM, U1(NN), U2(NN), 
D     &        V1(NN), V2(NN), HU1(NN), HV1(NN), VPOT(NN), VABC(NN), 
D     &        WORK)
         ENDDO
c       
         DO I=1,2,1
            NN = (I - ONE)*N + ONE
            CALL AU(DIM, N, NP, SHM, VPOT(NN), V2(NN), WORK(NN), VAR)
         ENDDO
         XD = XI(1)   
         DO I=2,N-1,1
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2*EP*EP
            OMGT = CSQC*EP
            ANT = -CSQC1*(XD - QC)/DELQ**2*EP
            J = I + N
            UAUX  = U2(I)
            U2(I) = UAUX + THREECST*WORK(I) - CSGAMMA1*UAUX
     &           - THREECST*CEP2*V2(I) + THREECST*ANT*V2(J)
     &           + THREECST*OMGT!*(V2(J+1) - V2(J-1))
     &           *(-V1(J+2) + EHT*V1(J+1) - EHT*V1(J-1) + V1(J-2)) 
     &           - CST*HU1(I) 
            HU1(I) = WORK(I) - GAMMA*U1(I)
     &           + CEP2*V1(I) + ANT*V1(J)
     &           + OMGT!*(V1(J+1) - V1(J-1))
     &           *(-V1(J+2) + EHT*V1(J+1) - EHT*V1(J-1) + V1(J-2))
            HV2(I) = - CEP2*U1(I) - ANT*U1(J)
     &           - OMGT!*(U1(J+1) - U1(J-1))
     &           *(-U1(J+2) + EHT*U1(J+1) - EHT*U1(J-1) + U1(J-2))
            U1(I) = UAUX

            UAUX  = U2(J)
            U2(J) = UAUX + THREECST*WORK(J) - CSGAMMA1*UAUX
     &           - THREECST*CEP2*V2(J) - THREECST*ANT*V2(I)
     &           - THREECST*OMGT!*(V2(I+1) - V2(I-1)) 
     &           *(-V1(I+2) + EHT*V1(I+1) - EHT*V1(I-1) + V1(I-2)) 
     &           - CST*HU1(J) 
            HU1(J) = WORK(J) - GAMMA*U1(J)
     &           + CEP2*V1(J) - ANT*V1(I)
     &           - OMGT!*(V1(I+1) - V1(I-1))
     &           *(-V1(I+2) + EHT*V1(I+1) - EHT*V1(I-1) + V1(I-2)) 
            HV2(J) = - CEP2*U1(J) + ANT*U1(I)
     &           + OMGT!*(U1(I+1) - U1(I-1))
     &           *(-U1(I+2) + EHT*U1(I+1) - EHT*U1(I-1) + U1(I-2))
            U1(J) = UAUX
         ENDDO
c
         DO I=1,2,1
            NN = (I - ONE)*N + ONE
            CALL AU(DIM, N, NP, SHM, VPOT(NN), U1(NN), WORK(NN), VAR) 
         ENDDO
         XD = XI(1)
         DO I=2,N-1,1
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2*EP*EP
            OMGT = CSQC*EP
            ANT  = CSQC1*(XD - QC)/DELQ**2*EP
            J = I + N
            VAUX = V2(I)
            V2(I) = VAUX - THREECST*WORK(I) - CSGAMMA1*VAUX
     &           + THREECST*CEP2*U1(I) - THREECST*ANT*U1(J)
     &           - THREECST*OMGT!*(U1(J+1) - U1(J-1)) 
     &           *(-U1(J+2) + EHT*U1(J+1) - EHT*U1(J-1) + U1(J-2))
     &           + CST*HV1(I)
            HV1(I) = WORK(I) + GAMMA*V1(I) + HV2(I)
            V1(I) = VAUX
            VAUX = V2(J)
            V2(J) = VAUX - THREECST*WORK(J) - CSGAMMA1*VAUX 
     &           + THREECST*CEP2*U1(J) + THREECST*ANT*U1(I)
     &           + THREECST*OMGT!*(U1(I+1) - U1(I-1)) 
     &           *(-U1(I+2) + EHT*U1(I+1) - EHT*U1(I-1) + U1(I-2))
     &           + CST*HV1(J) 
            HV1(J) = WORK(J) + GAMMA*V1(J) + HV2(I)
            V1(J) = VAUX
         ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         IF(FLG)THEN
            CRB2A = FCORREL(N, U2(1), U2(NP1)) 
     &           + FCORREL(N, V2(1), V2(NP1))
            CIB2A = FCORREL(N, U2(1), V2(NP1)) 
     &           - FCORREL(N, V2(1), U2(NP1))
            CRB2 = EFXDM*(CRB2A*COMGT - CIB2A*SOMGT)
            CIB2 = EFXDM*(CRB2A*SOMGT + CIB2A*COMGT)
c
            FLG = .FALSE.  
         ELSE   
            CRB3A = FCORREL(N, U2(1), U2(NP1)) 
     &           + FCORREL(N, V2(1), V2(NP1))
            CIB3A = FCORREL(N, U2(1), V2(NP1)) 
     &           - FCORREL(N, V2(1), U2(NP1))
            CRB3 = EFXDM*(CRB3A*COMGT - CIB3A*SOMGT)
            CIB3 = EFXDM*(CRB3A*SOMGT + CIB3A*COMGT)
c
            SAUXR = CRB1 + FOUR*CRB2 + CRB3
            SUMR = SUMR + DTD3*SAUXR
            CRB1 = CRB3
c
            SAUXI = CIB1 + FOUR*CIB2 + CIB3
            SUMI = SUMI + DTD3*SAUXI
            CIB1 = CIB3
c
            FLG = .TRUE.
         ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         MC = MC + ONE
         IF(MC.EQ.LC)MC = ZERO
c
         IF(PRTCRL(1:4).EQ.'.YES' .OR. 
     &        PRTCRL(1:8).EQ.'.PARTIAL' .AND. MC.EQ.ZERO)THEN
c    
            CRA = FCORREL(N, U0(1), U2(1)) 
     &           + FCORREL(N, V0(1), V2(1))
            CIA = FCORREL(N, U0(1), V2(1)) 
     &           - FCORREL(N, V0(1), U2(1))
            FNRA = FCORREL(N, U2(1), U2(1)) 
     &           + FCORREL(N, V2(1), V2(1))
c     FNI must be igual to zero
D            FNIA = FCORREL(N, U2(1), V2(1)) 
D     &           - FCORREL(N, V2(1), U2(1))
            FNTA = FNRA !+ FNIA
            CTA = CRA**2 + CIA**2
c     
            ERA = EIGENERG(DIM, N, NP, SHM, U2(1), VPOT(1), WORK(1), 
     &           VAR)
            EIA = EIGENERG(DIM, N, NP, SHM, V2(1), VPOT(1), WORK(1), 
     &           VAR)
            ETA  = (ERA + EIA)
            ER0A = EIGENERG(DIM, N, NP, SHM, U2(1), V(1), WORK(1), VAR)
            EI0A = EIGENERG(DIM, N, NP, SHM, V2(1), V(1), WORK(1), VAR)
            ET0A  = (ER0A + EI0A)  
c     
            CRB = FCORREL(N, U0(NP1), U2(NP1)) 
     &           + FCORREL(N, V0(NP1), V2(NP1))
            CIB = FCORREL(N, U0(NP1), V2(NP1)) 
     &           - FCORREL(N, V0(NP1), U2(NP1))
            FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &           + FCORREL(N, V2(NP1), V2(NP1))
c     FNI must be igual to zero
D            FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &           - FCORREL(N, V2(NP1), U2(NP1))
            FNTB = FNRB !+ FNIB
            CTB = CRB**2 + CIB**2
c     
            ERB = EIGENERG(DIM, N, NP, SHM, U2(NP1), VPOT(NP1),
     &           WORK(NP1), VAR)
            EIB = EIGENERG(DIM, N, NP, SHM, V2(NP1), VPOT(NP1), 
     &           WORK(NP1), VAR)
            ETB  = (ERB + EIB)
            ER0B = EIGENERG(DIM, N, NP, SHM, U2(NP1), V(NP1), 
     &           WORK(NP1), VAR)
            EI0B = EIGENERG(DIM, N, NP, SHM, V2(NP1), V(NP1), 
     &           WORK(NP1), VAR)
            ET0B  = (ER0B + EI0B)
c     
            ET = ETA + ETB
            ET0 = ET0A + ET0B
            FNT = FNTA + FNTB
c     
            IF(ABSORB .OR. GAMMA.NE.ZERO)THEN
               CONTINUE
            ELSE
               IF(ABS((FNTA-F0A)/F0A)*1.0D+2.GT.TOL)INFO = ONE
            ENDIF
c
            WRITE(*,1021)T, ET, ET0, CRA, CIA, CTA, CTB, FNTA, FNTB, FNT  
D            write(*,*)eta,etb
c            WRITE(*,1021)T, ET, ET0, CR, CI, (CT(I), I=1,NPR+1,1), FNT
c
            IF(MC.EQ.ZERO) MC1 = MC1 + ONE
            IF(MC1.EQ.NSHOT)THEN
               MC1 = ZERO
               MC2 = MC2 + ONE
               WRITE(CHNUM,'(I4)')MC2
               IF(MC2.LT.10)THEN
                  NEWNAM6 = 'ReImB_000'//CHNUM(4:4)//'.dat'
                  NEWNAM5 = 'veffB_000'//CHNUM(4:4)//'.dat'
                  NEWNAM4 = 'eigvcB_000'//CHNUM(4:4)//'.dat'
                  NEWNAM3 = 'ReImA_000'//CHNUM(4:4)//'.dat'
                  NEWNAM2 = 'veffA_000'//CHNUM(4:4)//'.dat'
                  NEWNAM1 = 'eigvcA_000'//CHNUM(4:4)//'.dat'
               ELSEIF(MC2.LT.100)THEN
                  NEWNAM6 = 'ReImB_00'//CHNUM(3:4)//'.dat'
                  NEWNAM5 = 'veffB_00'//CHNUM(3:4)//'.dat'
                  NEWNAM4 = 'eigvcB_00'//CHNUM(3:4)//'.dat'
                  NEWNAM3 = 'ReImA_00'//CHNUM(3:4)//'.dat'
                  NEWNAM2 = 'veffA_00'//CHNUM(3:4)//'.dat'
                  NEWNAM1 = 'eigvcA_00'//CHNUM(3:4)//'.dat'
               ELSEIF(MC2.LT.1000)THEN
                  NEWNAM6 = 'ReImB_0'//CHNUM(2:4)//'.dat'
                  NEWNAM5 = 'veffB_0'//CHNUM(2:4)//'.dat'
                  NEWNAM4 = 'eigvcB_0'//CHNUM(2:4)//'.dat'
                  NEWNAM3 = 'ReImA_0'//CHNUM(2:4)//'.dat'
                  NEWNAM2 = 'veffA_0'//CHNUM(2:4)//'.dat'
                  NEWNAM1 = 'eigvcA_0'//CHNUM(2:4)//'.dat'
               ELSEIF(MC2.LT.10000)THEN
                  NEWNAM6 = 'ReImB_'//CHNUM(1:4)//'.dat'
                  NEWNAM5 = 'veffB_'//CHNUM(1:4)//'.dat'
                  NEWNAM4 = 'eigvcB_'//CHNUM(1:4)//'.dat'
                  NEWNAM3 = 'ReImA_'//CHNUM(1:4)//'.dat'
                  NEWNAM2 = 'veffA_'//CHNUM(1:4)//'.dat'
                  NEWNAM1 = 'eigvcA_'//CHNUM(1:4)//'.dat'
               ELSE
                  WRITE(*,*)'Too many files to be printed! Stopping'
                  STOP
               ENDIF
               IF(PRTVEFF)THEN
                  CALL PRTPT(NEWNAM2, 9, ND, T, NP, XP, XI, SH, 
     &                 VPOT(1))
                  CALL PRTPT(NEWNAM5, 9, ND, T, NP, XP, XI, SH, 
     &                 VPOT(NP1))
               ENDIF
               IF(PRTEIGVC2)THEN
                  WRITE(22,*)'#set output "',NEWNAM1,'"'
                  WRITE(22,1031)NEWNAM1,t,NEWNAM3
                  DO I=1,N,1
                     WORK(I) = U2(I)*U2(I) + V2(I)*V2(I) + ET0A
                  ENDDO
                  CALL PRTPT(NEWNAM1, 9, ND, T, NP, XP, XI, SH, WORK)
                  CALL PRPT2(NEWNAM3, 8, ND, T, NP, XP, XI, SH, 
     &                 U2(1), V2(1))
                  DO I=NP1,TWON,1
                     WORK(I) = U2(I)*U2(I) + V2(I)*V2(I) + ET0B
                  ENDDO
                  CALL PRTPT(NEWNAM4, 9, ND, T, NP, XP, XI, SH, WORK)
                  CALL PRPT2(NEWNAM6, 8, ND, T, NP, XP, XI, SH, 
     &                 U2(NP1), V2(NP1))
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
      write(*,*)'omega, sumr, sumi',omega*shbar, sumr, sumi
c                        
      IF(PRTCRL(1:7).EQ.'.LASTWP')THEN
         NEWNAM4 = 'eigvcB_fwp.dat' 
         NEWNAM6 = 'ReImB_fwp.dat'
         NEWNAM1 = 'eigvcA_fwp.dat'
         NEWNAM3 = 'ReImA_fwp.dat'
         DO I=1,N,1
            WORK(I) = U2(I)*U2(I)+V2(I)*V2(I) 
         ENDDO
         CALL PRTPT(NEWNAM1, 9, ND, T-DT, NP, XP, XI, SH, WORK)
         CALL PRPT2(NEWNAM3, 8, ND, T-DT, NP, XP, XI, SH, 
     &        U2(1), V2(1))
         DO I=1,N,1
            WORK(I) = U2(I)*U2(I)+V2(I)*V2(I) 
         ENDDO
         CALL PRTPT(NEWNAM4, 9, ND, T-DT, NP, XP, XI, SH, WORK)
         CALL PRPT2(NEWNAM6, 8, ND, T-DT, NP, XP, XI, SH, 
     &        U2(NP1), V2(NP1))
c         
      ENDIF
c     ..
 1011 FORMAT(/,8X,'*t/fs*',6X,'*ET/a.u.*',6X,'*E0/a.u.*',6X,'*CIA*',
     &     6X,'*CRA*',6X,'*CTA*',6X,'*CTB*',6X,'*FTA*',6X,'*FTB*',6X,
     &     '*FT*')
c 1021 FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),3(F8.6,3X),1(E10.4,3X),
c     &     3(F8.6,3X))
 1021 FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),4(E10.4,3X),F12.7)
 1031 FORMAT('plot "',A15,'" title `',F9.4,1X,
     &     'fs` with lines lw 2.0, ','"',A14,
     &     '" notitle with lines lw 2.1, ',
     &     '"veffA_0001.dat" notitle with lines lt 3 lw 2.2')
      RETURN
      END
