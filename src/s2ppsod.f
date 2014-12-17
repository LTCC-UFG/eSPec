      SUBROUTINE S2PPSOD(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, EFC, 
     &     PRTCRL, TPABSOR, INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, 
     &     NPR, KP, NISPG, KL, E0, T0, TD, TP, OMG, SNI, DT, 
     &     TI, TF, TOL, DELQ, QC, AA, NP, ND, SH, SHM, 
     &     U1, U2, V1, V2, XI, XP, XF, SM, DM, VPOT, VABC, WORK, 
     &     VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL, TPABSOR
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR, KP
      INTEGER       NISPG, KL
      REAL*8        DT, TI, TF, TOL, E0, T0, TD, TP, OMG, SNI
      REAL*8        DELQ, QC, AA
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
*     (20/08/2005) First version written by Freddy based on PSOD
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO, THREE, FOUR, TWOPI, SHBAR, EHT, PI
      PARAMETER(
     &     ZERO  = +0.0D+0, 
     &     ONE   = +1.0D+0, 
     &     TWO   = +2.0D+0,
     &     THREE = +3.0D+0, 
     &     FOUR  = +4.0D+0,
     &     EHT   = +8.0D+0,
     &     PI    = +3.14159265358979323846D+0,
     &     TWOPI = +6.283185307179586476925287D+0,
     &     SHBAR = +6.58211889D-1 !ev*fs
     &     )
c     **
c     ** Local scalars 
C      LOGICAL       DERIV
      CHARACTER*4   CHNUM   
      CHARACTER*15  NEWNAM2, NEWNAM3, NEWNAM5, NEWNAM6
      CHARACTER*16  NEWNAM1, NEWNAM4, TPCRS
      INTEGER       I, J, MC, MC1, MC2, NP1, TWON
      REAL*8        CST, EINI, ET, ET0, FNT, SHFS, T, EFAUX
      REAL*8        CRA, CIA, CTA, ERA, EIA, ETA, ER0A, EI0A, ET0A
      REAL*8        FNRA, FNTA, F0
      REAL*8        CRB, CIB, CTB, ERB, EIB, ETB, ER0B, EI0B, ET0B
      REAL*8        FNRB, FNTB, CSQC, XD, OMGT, ANT, CSQC1, CSQC2
      REAL*8        EP, CEP2
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*3
cdel      INTEGER       
      REAL*8        U0(2*N), V0(2*N), V(2*N)
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
D      DERIV = .FALSE.
D      DERIV = .TRUE.
d      TPCRS = '.WFO'
       TPCRS = '.CORD'
d      TPCRS = '.DERIV'
d      TPCRS = '.PDERIV2'
c
      INFO  = ZERO
      T     = TI
      SHFS  = 2.4188843243D-2   !fs/a.u.(t) [hbar/E_{h}]
      CST   = TWO*DT/SHFS
      MC1   = ZERO
      J     = ZERO
      CSQC1 = ZERO
      CSQC2 = ZERO
      EP    = ZERO 
      CEP2  = ZERO
c      write(*,*)ep,cep2
c
      NP1   = N + ONE
      TWON  = TWO*N
      IF(TPCRS(1:6).EQ.'.DERIV')THEN
         CSQC = CST*AA/(2.0D0*SH(1))/SM(1)!*0.0d0
         CSQC1 = -TWO*CST*AA/SM(1)
      ELSEIF(TPCRS(1:8).EQ.'.PDERIV2')THEN
         CSQC  = CST*AA/(1.2D1*SH(1))/SM(1)!*0d0!*ONE/(DELQ*SQRT(PI))
         CSQC1 = CST*AA/SM(1)!*0d0!*ONE/(DELQ*SQRT(PI))!*0d0
         CSQC2 = CST*AA*AA/(TWO*SM(1))!*0d0!*ONE/(DELQ*SQRT(PI))
c     &        *ONE/(DELQ*SQRT(PI))
      ELSEIF(TPCRS(1:7).EQ.'.PDERIV')THEN
         CSQC  = CST*AA/(TWO*SH(1))/SM(1)!*ONE/(DELQ*SQRT(PI))
         CSQC1 = CST*AA/SM(1)!*ONE/(DELQ*SQRT(PI))
         CSQC2 = CST*AA*AA/(TWO*SM(1))!*ONE/(DELQ*SQRT(PI))
c     &        *ONE/(DELQ*SQRT(PI)) !*0d0
      ELSE
         CSQC = CST*AA/SM(1)
      ENDIF
      write(*,*)'oi',CSQC,CSQC1,CSQC2
      write(*,*)DELQ, QC, AA
c      read(*,*)
      write(*,*)sh(1),(xf(1)-xi(1))/(n-1),shm(1),csqc,shm(1)*30,AA
c
      DO I=1,TWON,1
         V(I)  = VPOT(I)
         U0(I) = U1(I)
         U2(I) = U0(I) 
         V0(I) = V1(I)
         V2(I) = V0(I) 
      ENDDO
c 
      FNRA = FCORREL(N, U2(1), U2(1)) 
     &     + FCORREL(N, V2(1), V2(1))
c     FNIA must be igual to zero
D      FNIA = FCORREL(N, U2(1), V2(1)) 
D     &     - FCORREL(N, V2(1), U2(1))
      FNTA = FNRA
c
      FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &     + FCORREL(N, V2(NP1), V2(NP1))
c     FNIB must be igual to zero
D      FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &     - FCORREL(N, V2(NP1), U2(NP1))
      FNTB = FNRB
      F0 = FNTA + FNTB
c     
      IF(PRTCRL(1:4).EQ.'.YES' .OR. PRTCRL(1:8).EQ.'.PARTIAL')THEN
         CRA = FCORREL(N, U0(1), U2(1)) 
     &        + FCORREL(N, V0(1), V2(1))
         CIA = FCORREL(N, U0(1), V2(1)) 
     &        - FCORREL(N, V0(1), U2(1))
         FNRA = FCORREL(N, U2(1), U2(1)) 
     &        + FCORREL(N, V2(1), V2(1))
c     FNIA must be igual to zero
c         FNIA = FCORREL(N, U2(1), V2(1)) 
c     &        - FCORREL(N, V2(1), U2(1))
         FNTA = FNRA
         CTA = CRA**2 + CIA**2
c
         ERA = EIGENERG(DIM, N, NP, SHM, U2(1), VPOT(1), WORK(1), VAR)
         EIA = EIGENERG(DIM, N, NP, SHM, V2(1), VPOT(1), WORK(1), VAR)
         ETA  = (ERA + EIA)
         ER0A = EIGENERG(DIM, N, NP, SHM, U2(1), V(1), WORK(1), VAR)
         EI0A = EIGENERG(DIM, N, NP, SHM, V2(1), V(1), WORK(1), VAR)
         ET0A  = (ER0A + EI0A)  
c
         CRB = FCORREL(N, U0(NP1), U2(NP1)) 
     &        + FCORREL(N, V0(NP1), V2(NP1))
         CIB = FCORREL(N, U0(NP1), V2(NP1)) 
     &        - FCORREL(N, V0(NP1), U2(NP1))
         FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &        + FCORREL(N, V2(NP1), V2(NP1))
c     FNIB must be igual to zero
c         FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
c     &        - FCORREL(N, V2(NP1), U2(NP1))
         FNTB = FNRB**2     
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
         ELSE
            WRITE(*,*)'Partial auto-correlation function:'
         ENDIF
         WRITE(*,*)'=================================='
         WRITE(*,1011)
         WRITE(*,1021)T, ET, ET0, CRA, CIA, CTA, CRB, CIB, CTB, FNTA, 
     &        FNTB, FNT 
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
      IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, (T+DT)/TWO, OMG, 
     &     SNI, KL, DM(1), V(1), VPOT(1))
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/FOUR, NP, SHM, 
     &     U1(1), U2(1), V1(1), V2(1), 
     &     VPOT(1), VABC(1), WORK, VAR)
      IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, (T+DT)/TWO, OMG, 
     &     SNI, KL, DM(NP1), V(NP1), VPOT(NP1))
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/FOUR, NP, SHM, 
     &     U1(NP1), U2(NP1), V1(NP1), V2(NP1), 
     &     VPOT(NP1), VABC(1), WORK, VAR)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(TPCRS(1:6).EQ.'.DERIV')THEN
c     i=1
         XD = XI(1)
         OMGT = CSQC/FOUR*EXP(-((XD - QC)/DELQ)**2)
         ANT = (XD - QC)/DELQ**2/FOUR*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &        - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
         U2(1) = U2(1) + OMGT*V1(NP1+1)
         V2(1) = V2(1) - OMGT*U1(NP1+1)
         U2(NP1) = U2(NP1) - OMGT*V1(2)
     &        - ANT*V1(1)
         V2(NP1) = V2(NP1) + OMGT*U1(2)
     &        + ANT*U1(1)
         write(2,*)xd,OMGT,u1(1),(u1(3) + u1(2))/(6.0D0*SH(1))
c     i=2:n-1
         DO I=3,N-1,1
            J = I + N
            XD = XD + SH(1)
            OMGT = CSQC/FOUR*EXP(-((XD - QC)/DELQ)**2)
            ANT = (XD - QC)/DELQ**2/FOUR*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &           - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
            U2(I) = U2(I) + OMGT*(V1(J+1) - V1(J-1))
            V2(I) = V2(I) - OMGT*(U1(J+1) - U1(J-1))
            U2(J) = U2(J) - OMGT*(V1(I+1) - V1(I-1))
     &           - ANT*V1(I)
            V2(J) = V2(J) + OMGT*(U1(I+1) - U1(I-1))
     &           + ANT*U1(I)
            write(2,*)xd,OMGT,u1(i),(u1(i+2) + u1(i+1) - u1(i-1) 
     &           - u1(i-2))/(6.0D0*SH(1))
         ENDDO
c     i=n
         XD = XD + SH(1) 
         OMGT = CSQC/FOUR*EXP(-((XD - QC)/DELQ)**2)  
         ANT = (XD - QC)/DELQ**2/FOUR*CSQC1*EXP(-((XD - QC)/DELQ)**2)
     &        - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
         U2(N) = U2(N) - OMGT*V1(TWON-1)
         V2(N) = V2(N) + OMGT*U1(TWON-1)
         U2(TWON) = U2(TWON) + OMGT*V1(N-1)
     &        - ANT*V1(N)
         V2(TWON) = V2(TWON) - OMGT*U1(N-1)
     &        + ANT*U1(N)
         write(2,*)xd,OMGT,u1(n),(- u1(n-1) - u1(n-2))/(6.0D0*SH(1))
         write(2,*)
      ELSEIF(TPCRS(1:8).EQ.'.PDERIV2')THEN
c     i=1
         XD = XI(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/FOUR*EP*EP
         OMGT = CSQC/FOUR*EP
         ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
         U2(1) = U2(1) + CEP2*V1(1) + ANT*V1(NP1) 
     &        + OMGT*(-V1(NP1+2) + EHT*V1(NP1+1)) 
         V2(1) = V2(1) - CEP2*U1(1) - ANT*U1(NP1)
     &        - OMGT*(-U1(NP1+2) + EHT*U1(NP1+1)) 
         U2(NP1) = U2(NP1) + CEP2*V1(NP1) - ANT*V1(1)
     &        - OMGT*(-V1(3) + EHT*V1(2)) 
         V2(NP1) = V2(NP1) - CEP2*U1(NP1) + ANT*U1(1)
     &        + OMGT*(-U1(3) + EHT*U1(2)) 
c     i=2         
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/FOUR*EP*EP
         OMGT = CSQC/FOUR*EP
         ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
         U2(2) = U2(2) + CEP2*V1(2) + ANT*V1(J)
     &        + OMGT*(-V1(NP1+3) + EHT*V1(NP1+2) - EHT*V1(NP1))
         V2(2) = V2(2) - CEP2*U1(2) - ANT*U1(J)
     &        - OMGT*(-U1(NP1+3) + EHT*U1(NP1+2) - EHT*U1(NP1))
         U2(NP1+1) = U2(NP1+1) + CEP2*V1(NP1+1) - ANT*V1(I)
     &        - OMGT*(-V1(4) + EHT*V1(3) - EHT*V1(1))
         V2(NP1+1) = V2(NP1+1) - CEP2*U1(NP1+1) + ANT*U1(I)
     &        + OMGT*(-U1(4) + EHT*U1(3) - EHT*U1(1))
c     i=3:n-2
         DO I=3,N-2,1
            J = I + N
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2/FOUR*EP*EP
            OMGT = CSQC/FOUR*EP
            ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
            U2(I) = U2(I) + CEP2*V1(I) + ANT*V1(J)
     &        + OMGT*(-V1(J+2) + EHT*V1(J+1) - EHT*V1(J-1) + V1(J-2)) 
            V2(I) = V2(I) - CEP2*U1(I) - ANT*U1(J)
     &        - OMGT*(-U1(J+2) + EHT*U1(J+1) - EHT*U1(J-1) + U1(J-2)) 
            U2(J) = U2(J) + CEP2*V1(J) - ANT*V1(I)
     &        - OMGT*(-V1(I+2) + EHT*V1(I+1) - EHT*V1(I-1) + V1(I-2)) 
            V2(J) = V2(J) - CEP2*U1(J) + ANT*U1(I)
     &        + OMGT*(-U1(I+2) + EHT*U1(I+1) - EHT*U1(I-1) + U1(I-2)) 
         ENDDO
c     i=N-1 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/FOUR*EP*EP
         OMGT = CSQC/FOUR*EP
         ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
         U2(N-1) = U2(N-1) + CEP2*V1(N-1) + ANT*V1(J)
     &        + OMGT*(EHT*V1(TWON) - EHT*V1(TWON-2) + V1(TWON-3))
         V2(N-1) = V2(N-1) - CEP2*U1(N-1) - ANT*U1(J)
     &        - OMGT*(EHT*U1(TWON) - EHT*U1(TWON-2) + U1(TWON-3))
         U2(TWON-1) = U2(TWON-1) + CEP2*V1(TWON-1) - ANT*V1(I)
     &        - OMGT*(EHT*V1(N) - EHT*V1(N-2) + V1(N-3))
         V2(TWON-1) = V2(TWON-1) - CEP2*U1(TWON-1) + ANT*U1(I)
     &        + OMGT*(EHT*U1(N) - EHT*U1(N-2) + U1(N-3))
c     i=N 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/FOUR*EP*EP
         OMGT = CSQC/FOUR*EP
         ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
         U2(N) = U2(N) + CEP2*V1(N) + ANT*V1(J)
     &        + OMGT*(-EHT*V1(TWON-1) + V1(TWON-2))
         V2(N) = V2(N) - CEP2*U1(N) - ANT*U1(J)
     &        - OMGT*(-EHT*U1(TWON-1) + U1(TWON-2))
         U2(TWON) = U2(TWON) + CEP2*V1(TWON) - ANT*V1(I)
     &        - OMGT*(-EHT*V1(N-1) + V1(N-2))
         V2(TWON) = V2(TWON) - CEP2*U1(TWON) + ANT*U1(I)
     &        + OMGT*(-EHT*U1(N-1) + U1(N-2))
      ELSEIF(TPCRS(1:7).EQ.'.PDERIV')THEN
c     i=1
         XD = XI(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/FOUR*EP*EP
         OMGT = CSQC/FOUR*EP
         ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
         U2(1) = U2(1) + CEP2*V1(1) 
     &        + OMGT*V1(NP1+1) + ANT*V1(NP1)
         V2(1) = V2(1) - CEP2*U1(1) 
     &        - OMGT*(U1(NP1+2) + U1(NP1+1)) - ANT*U1(NP1)
         U2(NP1) = U2(NP1) + CEP2*V1(NP1)
     &        - OMGT*(V1(3) + V1(2)) - ANT*V1(1)
         V2(NP1) = V2(NP1) - CEP2*U1(NP1)
     &        + OMGT*(U1(3) + U1(2)) + ANT*U1(1)
c     i=2:n-1
         DO I=2,N-1,1
            J = I + N
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2/FOUR*EP*EP
            OMGT = CSQC/FOUR*EP
            ANT = -CSQC1/FOUR*(XD - QC)/DELQ**2*EP
            U2(I) = U2(I) + CEP2*V1(I)
     &        + OMGT*(V1(J+1) - V1(J-1)) + ANT*V1(J)
            V2(I) = V2(I) - CEP2*U1(I)
     &        - OMGT*(U1(J+1) - U1(J-1)) - ANT*U1(J)
            U2(J) = U2(J) + CEP2*V1(J)
     &        - OMGT*(V1(I+1) - V1(I-1)) - ANT*V1(I)
            V2(J) = V2(J) - CEP2*U1(J)
     &        + OMGT*(U1(I+1) - U1(I-1)) + ANT*U1(I)
         ENDDO
c     i=n
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/FOUR*EP*EP
         OMGT = CSQC/FOUR*EP
         ANT = CSQC1/FOUR*(XD - QC)/DELQ**2*EP
         U2(N) = U2(N) + CEP2*V1(N)
     &        - OMGT*V1(TWON-1) + ANT*V1(J)
         V2(N) = V2(N) - CEP2*U1(N)
     &        + OMGT*U1(TWON-1) - ANT*U1(J)
         U2(TWON) = U2(TWON) + CEP2*V1(TWON)
     &        + OMGT*V1(N-1) - ANT*V1(I)
         V2(TWON) = V2(TWON) - CEP2*V1(TWON)
     &        - OMGT*U1(N-1) + ANT*U1(I)
      ELSEIF(TPCRS(1:5).EQ.'.CORD')THEN
         XD = XI(1)
         DO I=1,N,1
            J = N + I
            XD = XD + SH(1)
            OMGT = CSQC/FOUR*XD
D            OMGT = CSQC/FOUR*EXP(-((XD - QC)/DELQ)**2)*XD
            U2(I) = U2(I) + OMGT*V1(J)
            V2(I) = V2(I) - OMGT*U1(J)
            U2(J) = U2(J) + OMGT*V1(I)
            V2(J) = V2(J) - OMGT*U1(I)
         ENDDO
      ELSE
         XD = XI(1)
         DO I=1,N,1
            J = N + I
            XD = XD + SH(1)
            OMGT = CSQC/FOUR*EXP(-((XD - QC)/DELQ)**2)
            U2(I) = U2(I) + OMGT*V1(J)
            V2(I) = V2(I) - OMGT*U1(J)
            U2(J) = U2(J) + OMGT*V1(I)
            V2(J) = V2(J) - OMGT*U1(I)
         ENDDO
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     .. Second half step in the propagation towards the point T+DT (SOD) 
      IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, SNI, 
     &     KL, DM(1), V(1), VPOT(1))
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, NP, SHM, 
     &     U1(1), U2(1), V1(1), V2(1), VPOT(1), VABC(1), WORK, VAR)
      IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, SNI, 
     &     KL, DM(NP1), V(NP1), VPOT(NP1))
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, NP, SHM, 
     &     U1(NP1), U2(NP1), V1(NP1), V2(NP1),  
     &     VPOT(NP1), VABC(1), WORK, VAR)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(TPCRS(1:6).EQ.'.DERIV')THEN
c     i=1
         XD = XI(1)
         OMGT = CSQC/TWO*EXP(-((XD - QC)/DELQ)**2)
         ANT = (XD - QC)/DELQ**2/TWO*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &        - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
         U2(1) = U2(1) + OMGT*V1(NP1+1)
         V2(1) = V2(1) - OMGT*U1(NP1+1)
         U2(NP1) = U2(NP1) - OMGT*V1(2)
     &        - ANT*V1(1)
         V2(NP1) = V2(NP1) + OMGT*U1(2)
     &        + ANT*U1(1)
c     i=2:n-1
         DO I=2,N-1,1
            J = I + N
            XD = XD + SH(1)
            OMGT = CSQC/TWO*EXP(-((XD - QC)/DELQ)**2)
            ANT = (XD - QC)/DELQ**2/TWO*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &           - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
C            write(*,*)ant,omgt!,csqc1
            U2(I) = U2(I) + OMGT*(V1(J+1) - V1(J-1))
            V2(I) = V2(I) - OMGT*(U1(J+1) - U1(J-1))
            U2(J) = U2(J) - OMGT*(V1(I+1) - V1(I-1))
     &        - ANT*V1(I)
            V2(J) = V2(J) + OMGT*(U1(I+1) - U1(I-1))
     &        + ANT*U1(I) 
            write(2,*)xd,OMGT,u1(i),(u1(i+2) + u1(i+1) - u1(i-1) 
     &           - u1(i-2))/(6.0D0*SH(1))
         ENDDO
         write(2,*)
c     i=n
         XD = XD + SH(1) 
         OMGT = CSQC/TWO*EXP(-((XD - QC)/DELQ)**2)  
         ANT = (XD - QC)/DELQ**2/TWO*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &        - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
         U2(N) = U2(N) - OMGT*V1(TWON-1)
         V2(N) = V2(N) + OMGT*U1(TWON-1)
         U2(TWON) = U2(TWON) + OMGT*V1(N-1)
     &        - ANT*V1(N)
         V2(TWON) = V2(TWON) - OMGT*U1(N-1)
     &        + ANT*U1(N)
      ELSEIF(TPCRS(1:8).EQ.'.PDERIV2')THEN
c     i=1
         XD = XI(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/TWO*EP*EP
         OMGT = CSQC/TWO*EP
         ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
         U2(1) = U2(1) + CEP2*V1(1) + ANT*V1(NP1) 
     &        + OMGT*(-V1(NP1+2) + EHT*V1(NP1+1)) 
         V2(1) = V2(1) - CEP2*U1(1) - ANT*U1(NP1)
     &        - OMGT*(-U1(NP1+2) + EHT*U1(NP1+1)) 
         U2(NP1) = U2(NP1) + CEP2*V1(NP1) - ANT*V1(1)
     &        - OMGT*(-V1(3) + EHT*V1(2)) 
         V2(NP1) = V2(NP1) - CEP2*U1(NP1) + ANT*U1(1)
     &        + OMGT*(-U1(3) + EHT*U1(2)) 
c     i=2         
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/TWO*EP*EP
         OMGT = CSQC/TWO*EP
         ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
         U2(2) = U2(2) + CEP2*V1(2) + ANT*V1(J)
     &        + OMGT*(-V1(NP1+3) + EHT*V1(NP1+2) - EHT*V1(NP1))
         V2(2) = V2(2) - CEP2*U1(2) - ANT*U1(J)
     &        - OMGT*(-U1(NP1+3) + EHT*U1(NP1+2) - EHT*U1(NP1))
         U2(NP1+1) = U2(NP1+1) + CEP2*V1(NP1+1) - ANT*V1(I)
     &        - OMGT*(-V1(4) + EHT*V1(3) - EHT*V1(1))
         V2(NP1+1) = V2(NP1+1) - CEP2*U1(NP1+1) + ANT*U1(I)
     &        + OMGT*(-U1(4) + EHT*U1(3) - EHT*U1(1))
c     i=3:n-2
         DO I=3,N-2,1
            J = I + N
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2/TWO*EP*EP
            OMGT = CSQC/TWO*EP
            ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
            U2(I) = U2(I) + CEP2*V1(I) + ANT*V1(J)
     &        + OMGT*(-V1(J+2) + EHT*V1(J+1) - EHT*V1(J-1) + V1(J-2)) 
            V2(I) = V2(I) - CEP2*U1(I) - ANT*U1(J)
     &        - OMGT*(-U1(J+2) + EHT*U1(J+1) - EHT*U1(J-1) + U1(J-2)) 
            U2(J) = U2(J) + CEP2*V1(J) - ANT*V1(I)
     &        - OMGT*(-V1(I+2) + EHT*V1(I+1) - EHT*V1(I-1) + V1(I-2)) 
            V2(J) = V2(J) - CEP2*U1(J) + ANT*U1(I)
     &        + OMGT*(-U1(I+2) + EHT*U1(I+1) - EHT*U1(I-1) + U1(I-2)) 
         ENDDO
c     i=N-1 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/TWO*EP*EP
         OMGT = CSQC/TWO*EP
         ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
         U2(N-1) = U2(N-1) + CEP2*V1(N-1) + ANT*V1(J)
     &        + OMGT*(EHT*V1(TWON) - EHT*V1(TWON-2) + V1(TWON-3))
         V2(N-1) = V2(N-1) - CEP2*U1(N-1) - ANT*U1(J)
     &        - OMGT*(EHT*U1(TWON) - EHT*U1(TWON-2) + U1(TWON-3))
         U2(TWON-1) = U2(TWON-1) + CEP2*V1(TWON-1) - ANT*V1(I)
     &        - OMGT*(EHT*V1(N) - EHT*V1(N-2) + V1(N-3))
         V2(TWON-1) = V2(TWON-1) - CEP2*U1(TWON-1) + ANT*U1(I)
     &        + OMGT*(EHT*U1(N) - EHT*U1(N-2) + U1(N-3))
c     i=N 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/TWO*EP*EP
         OMGT = CSQC/TWO*EP
         ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
         U2(N) = U2(N) + CEP2*V1(N) + ANT*V1(J)
     &        + OMGT*(-EHT*V1(TWON-1) + V1(TWON-2))
         V2(N) = V2(N) - CEP2*U1(N) - ANT*U1(J)
     &        - OMGT*(-EHT*U1(TWON-1) + U1(TWON-2))
         U2(TWON) = U2(TWON) + CEP2*V1(TWON) - ANT*V1(I)
     &        - OMGT*(-EHT*V1(N-1) + V1(N-2))
         V2(TWON) = V2(TWON) - CEP2*U1(TWON) + ANT*U1(I)
     &        + OMGT*(-EHT*U1(N-1) + U1(N-2))
      ELSEIF(TPCRS(1:7).EQ.'.PDERIV')THEN
c     i=1
         XD = XI(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/TWO*EP*EP
         OMGT = CSQC/TWO*EP
         ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
         U2(1) = U2(1) + CEP2*V1(1) 
     &        + OMGT*V1(NP1+1) + ANT*V1(NP1)
         V2(1) = V2(1) - CEP2*U1(1) 
     &        - OMGT*U1(NP1+1) - ANT*U1(NP1)
         U2(NP1) = U2(NP1) + CEP2*V1(NP1)
     &        - OMGT*V1(2) - ANT*V1(1)
         V2(NP1) = V2(NP1) - CEP2*U1(NP1)
     &        + OMGT*U1(2) + ANT*U1(1)
c     i=2:n-1
         DO I=2,N-1,1
            J = I + N
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2/TWO*EP*EP
            OMGT = CSQC/TWO*EP
            ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
            U2(I) = U2(I) + CEP2*V1(I)
     &        + OMGT*(V1(J+1) - V1(J-1)) + ANT*V1(J)
            V2(I) = V2(I) - CEP2*U1(I)
     &        - OMGT*(U1(J+1) - U1(J-1)) - ANT*U1(J)
            U2(J) = U2(J) + CEP2*V1(J)
     &        - OMGT*(V1(I+1) - V1(I-1)) - ANT*V1(I)
            V2(J) = V2(J) - CEP2*U1(J)
     &        + OMGT*(U1(I+1) - U1(I-1)) + ANT*U1(I)
         ENDDO
c     i=n
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2/TWO*EP*EP
         OMGT = CSQC/TWO*EP
         ANT = -CSQC1/TWO*(XD - QC)/DELQ**2*EP
         U2(N) = U2(N) + CEP2*V1(N)
     &        - OMGT*V1(TWON-1) + ANT*V1(J)
         V2(N) = V2(N) - CEP2*U1(N)
     &        + OMGT*U1(TWON-1) - ANT*U1(J)
         U2(TWON) = U2(TWON) + CEP2*V1(TWON)
     &        + OMGT*V1(N-1) - ANT*V1(I)
         V2(TWON) = V2(TWON) - CEP2*U1(TWON)
     &        - OMGT*U1(N-1) + ANT*U1(I)
      ELSEIF(TPCRS(1:5).EQ.'.CORD')THEN
         XD = XI(1)
         DO I=1,N,1
            J = N + I
            XD = XD + SH(1)
            OMGT = CSQC/TWO*XD
D            OMGT = CSQC/TWO*EXP(-((XD - QC)/DELQ)**2)*XD
            U2(I) = U2(I) + OMGT*V1(J)
            V2(I) = V2(I) - OMGT*U1(J)
            U2(J) = U2(J) + OMGT*V1(I)
            V2(J) = V2(J) - OMGT*U1(I)
         ENDDO
      ELSE
         XD = XI(1)
         DO I=1,N,1
            J = N + I
            XD = XD + SH(1)
            OMGT = CSQC/TWO*EXP(-((XD - QC)/DELQ)**2)
            U2(I) = U2(I) + OMGT*V1(J)
            V2(I) = V2(I) - OMGT*U1(J)
            U2(J) = U2(J) + OMGT*V1(I)
            V2(J) = V2(J) - OMGT*U1(I)
         ENDDO
      ENDIF
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
c     FNIA must be igual to zero
D     FNIA = FCORREL(N, U2(1), V2(1)) 
D     &        - FCORREL(N, V2(1), U2(1))
         FNTA = FNRA
         CTA = CRA**2 + CIA**2
c
         ERA = EIGENERG(DIM, N, NP, SHM, U2(1), VPOT(1), WORK(1), VAR)
         EIA = EIGENERG(DIM, N, NP, SHM, V2(1), VPOT(1), WORK(1), VAR)
         ETA  = (ERA + EIA)
         ER0A = EIGENERG(DIM, N, NP, SHM, U2(1), V(1), WORK(1), VAR)
         EI0A = EIGENERG(DIM, N, NP, SHM, V2(1), V(1), WORK(1), VAR)
         ET0A  = (ER0A + EI0A)  
c
         CRB = FCORREL(N, U0(NP1), U2(NP1)) 
     &        + FCORREL(N, V0(NP1), V2(NP1))
         CIB = FCORREL(N, U0(NP1), V2(NP1)) 
     &        - FCORREL(N, V0(NP1), U2(NP1))
         FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &        + FCORREL(N, V2(NP1), V2(NP1))
c     FNIB must be igual to zero
D     FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &        - FCORREL(N, V2(NP1), U2(NP1))
         FNTB = FNRB
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
         WRITE(*,1021)T, ET, ET0, CRA, CIA, CTA, CTB, CIB, CTB, FNTA, 
     &        FNTB, FNT 
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
      DO I=1,TWON,1
         U1(I) = U0(I)
         V1(I) = V0(I)
D         write(*,*)u2(i),v2(i),u1(i),u0(i),v1(i),v0(i)
D         read(*,*)        
      ENDDO
c
      DO T=TI+2*DT,TF+DT/2,DT
c     .. Steps in the propagation towards T (SOD).
         IF(E0.NE.ZERO)THEN    
            CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, 
     &           DM(1), V(1), VPOT(1))
         ENDIF 
         CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST, NP, SHM, 
     &        U1(1), U2(1), V1(1), V2(1), 
     &        VPOT(1), VABC(1), WORK, VAR)
c
         IF(E0.NE.ZERO)THEN
            CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, 
     &           DM(NP1), V(NP1), VPOT(NP1))
         ENDIF
         CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST, NP, SHM, 
     &        U1(NP1), U2(NP1), V1(NP1), V2(NP1),  
     &        VPOT(NP1), VABC(1), WORK, VAR)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(TPCRS(1:6).EQ.'.DERIV')THEN
c     i=1
         XD = XI(1)
         OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)
         ANT = (XD - QC)/DELQ**2*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &        - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
         U2(1) = U2(1) + OMGT*V1(NP1+1)
         V2(1) = V2(1) - OMGT*U1(NP1+1)
         U2(NP1) = U2(NP1) - OMGT*V1(2)
     &        - ANT*V1(1)
         V2(NP1) = V2(NP1) + OMGT*U1(2)
     &        + ANT*U1(1)
c     i=2:n-2
         DO I=2,N-1,1
            J = I + N
            XD = XD + SH(1)
            OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)
            ANT = (XD - QC)/DELQ**2*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &           - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
d            write(*,*)ant,omgt!,csqc1
            U2(I) = U2(I) + OMGT*(V1(J+1) - V1(J-1))
            V2(I) = V2(I) - OMGT*(U1(J+1) - U1(J-1))
            U2(J) = U2(J) - OMGT*(V1(I+1) - V1(I-1))
     &        - ANT*V1(I)
            V2(J) = V2(J) + OMGT*(U1(I+1) - U1(I-1))
     &        + ANT*U1(I)
         ENDDO
c     i=n
         XD = XD + SH(1) 
         OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)  
         ANT = (XD - QC)/DELQ**2*CSQC1*EXP(-((XD - QC)/DELQ)**2)
d     &        - AA*OMGT*EXP(-((XD - QC)/DELQ)**2)
         U2(N) = U2(N) - OMGT*V1(TWON-1)
         V2(N) = V2(N) + OMGT*U1(TWON-1)
         U2(TWON) = U2(TWON) + OMGT*V1(N-1)
     &        - ANT*V1(N)
         V2(TWON) = V2(TWON) - OMGT*U1(N-1)
     &        + ANT*U1(N)
      ELSEIF(TPCRS(1:8).EQ.'.PDERIV2')THEN
c     i=1
         XD = XI(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = -CSQC1*(XD - QC)/DELQ**2*EP
         U2(1) = U2(1) + CEP2*V1(1) + ANT*V1(NP1) 
     &        + OMGT*(-V1(NP1+2) + EHT*V1(NP1+1)) 
         V2(1) = V2(1) - CEP2*U1(1) - ANT*U1(NP1)
     &        - OMGT*(-U1(NP1+2) + EHT*U1(NP1+1)) 
         U2(NP1) = U2(NP1) + CEP2*V1(NP1) - ANT*V1(1)
     &        - OMGT*(-V1(3) + EHT*V1(2)) 
         V2(NP1) = V2(NP1) - CEP2*U1(NP1) + ANT*U1(1)
     &        + OMGT*(-U1(3) + EHT*U1(2)) 
c     i=2         
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = -CSQC1*(XD - QC)/DELQ**2*EP
         U2(2) = U2(2) + CEP2*V1(2) + ANT*V1(J)
     &        + OMGT*(-V1(NP1+3) + EHT*V1(NP1+2) - EHT*V1(NP1))
         V2(2) = V2(2) - CEP2*U1(2) - ANT*U1(J)
     &        - OMGT*(-U1(NP1+3) + EHT*U1(NP1+2) - EHT*U1(NP1))
         U2(NP1+1) = U2(NP1+1) + CEP2*V1(NP1+1) - ANT*V1(I)
     &        - OMGT*(-V1(4) + EHT*V1(3) - EHT*V1(1))
         V2(NP1+1) = V2(NP1+1) - CEP2*U1(NP1+1) + ANT*U1(I)
     &        + OMGT*(-U1(4) + EHT*U1(3) - EHT*U1(1))
c     i=3:n-2
         DO I=3,N-2,1
            J = I + N
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2*EP*EP
            OMGT = CSQC*EP
d            ANT  = CSQC*(EXP(-((XD + SH(1) - QC)/DELQ)**2) 
d     &           - EXP(-((XD - SH(1) - QC)/DELQ)**2))/TWO 
c            ANT  = CSQC*(
c     &           - EXP(-((XD + 2*SH(1) - QC)/DELQ)**2) 
c     &           + EHT*EXP(-((XD + SH(1) - QC)/DELQ)**2)
c     &           - EHT*EXP(-((XD - SH(1) - QC)/DELQ)**2)
c     &           + EXP(-((XD - 2*SH(1) - QC)/DELQ)**2)
c     &           )/TWO
            ANT = -CSQC1*(XD - QC)/DELQ**2*EP
c            write(*,*)ANT,-CSQC1*(XD - QC)/DELQ**2*EP,csqc1
            U2(I) = U2(I) + CEP2*V1(I) + ANT*V1(J)
     &        + OMGT*(-V1(J+2) + EHT*V1(J+1) - EHT*V1(J-1) + V1(J-2)) 
            V2(I) = V2(I) - CEP2*U1(I) - ANT*U1(J)
     &        - OMGT*(-U1(J+2) + EHT*U1(J+1) - EHT*U1(J-1) + U1(J-2)) 
            U2(J) = U2(J) + CEP2*V1(J) - ANT*V1(I)
     &        - OMGT*(-V1(I+2) + EHT*V1(I+1) - EHT*V1(I-1) + V1(I-2)) 
            V2(J) = V2(J) - CEP2*U1(J) + ANT*U1(I)
     &        + OMGT*(-U1(I+2) + EHT*U1(I+1) - EHT*U1(I-1) + U1(I-2)) 
         ENDDO
c     i=n-1 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = -CSQC1*(XD - QC)/DELQ**2*EP
         U2(N-1) = U2(N-1) + CEP2*V1(N-1) + ANT*V1(J)
     &        + OMGT*(EHT*V1(TWON) - EHT*V1(TWON-2) + V1(TWON-3))
         V2(N-1) = V2(N-1) - CEP2*U1(N-1) - ANT*U1(J)
     &        - OMGT*(EHT*U1(TWON) - EHT*U1(TWON-2) + U1(TWON-3))
         U2(TWON-1) = U2(TWON-1) + CEP2*V1(TWON-1) - ANT*V1(I)
     &        - OMGT*(EHT*V1(N) - EHT*V1(N-2) + V1(N-3))
         V2(TWON-1) = V2(TWON-1) - CEP2*U1(TWON-1) + ANT*U1(I)
     &        + OMGT*(EHT*U1(N) - EHT*U1(N-2) + U1(N-3))
c     i=n 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = -CSQC1*(XD - QC)/DELQ**2*EP
         U2(N) = U2(N) + CEP2*V1(N) + ANT*V1(J)
     &        + OMGT*(-EHT*V1(TWON-1) + V1(TWON-2))
         V2(N) = V2(N) - CEP2*U1(N) - ANT*U1(J)
     &        - OMGT*(-EHT*U1(TWON-1) + U1(TWON-2))
         U2(TWON) = U2(TWON) + CEP2*V1(TWON) - ANT*V1(I)
     &        - OMGT*(-EHT*V1(N-1) + V1(N-2))
         V2(TWON) = V2(TWON) - CEP2*U1(TWON) + ANT*U1(I)
     &        + OMGT*(-EHT*U1(N-1) + U1(N-2))
      ELSEIF(TPCRS(1:7).EQ.'.PDERIV')THEN
c     i=1
         XD = XI(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = -CSQC1*(XD - QC)/DELQ**2*EP
c         ANT = -CSQC1*(EXP(-((XD + SH(1) - QC)/DELQ)**2) 
c     &           - EXP(-((XD - SH(1) - QC)/DELQ)**2))/(2*SH(1)) !(XD - QC)/DELQ*EP
         U2(1) = U2(1) + CEP2*V1(1) 
     &        + OMGT*V1(NP1+1) + ANT*V1(NP1)
         V2(1) = V2(1) - CEP2*U1(1) 
     &        - OMGT*U1(NP1+1) - ANT*U1(NP1)
         U2(NP1) = U2(NP1) + CEP2*V1(NP1)
     &        - OMGT*V1(2) - ANT*V1(1)
         V2(NP1) = V2(NP1) - CEP2*U1(NP1)
     &        + OMGT*U1(2) + ANT*U1(1)
c     i=2:n-1
         DO I=2,N-1,1
            J = I + N
            XD = XD + SH(1)
            EP = EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2*EP*EP
            OMGT = CSQC*EP
c            ANT  = CSQC*(EXP(-((XD + SH(1) - QC)/DELQ)**2) 
c     &           - EXP(-((XD - SH(1) - QC)/DELQ)**2))/TWO 
            ANT = -CSQC1*(XD - QC)/DELQ**2*EP
c            write(*,*)ANT,-CSQC1*(XD - QC)/DELQ**2*EP
            U2(I) = U2(I) + CEP2*V1(I)
     &        + OMGT*(V1(J+1) - V1(J-1)) + ANT*V1(J)
            V2(I) = V2(I) - CEP2*U1(I)
     &        - OMGT*(U1(J+1) - U1(J-1)) - ANT*U1(J)
            U2(J) = U2(J) + CEP2*V1(J)
     &        - OMGT*(V1(I+1) - V1(I-1)) - ANT*V1(I)
            V2(J) = V2(J) - CEP2*U1(J)
     &        + OMGT*(U1(I+1) - U1(I-1)) + ANT*U1(I)
         ENDDO
c     i=n 
         XD = XD + SH(1)
         EP = EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = -CSQC1*(XD - QC)/DELQ**2*EP
         U2(N) = U2(N) + CEP2*V1(N)
     &        - OMGT*V1(TWON-1) + ANT*V1(J)
         V2(N) = V2(N) - CEP2*U1(N)
     &        + OMGT*U1(TWON-1) - ANT*U1(J)
         U2(TWON) = U2(TWON) + CEP2*V1(TWON)
     &        + OMGT*V1(N-1) - ANT*V1(I)
         V2(TWON) = V2(TWON) - CEP2*U1(TWON)
     &        - OMGT*U1(N-1) + ANT*U1(I)
      ELSEIF(TPCRS(1:5).EQ.'.CORD')THEN
         XD = XI(1)
         DO I=1,N,1
            J = N + I
            XD = XD + SH(1)
            OMGT = CSQC*XD
D            OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)
            U2(I) = U2(I) + OMGT*V1(J)
            V2(I) = V2(I) - OMGT*U1(J)
            U2(J) = U2(J) + OMGT*V1(I)
            V2(J) = V2(J) - OMGT*U1(I)
         ENDDO
      ELSE
         XD = XI(1)
         DO I=1,N,1
            J = N + I
            XD = XD + SH(1)
            OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)
            U2(I) = U2(I) + OMGT*V1(J)
            V2(I) = V2(I) - OMGT*U1(J)
            U2(J) = U2(J) + OMGT*V1(I)
            V2(J) = V2(J) - OMGT*U1(I)
         ENDDO
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c     FNIA must be igual to zero
D            FNIA = FCORREL(N, U2(1), V2(1)) 
D     &           - FCORREL(N, V2(1), U2(1))
            FNTA = FNRA
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
c     FNIB must be igual to zero
D            FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &           - FCORREL(N, V2(NP1), U2(NP1))
            FNTB = FNRB
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
            IF(ABSORB)THEN
               CONTINUE
            ELSE
               IF(ABS((FNT-F0)/F0)*1.0D+2.GT.TOL)INFO = ONE
            ENDIF
c
D            write(*,*)FCORREL(N, U2(NP1), U2(NP1)),
D     &           FCORREL(N, V2(NP1), V2(NP1)),
D     &           FCORREL(N, U2(1), U2(1)), 
D     &           FCORREL(N, V2(1), V2(1))
            WRITE(*,1021)T, ET, ET0, CRA, CIA, CTA, CRB, CIB, CTB, 
     &           FNTA, FNTB, FNT  
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
                  CALL PRTPT(NEWNAM4, 9, ND, T, NP, XP, XI, SH, 
     &                 WORK(NP1))
                  CALL PRPT2(NEWNAM6, 8, ND, T, NP, XP, XI, SH, 
     &                 U2(NP1), V2(NP1))
               ENDIF
            ENDIF
            IF(PRTPULS)THEN
               EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
               WRITE(21,*)T,EFAUX         !/SQRT(E0)/2.1132D-9
            ENDIF 
         ENDIF
      ENDDO
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
      ENDIF
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

