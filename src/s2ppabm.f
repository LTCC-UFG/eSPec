      SUBROUTINE S2PPABM(ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS, DIM, EFC, 
     &     PRTCRL, TPABSOR, INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, 
     &     NPR, KP, NISPG, KL, KLX, E0, T0, TD, TP, OMG, SNI, DT, 
     &     TI, TF, TOL, GAMMA, OMEGA, DMX, T0X, TPX, NP, ND, SH, SHM, 
     &     U1, U2, V1, V2, XI, XP, XF, DM, VPOT, VABC, WORK, VAR, LNZVC) 
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, PRTVEFF, PRTEIGVC2, PRTPULS 
      CHARACTER*(*) DIM, EFC, PRTCRL, TPABSOR
      INTEGER       INFO, LC, MP, MPI, MPR, MXDCT, N, NSHOT, NPR, KP
      INTEGER       NISPG, KL, KLX
      REAL*8        DT, TI, TF, TOL, E0, T0, TD, TP, OMG, SNI
      REAL*8        GAMMA, OMEGA, DMX, T0X, TPX
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*)
      INTEGER       NP(*), ND(*)
      REAL*8        SH(*), SHM(*), U1(*), U2(*), V1(*), V2(*), XI(*)
      REAL*8        XP(*), DM(*), XF(*), VPOT(*), VABC(*), WORK(*)
      REAL*8        VAR(*)
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
*     (22/03/2005) First version written by Freddy based on ppsod.f
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
      LOGICAL       FLG
      CHARACTER*4   CHNUM   
      CHARACTER*15  NEWNAM2, NEWNAM3, NEWNAM5, NEWNAM6
      CHARACTER*16  NEWNAM1, NEWNAM4
      INTEGER       I, J, MC, MC1, MC2, NP1, KL2, TWON
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
c     ** Local arrays 
cdel      LOGICAL
cdel      CHARACTER*3
cdel      INTEGER       
      REAL*8        U0(2*N), V0(2*N), V(2*N), HU1(2*N), HV1(2*N)
      REAL*8        UU(2*N), VV(2*N)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
      REAL*8        FCORREL, EIGENERG, EF    
c     **
c     ** External subroutines 
      EXTERNAL      PSODAUX, PABMAUX, PRTPT
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
      CSGAMMA2 = CST*GAMMA
      CSGAMMA1 = THREE*CSGAMMA2
      CSOMEGA  = OMEGA
      MC1   = ZERO
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
      DO I=1,N,1
         V(I) = VPOT(I)
         U0(I) = U1(I)
         U2(I) = U0(I) 
         V0(I) = V1(I)
         V2(I) = V0(I) 
      ENDDO
      DO I=NP1,TWON,1
         V(I) = VPOT(I)
         U1(I) = ZERO
         U2(I) = ZERO 
         V1(I) = ZERO
         V2(I) = ZERO 
      ENDDO
c 
      FNRA = FCORREL(N, U2(1), U2(1)) 
     &     + FCORREL(N, V2(1), V2(1))
c     FNI must be igual to zero
D      FNIA = FCORREL(N, U2(1), V2(1)) 
D     &     - FCORREL(N, V2(1), U2(1))
      FNTA = FNRA
      F0A = FNTA
c
      FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &     + FCORREL(N, V2(NP1), V2(NP1))
c     FNI must be igual to zero
D      FNIB = FCORREL(N, U2(NP1), V2(NP1)) 
D     &     - FCORREL(N, V2(NP1), U2(NP1))
      FNTB = FNRB
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
            CALL PRTPT(NEWNAM4, 9, ND, T, NP, XP, XI, SH, WORK(NP1))
            CALL PRPT2(NEWNAM6, 8, ND, T, NP, XP, XI, SH, 
     &           U2(NP1), V2(NP1))
         ENDIF
         IF(PRTPULS)THEN
            OPEN(UNIT=21,STATUS='UNKNOWN',FILE='pulses.dat')
            EFAUX =  EF(EFC, E0, T0, TD, TP, OMG, SNI, KL, T)
            WRITE(21,*)T, EFAUX 
         ENDIF 
      ENDIF
c
c     .. First half step in the propagation (Taylor)
      IF(E0.NE.ZERO)CALL CHPOT(EFC, N, E0, T0, TD, TP, (T+DT)/TWO, OMG, 
     &     SNI, KL, DM(1), V(1), VPOT(1))
      CALL PSODAUX(ABSORB, TPABSOR, DIM, N, CST/FOUR, NP, SHM, 
     &     U1(1), U2(1), V1(1), V2(1), 
     &     VPOT(1), VABC(1), WORK, VAR)
c
c     .. Second half step in the propagation (ABM)
      CALL AU(DIM, N, NP, SHM, VPOT, V1, HU1, VAR)
      CALL AU(DIM, N, NP, SHM, VPOT, U1, HV1, VAR)
      CALL CHPOT(EFC, N, E0, T0, TD, TP, T+DT, OMG, SNI, KL, 
     &     DM(1), V(1), VPOT(1))
      CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST/TWO, THREECST/TWO, 
     &     ZERO, ZERO, NP, SHM, U1(1), U2(1), V1(1), 
     &     V2(1), HU1(1), HV1(1), VPOT(1), VABC(1), WORK, VAR)
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
               CALL PRTPT(NEWNAM4, 9, ND, T, NP, XP, XI, SH, WORK(NP1))
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
      DO I=1,N,1
         U1(I) = U0(I)
         V1(I) = V0(I)
      ENDDO 
      CALL AU(DIM, N, NP, SHM, VPOT, V0, HU1, VAR)
      CALL AU(DIM, N, NP, SHM, VPOT, U0, HV1, VAR)
c
      SUMR = ZERO
      SUMI = ZERO
      FLG = .FALSE.
D      write(*,*)'CSGAMMA',csgamma1,csgamma2,PRTCRL(1:8),lc
D      write(*,*)'OMEGA,DMX',OMEGA,DMX
c      
      DO T=TI+2*DT,TF+DT/2,DT
c
         EFX     = EXP(-((T - T0X)/TPA)**KL2)
         EFXDM   = EFX*DMX 
         CSEFXDM = CST*EFXDM
         OMGT  = CSOMEGA*T
         COMGT = CSEFXDM*DCOS(OMGT)
         SOMGT = CSEFXDM*DSIN(OMGT)
         COMGT1 = CSEFXDM*DCOS(OMGT)/TWO
         SOMGT1 = CSEFXDM*DSIN(OMGT)/TWO
         COMGT3 = THREE/TWO*CSEFXDM*DCOS(OMGT)
         SOMGT3 = THREE/TWO*CSEFXDM*DSIN(OMGT)
c
         DO I=1,N,1
            J = N + I
            UU(I) = COMGT1*V1(J) + SOMGT1*U1(J)
            VV(I) = - COMGT1*U1(J) + SOMGT1*V1(J)
            UU(J) = COMGT1*V1(I) - SOMGT1*U1(I)
            VV(J) = - COMGT1*U1(I) - SOMGT1*V1(I)
         ENDDO     
c
         CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, 
     &        DM(1), V(1), VPOT(1)) 
         CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST, THREECST, ZERO, 
     &        ZERO, NP, SHM, U1(1), U2(1), V1(1), V2(1), 
     &        HU1(1), HV1(1), VPOT(1), VABC(1), WORK, VAR)
c
         CALL CHPOT(EFC, N, E0, T0, TD, TP, T, OMG, SNI, KL, 
     &        DM(NP1), V(NP1), VPOT(NP1))

         CALL PABMAUX(ABSORB, TPABSOR, DIM, N, CST, THREECST, CSGAMMA1, 
     &        CSGAMMA2, NP, SHM, U1(NP1), U2(NP1), V1(NP1), V2(NP1),  
     &        HU1(NP1), HV1(NP1), VPOT(NP1), VABC(NP1), WORK, VAR)
c
D         DO I=1,N,1
D            J = N + I
D            U2(I) = U2(I) - COMGT*V1(J) - SOMGT*U1(J)
D            V2(I) = V2(I) + COMGT*U1(J) - SOMGT*V1(J)
D            U2(J) = U2(J) - COMGT*V1(I) + SOMGT*U1(I)
D            V2(J) = V2(J) + COMGT*U1(I) + SOMGT*V1(I)
D         ENDDO
         DO I=1,N,1
            J = N + I
            U2(I) = U2(I) - COMGT3*V1(J) - SOMGT3*U1(J) + UU(I)
            V2(I) = V2(I) + COMGT3*U1(J) - SOMGT3*V1(J) + VV(I)
            U2(J) = U2(J) - COMGT3*V1(I) + SOMGT3*U1(I) + UU(J)
            V2(J) = V2(J) + COMGT3*U1(I) + SOMGT3*V1(I) + VV(J)
         ENDDO         
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
            CRB = FCORREL(N, U0(1), U2(NP1)) 
     &           + FCORREL(N, V0(1), V2(NP1))
            CIB = FCORREL(N, U0(1), V2(NP1)) 
     &           - FCORREL(N, V0(1), U2(NP1))
            FNRB = FCORREL(N, U2(NP1), U2(NP1)) 
     &           + FCORREL(N, V2(NP1), V2(NP1))
c     FNI must be igual to zero
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
 1021 FORMAT(1X,F14.6,3X,2(E12.6,3X),2(F8.5,3X),7(E10.4,3X))
c
 1031 FORMAT('plot "',A15,'" title `',F9.4,1X,
     &     'fs` with lines lw 2.0, ','"',A14,
     &     '" notitle with lines lw 2.1, ',
     &     '"veffA_0001.dat" notitle with lines lt 3 lw 2.2')

      RETURN
      END
