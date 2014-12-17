      SUBROUTINE OFFDIAG() 
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
      REAL*8        LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
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
      EXTERNAL      PSODAUX, PRTPT
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT 
c     .. Start program      
c     ..
c     .. Starting values
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF(TPCRS(1:6).EQ.'.DERIV')THEN
c     i=1
         XD = XI(1)
         OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)
         ANT = (XD - QC)/DELQ*CSQC1*EXP(-((XD - QC)/DELQ)**2)
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
            ANT = (XD - QC)/DELQ*CSQC1*EXP(-((XD - QC)/DELQ)**2)
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
         ANT = (XD - QC)/DELQ*CSQC1*EXP(-((XD - QC)/DELQ)**2)
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
         ANT = CSQC1*(XD - QC)/DELQ*EP
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
         ANT = CSQC1*(XD - QC)/DELQ*EP
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
            ANT = CSQC1*(XD - QC)/DELQ*EP
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
         ANT = CSQC1*(XD - QC)/DELQ*EP
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
         ANT = CSQC1*(XD - QC)/DELQ*EP
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
         EP = XD!EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = CSQC1!*(XD - QC)/DELQ*EP
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
            EP = XD!EXP(-((XD - QC)/DELQ)**2)
            CEP2 = CSQC2*EP*EP
            OMGT = CSQC*EP
            ANT = CSQC1!*(XD - QC)/DELQ*EP
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
         EP = XD!EXP(-((XD - QC)/DELQ)**2)
         CEP2 = CSQC2*EP*EP
         OMGT = CSQC*EP
         ANT = CSQC1!*(XD - QC)/DELQ*EP
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
D            OMGT = CSQC*EXP(-((XD - QC)/DELQ)**2)*XD
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

