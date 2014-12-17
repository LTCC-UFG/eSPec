c     .. 1 dimensional potentials
c     ===========================
      FUNCTION POT1D(POTCH, CST, RX)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
cdel      CHARACTER*1
      CHARACTER*(*) POTCH         
cdel      INTEGER  
      REAL*8        RX
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*) 
cdel      INTEGER       
      REAL*8        CST(*)      
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
*     (03/02/2003) First version POT1D written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO !, A0A, A0CGS, CCGS, TAU, TWOPI
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, TWO = +2.0D+0) 
!     &     A0A = +5.291772083D-1, 
!     &     A0CGS = +5.291772083D-9, CCGS = +2.99792458D+10, 
!     &     TAU = 2.4188843243D-17, TWOPI = +6.283185307D0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
      REAL*8        POT1D
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
      INTRINSIC     SQRT
c     .. Start program      
      IF(POTCH(1:4).EQ.'.OHS')THEN
c     .. One-dimensional harmonic oscillator
         POT1D = 5.0D-1*CST(2)*(RX - CST(1))**2
      ELSEIF(POTCH(1:5).EQ.'.DOHS')THEN
c     .. One-dimensional double minimum harmonic oscillator
         POT1D = 5.0D-1*CST(2)*(RX - CST(1))**2 
     &        + CST(5)*DEXP((RX - CST(3))**2/CST(4)**2)
      ELSEIF(POTCH(1:7).EQ.'.MORSEM')THEN
c     .. Morse potential
c     G.E. Kellerhals, N. Sathyamurthy, and L. M. Raff, J. Chem. Phys, 64 (1976) 818 
c     pot1d=D*(1-exp(-B*(rij-r01)))**2
D         write(*,*)potch
D         write(*,*)cst(3),cst(2),cst(1)
D         read(*,*)
         POT1D = CST(3)*(
     &          DEXP(-TWO*CST(2)*(RX - CST(1))) 
     &        - TWO*DEXP(-CST(2)*(RX - CST(1)))
     &        )
d         write(*,*)'.MORSEM',RX,POT1D
D         write(*,*)'(RX - CST(1))',(RX - CST(1))
d         write(*,*)'DEXP',DEXP(-TWO*CST(2)*(RX - CST(1))), 
d     &        DEXP(-CST(2)*(RX - CST(1)))
      ELSEIF(POTCH(1:12).EQ.'.ANTI-MORSEM')THEN
c     .. Anti-Morse potential
c     G.E. Kellerhals, N. Sathyamurthy, and L. M. Raff, J. Chem. Phys, 64 (1976) 818 
c     pot1d=D*(1-exp(-B*(rij-r01)))**2
D         write(*,*)potch
D         write(*,*)cst(3),cst(2),cst(1)
D         read(*,*)
         POT1D = CST(3)/TWO*(
     &          DEXP(-TWO*CST(2)*(RX - CST(1))) 
     &        + TWO*DEXP(-CST(2)*(RX - CST(1)))
     &        )
d         write(*,*)'.ANTI-MORSEM',RX,POT1D
D         write(*,*)'(RX - CST(1))',(RX - CST(1))
d         write(*,*)'DEXP2',-TWO*CST(2)*(RX - CST(1)),
d     &        -CST(2)*(RX - CST(1)) 
      ELSEIF(POTCH(1:10).EQ.'.MOD-MORSE')THEN
c     .. Modified Morse potential
         POT1D = CST(3)*(CST(4) - DEXP(-CST(2)*(RX - CST(1))))**2 
     &        + CST(5)
      ELSEIF(POTCH(1:6).EQ.'.MORSE')THEN
c     .. Morse potential
c     P. M. Morse, Phys. Rev., 34 (1929) 57  
c     pot1d=D*(1-exp(-B*(rij-r01)))**2
         POT1D = CST(3)*(ONE - DEXP(-CST(2)*(RX - CST(1))))**2
      ELSEIF(POTCH(1:11).EQ.'.ANTI-MORSE')THEN
c     .. Anti-Morse potential
c     pot1d=D*(1-exp(-B*(rij-r01)))**2
         POT1D = CST(3)/TWO*(ONE + DEXP(-CST(2)*(RX - CST(1))))**2 
      ELSEIF(POTCH(1:7).EQ.'.DMORSE')THEN
c     .. Double minimum Morse potential
c     pot1d=D*(1-exp(-B*(rij-r01)))**2
c     &     +A*exp(-C*(rij-r02)**2)
         POT1D = CST(3)*(ONE - DEXP(-CST(4)*(RX - CST(1))))**2
     &        + CST(5)*DEXP(-CST(6)*(RX - CST(2))**2)
      ELSEIF(POTCH(1:5).EQ.'.WELL')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2))THEN
            POT1D = ZERO
         ELSE
            POT1D = CST(3)
         ENDIF
      ELSEIF(POTCH(1:8).EQ.'.BARRIER')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2))THEN
            POT1D = CST(3)
         ELSE
            POT1D = ZERO
         ENDIF
      ELSEIF(POTCH(1:6).EQ.'.CONST' .OR. POTCH(1:5).EQ.'.FREE')THEN
         POT1D = CST(1)
      ELSE
         WRITE(*,1011)
         WRITE(*,*)POTCH
         STOP
      ENDIF
 1011 FORMAT('<<<>>> Potential function error <<<>>>')
c
      RETURN
      END

c     .. 2 dimensional potentials
c     ===========================
      FUNCTION POT2D(POTCH, SM, CST, RX, RY)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) POTCH      
cdel      INTEGER  
      REAL*8        RX, RY
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*) 
cdel      INTEGER 
      REAL*8        SM(*), CST(*)
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
*     (03/02/2003) First version POT2D written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE, TWO
      PARAMETER     (ZERO = +0.0D+0, 
     &     ONE = +1.0D+0, 
     &     TWO = +2.0D+0
     &     )
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
      REAL*8        POT2D, V1A, V1B, V2, RXA, RYA, RZA, VAB, VBC, RZ
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
      REAL*8        POT3D
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
c      write(*,*)potch
      IF(POTCH(1:4).EQ.'.OHS')THEN
c     .. Two-dimensional harmonic oscillator
         POT2D = 5.0D-1*CST(2)*(RX - CST(1))**2 
     &        + 5.0D-1*CST(4)*(RY - CST(3))**2
      ELSEIF(POTCH(1:6).EQ.'.MORSE')THEN
c     .. Two-dimensional Morse potential 
c     POT1D = D*(1.0d+0-exp(-B*(rij-r01)))**2 + D1*(1.0d+0-exp(-B1*(rij-r01)))**2   
         POT2D = CST(3)*(ONE - DEXP(-CST(2)*(RX - CST(1))))**2
     &        + CST(6)*(ONE - DEXP(-CST(5)*(RY - CST(4))))**2
      ELSEIF(POTCH(1:4).EQ.'.LSM')THEN
c     .. Modified Lippincott-Schroeder potential 
c     [3] H. J. Bakker and H. -K Nienhuys, Science, 297 (2002) 588
c     POT2D = D1A*(ONE - DEXP(-AL1A*(R1 - R1E)**2/(TWO*R1))) + D1B*(ONE - DEXP(-AL1B*(R2 - R1 - R1E)**2/(TWO*(R2 - R1E)))) + D2*(ONE - DEXP(-AL2*(R2 - R2E)))**2
         V1A = CST(3)*(ONE - DEXP(-CST(2)*(RX - CST(1))**2/(TWO*RX 
     &        + 1.0D-99)))
         V1B = CST(5)*(ONE - DEXP(-CST(4)*(RY - RX - CST(1))**2
     &        /(TWO*(RY - CST(1)) + 1.0D-99)))
         V2 = CST(2)*(ONE - DEXP(-CST(7)*(RY - CST(6))))**2
         POT2D = V1A + V1B + V2 
      ELSEIF(POTCH(1:3).EQ.'.LS')THEN
c     .. Lippincott-Schroeder potential
c     [2] H. J. Bakker and H. -K Nienhuys, Science, 297 (2002) 588
c     POT2D = D1A*(ONE - DEXP(-AL1A*(R1 - R1E)**2/(TWO*R1))) 
c     + D1B*(ONE - DEXP(-AL1B*(R2 - R1 - R1E)**2/(TWO*(R2 - R1E)))) 
         V1A = CST(3)*(ONE - DEXP(-CST(2)*(RX - CST(1))**2/(TWO*RX 
     &        + 1.0D-99)))
         V1B = CST(5)*(ONE - DEXP(-CST(4)*(RY - RX - CST(1))**2
     &        /(TWO*(RY - CST(1)) + 1.0D-99)))
         POT2D = V1A + V1B 
      ELSEIF(POTCH(1:7).EQ.'.LBH_CC')THEN
c     .. Leforestier-Bergeron-Hiberty potential 
c     C. Leforestier, G. Bergeron, P.C. Hiberty, Chem. Phys. Lett., 84 (1981) 385
         RYA = RY/SQRT(SM(2)*SM(3)*(SM(1) + SM(2) + SM(3))/
     &        ((SM(2) + SM(3))**2*SM(1)))
         RXA = RX + SM(2)*RYA/(SM(2) + SM(3)) - RYA 
c        
         VAB = CST(1)*DEXP(-ONE/CST(2)*RXA)
         VBC = CST(3)*(DEXP(-TWO*CST(4)*(RYA - CST(5))) 
     &        - TWO*DEXP(-CST(4)*(RYA - CST(5))))
         POT2D = VAB + VBC 
      ELSEIF(POTCH(1:4).EQ.'.LBH')THEN
c     .. Leforestier-Bergeron-Hiberty potential 
c     C. Leforestier, G. Bergeron, P.C. Hiberty, Chem. Phys. Lett., 84 (1981) 385
         VAB = CST(1)*DEXP(-ONE/CST(2)*RX)
         VBC = CST(3)*(DEXP(-TWO*CST(4)*(RY - CST(5))) 
     &        - TWO*DEXP(-CST(4)*(RY - CST(5))))
         POT2D = VAB + VBC 


      ELSEIF(POTCH(1:8).EQ.'.LEPS_CC')THEN
c     .. LEPS potential
c     G.E. Kellerhals, N. Sathyamurthy, and L. M. Raff, J. Chem. Phys, 64 (1976) 818 
c         write(*,*)sm(1),sm(2),sm(3)
         RYA = RY/SQRT(SM(2)*SM(3)*(SM(1) + SM(2) + SM(3))/
     &        ((SM(2) + SM(3))**2*SM(1)))
         RXA = RX - SM(2)*RYA/(SM(2) + SM(3))
         RZA = SQRT(RXA**2 + RYA**2 - TWO*RXA*RYA*COS(CST(13)))
c
D         RYA = RY/SQRT(SM(2)*SM(3)*(SM(1) + SM(2) + SM(3))/
D     &        ((SM(2) + SM(3))**2*SM(1)))
D         RXA = RX + SM(2)*RYA/(SM(2) + SM(3)) - RYA
D         RZA = RXA + RYA 
D         RXA = RY/SQRT(SM(2)*SM(3)*(SM(1) + SM(2) + SM(3))/
D     &        ((SM(2) + SM(3))**2*SM(1)))
D         RYA = RX + SM(2)*RXA/(SM(2) + SM(3)) - RXA
c
c
         POT2D = POT3D(POTCH, CST, RXA, RYA, RZA)
c
      ELSEIF(POTCH(1:5).EQ.'.LEPS')THEN
c     .. LEPS potential
c     G.E. Kellerhals, N. Sathyamurthy, and L. M. Raff, J. Chem. Phys, 64 (1976) 818 
         RZ = SQRT(RX**2 + RY**2 - TWO*RX*RY*COS(CST(13)))
c
         POT2D = POT3D(POTCH, CST, RX, RY, RZ)
c
      ELSEIF(POTCH(1:7).EQ.'.WELLXY')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2))THEN
            POT2D = ZERO
         ELSE
            POT2D = CST(3)
         ENDIF
         IF(RY.GE.CST(4) .AND. RY.LE.CST(5))THEN
            POT2D = ZERO + POT2D
         ELSE
            POT2D = CST(6)+ POT2D
         ENDIF 
      ELSEIF(POTCH(1:6).EQ.'.WELLX')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2))THEN
            POT2D = ZERO
         ELSE
            POT2D = CST(3)
         ENDIF
      ELSEIF(POTCH(1:6).EQ.'.WELLY')THEN
         IF(RY.GE.CST(1) .AND. RY.LE.CST(2))THEN
            POT2D = ZERO
         ELSE
            POT2D = CST(3)
         ENDIF 
      ELSEIF(POTCH(1:5).EQ.'.WELL')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2) .AND.
     &        RY.GE.CST(3) .AND. RY.LE.CST(4))THEN
            POT2D = ZERO
         ELSE
            POT2D = CST(5)
         ENDIF
      ELSEIF(POTCH(1:10).EQ.'.BARRIERXY')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2))THEN
            POT2D = CST(3)
         ELSE
            POT2D = ZERO
         ENDIF
         IF(RY.GE.CST(4) .AND. RY.LE.CST(5))THEN
            POT2D = CST(6) + POT2D
         ELSE
            POT2D = ZERO + POT2D
         ENDIF 
      ELSEIF(POTCH(1:9).EQ.'.BARRIERX')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2))THEN
            POT2D = CST(3)
         ELSE
            POT2D = ZERO
         ENDIF
      ELSEIF(POTCH(1:9).EQ.'.BARRIERY')THEN  
         IF(RY.GE.CST(1) .AND. RY.LE.CST(2))THEN
            POT2D = CST(3)
         ELSE
            POT2D = ZERO
         ENDIF      
      ELSEIF(POTCH(1:8).EQ.'.BARRIER')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2) .AND.
     &        RY.GE.CST(3) .AND. RY.LE.CST(4))THEN
            POT2D = CST(5)
         ELSE
            POT2D = ZERO
         ENDIF
      ELSEIF(POTCH(1:6).EQ.'.CONST' .OR. POTCH(1:5).EQ.'.FREE')THEN
         POT2D = CST(1)
      ELSE
         WRITE(*,1011)
         WRITE(*,*)POTCH
         STOP
      ENDIF
 1011 FORMAT('<<<>>> Potential function error <<<>>>')
c     
      RETURN
      END

c     .. 3 dimensional potentials
c     =========================
      FUNCTION POT3D(POTCH, CST, RX, RY, RZ)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) POTCH       
cdel      INTEGER  
      REAL*8        RX, RY, RZ
c     **
c     ** Array arguments
cdel      LOGICAL
cdel      CHARACTER*(*) POTCH      
cdel      INTEGER       CST(*)
      REAL*8        CST(*)      
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
*     (03/02/2003) First version POT3D written by Freddy. 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
      REAL*8        POT3D, A, B, C, Q1, Q2, Q3, J1, J2, J3

c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
      REAL*8        POT1D
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program      
      IF(POTCH(1:4).EQ.'.OHS')THEN
c     .. three-dimensional harmonic oscillator
         POT3D = 5.0D-1*CST(2)*(RX - CST(1))**2 
     &        + 5.0D-1*CST(4)*(RY - CST(3))**2 
     &        + 5.0D-1*CST(6)*(RZ - CST(5))**2
      ELSEIF(POTCH(1:5).EQ.'.WELL')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2) .AND.
     &        RY.GE.CST(3) .AND. RY.LE.CST(4) .AND.
     &        RZ.GE.CST(5) .AND. RZ.LE.CST(6))THEN
            POT3D = ZERO
         ELSE
            POT3D = CST(7)
         ENDIF
      ELSEIF(POTCH(1:8).EQ.'.BARRIER')THEN
         IF(RX.GE.CST(1) .AND. RX.LE.CST(2) .AND.
     &        RY.GE.CST(3) .AND. RY.LE.CST(4) .AND.
     &        RY.GE.CST(5) .AND. RY.LE.CST(6))THEN
            POT3D = CST(7)
         ELSE
            POT3D = ZERO
         ENDIF
      ELSEIF(POTCH(1:5).EQ.'.LEPS')THEN
c     .. LEPS potential
c     G.E. Kellerhals, N. Sathyamurthy, and L. M. Raff, J. Chem. Phys, 64 (1976) 818 
c
D         write(*,*)cst(1),cst(2),cst(3)
D         write(*,*)cst(4),cst(5),cst(6)
D         write(*,*)cst(7),cst(8),cst(9)
D         write(*,*)cst(10),cst(11),cst(12)
D         read(*,*)
d         write(*,*)'rx,ry,rz',rx,ry,rz
         A = CST(10)
         B = CST(11)
         C = CST(12)
c
         Q1 = 5.0D-1*(
     &        POT1D('.MORSEM', CST(1), RX)
     &        + (1 - A)/(1 + A)*POT1D('.ANTI-MORSEM', CST(1), RX)
     &        )

         Q2 = 5.0D-1*(
     &        POT1D('.MORSEM', CST(4), RY) 
     &        + (1 - B)/(1 + B)*POT1D('.ANTI-MORSEM', CST(4), RY)
     &        )

         Q3 = 5.0D-1*(
     &        POT1D('.MORSEM', CST(7), RZ) 
     &        + (1 - C)/(1 + C)*POT1D('.ANTI-MORSEM', CST(7), RZ)
     &        )
c
         J1 = 5.0D-1*(
     &        POT1D('.MORSEM', CST(1), RX) 
     &        - (1 - A)/(1 + A)*POT1D('.ANTI-MORSEM', CST(1), RX)
     &        )

         J2 = 5.0D-1*(
     &        POT1D('.MORSEM', CST(4), RY) 
     &        - (1 - B)/(1 + B)*POT1D('.ANTI-MORSEM', CST(4), RY)
     &        )

         J3 = 5.0D-1*(
     &        POT1D('.MORSEM', CST(7), RZ) 
     &        - (1 - C)/(1 + C)*POT1D('.ANTI-MORSEM', CST(7), RZ)
     &        )
c
d         write(*,*)'j1,j2,j3',j1,j2,j3
d         write(*,*)'q1,q2,q3',q1,q2,q3
         POT3D = Q1 + Q2 + Q3 - SQRT(
     &        J1**2 + J2**2 + J3**2 
     &        - J1*J2 - J2*J3 - J1*J3
     &        )
d         write(*,*)'pot3d',pot3d,Q1 + Q2 + Q3,J1**2 + J2**2 + J3**2
d     &        - J1*J2 - J2*J3 - J1*J3
d         read(*,*)
c
      ELSEIF(POTCH(1:6).EQ.'.CONST' .OR. POTCH(1:5).EQ.'.FREE')THEN
         POT3D = CST(1)
      ELSE
         WRITE(*,1011)
         WRITE(*,*)POTCH
         STOP
      ENDIF
 1011 FORMAT('<<<>>> Potential function error <<<>>>')
c
      RETURN
      END

