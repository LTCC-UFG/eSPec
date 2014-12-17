      PROGRAM Ginp
      IMPLICIT REAL*8(A-h,O-Z)
c
      READ(1,*)N, DW, WI
c
      OMG = (N - 1)*DW + WI
c
cxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xx
c23456789012345678901234567890123456789012345678901234567890123456789012
      WRITE(*,11)'*** eSPec input file ***       '
      WRITE(*,11)'========================'
      WRITE(*,11)'**MAIN'
      WRITE(*,11)'*TITLE'
      WRITE(*,11)'Input to perform IXPP calculations in two coupled fina
     &l states'
      WRITE(*,11)'*DIMENSION'
      WRITE(*,11)'.1D'
      WRITE(*,11)'256/'
      WRITE(*,11)'*POTENTIAL'
      WRITE(*,11)'.FILE'
      WRITE(*,11)'data.inp'
      WRITE(*,11)'*MASS'
      WRITE(*,11)'1.0/'
      WRITE(*,11)'*TPCALC'
      WRITE(*,11)'.PROPAGATION'
      WRITE(*,11)'*INIEIGVC'
      WRITE(*,11)'.CALC'
      WRITE(*,11)'*CHANGE'
      WRITE(*,11)'.YES'
      WRITE(*,11)'*PRTCRL'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTPOT'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTEIGVC'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTVEFF'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTEIGVC2'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTPULSE'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTDIPOLE'
      WRITE(*,11)'.NO'
      WRITE(*,11)''
      WRITE(*,11)'**TI'
      WRITE(*,11)'*TPDIAG'
      WRITE(*,11)'.LANCZS   !.LANCZS; MTRXDIAG '
      WRITE(*,11)'*NIST'
      WRITE(*,11)'1 0/'
      WRITE(*,11)'*ABSTOL'
      WRITE(*,11)'1D-6/'
      WRITE(*,11)''
      WRITE(*,11)'**TD'
      WRITE(*,11)'*PROPAG'
      WRITE(*,11)'.P2PABM'
      WRITE(*,11)'0.0 1000.0 1D-4 10/'
      WRITE(*,11)'.PULSE'
      WRITE(*,11)'100D0 500D0 1.0D0/'
      WRITE(*,11)'.GAMMA_OMEG_TRNSDIPOL'
      WRITE(*,*)0.1D-5, OMG, 1.0D-4
      WRITE(*,11)'*PULSE_AND_DIPOLE'
      WRITE(*,11)'.GGAUS '
c     .. CO IR pulse parameters
c     WRITE(*,11)'5D+12  100.0  0.193  500.0  0.0000  1.0/ '
c     .. 
c     .. N2 IR pulse parameters 
c     WRITE(*,11)'5D+12  100.0  0.236  500.0  0.0000  1.0/ '
c     WRITE(*,11)'5D+13  100.0  0.236  500.0  0.0000  1.0/ '
      WRITE(*,11)'5D+14  100.0  0.236  500.0  0.0000  1.0/ '
c     WRITE(*,11)'5D+15  100.0  0.236  500.0  0.0000  1.0/ '
c     WRITE(*,11)'5D+14  100.0  0.150  500.0  0.0000  1.0/ '
c     WRITE(*,11)'5D+14  100.0  0.202  500.0  0.0000  1.0/ '
c     WRITE(*,11)'5D+14  100.0  0.300  500.0  0.0000  1.0/ '
c     .. 
      WRITE(*,11)'.READ'
c     .. CO dipole monent 
c     WRITE(*,11)'CO-CO+_mdipol.inp'
c     .. N2 dipole monent 
      WRITE(*,11)'dipol_N2.inp'
      WRITE(*,11)'*PRPGSTATE'
      WRITE(*,11)'0/'
      WRITE(*,11)'*TPTRANS'
      WRITE(*,11)'.ONE'
      WRITE(*,11)'*PRPTOL'
      WRITE(*,11)'3D0/'
      WRITE(*,11)'*NPROJECTIONS'
      WRITE(*,11)'6  1'
      WRITE(*,11)''
      WRITE(*,11)'**END'
c     ..
 11   FORMAT(A)
      RETURN
      END
