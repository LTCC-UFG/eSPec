      PROGRAM Ginp
      IMPLICIT REAL*8(A-h,O-Z)
c
      PST = 5D-5
      READ(1,*)N, DT
c
      TF = (N - 1)*DT 
c      
      WRITE(*,11)'*** eSPec input file ***       '
      WRITE(*,11)'========================'
      WRITE(*,11)'**MAIN'
      WRITE(*,11)'*TITLE'
      WRITE(*,11)'Input to perform IXPP calculations in two coupled fina
     &l states'
      WRITE(*,11)'*DIMENSION'
      WRITE(*,11)'.1D'
      WRITE(*,11)'310/'
      WRITE(*,11)'*POTENTIAL'
      WRITE(*,11)'.FILE'
      WRITE(*,11)'data.inp'
      WRITE(*,11)'*MASS'
      WRITE(*,11)'1.0D0/'
      WRITE(*,11)'*TPCALC'
      WRITE(*,11)'.PROPAGATION'
      WRITE(*,11)'*INIEIGVC'
      WRITE(*,11)'.GETC'
      WRITE(*,11)'*CHANGE'
      WRITE(*,11)'.NO'
      WRITE(*,11)'*PRTCRL'
      WRITE(*,11)'.LASTWP'
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
      WRITE(*,11)'**TD'
      WRITE(*,11)'*PROPAG'
      WRITE(*,11)'.P2PSOD' 
      WRITE(*,*)'0.0 ', TF, -PST,'/'
      WRITE(*,11)'.AA_QC_DELQ(au)'
      WRITE(*,11)'20.0  2.18499584  0.24/' 
      WRITE(*,11)'*TPTRANS'
      WRITE(*,11)'.ONE'
      WRITE(*,11)'*PRPTOL'
      WRITE(*,11)'3D0/'
      WRITE(*,11)'*NSHOT'
      WRITE(*,11)'2/'
      WRITE(*,11)'**END'
      WRITE(*,11)''
c     ..
 11   FORMAT(A)
      RETURN
      END
cxxxxxxxx1xxxxxxxxx2xxxxxxxxx3xxxxxxxxx4xxxxxxxxx5xxxxxxxxx6xxxxxxxxx7xx
c23456789012345678901234567890123456789012345678901234567890123456789012
