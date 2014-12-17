      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*15 INPUT, OUTPUT
      REAL*8       VEC(100)
c
      WRITE(*,*)'Input file'
      READ(*,*)INPUT
      WRITE(*,*)'Number of vectors'
      READ(*,*)N
      WRITE(*,*)'Out eigvec'
      READ(*,*)IVC
      WRITE(*,*)'Name output file'
      READ(*,*)OUTPUT
c
      NC = ICHLENGTH(INPUT, 0)
      OPEN(1, STATUS='OLD', FILE=INPUT(1:NC))
      NC = ICHLENGTH(OUTPUT, 0)
      OPEN(2, STATUS='UNKNOWN', FILE=OUTPUT(1:NC))
      L = -1
 5    CONTINUE
      L = L + 1
      IF(L.EQ.0)THEN
         NC = NC + 1
         WRITE(2,*)'|-->  DATA',NC
      ENDIF
      K = 0
 10   CONTINUE
      READ(1,*,ERR=5,END=20)(VEC(I), I=1,N,1)
      K = K + 1
      L = -1
      IF(K.EQ.1)VECAUX = VEC(1)
      IF(VEC(1).NE.VECAUX)THEN
         WRITE(2,*)''
         VECAUX = VEC(1)
         K = 0
      ENDIF
      WRITE(2,*)VEC(1),VEC(2),VEC(IVC)
      GOTO 10
 20   CONTINUE
c
      STOP
 9999 WRITE(*,*)'<<<>>> Input file error <<<>>>' 
      END
