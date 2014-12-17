      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8  VEC1(10000), VEC2(10000)
c
      WRITE(*,*)'De o valor de gamma'
      READ(*,*)GAMMA
c
      DO N=1,10000,1
         READ(1,*,END=99)VEC1(N), VEC2(N) 
      ENDDO
 99   CONTINUE
      N = N - 1
      write(*,*)'nf=',n
c     
      DO I=1,N,1
         SUM = 0.0D0
         DO J=1,N,1
            SUM = SUM + VEC2(I)*G(VEC1(I), VEC1(J), GAMMA)
         ENDDO
         WRITE(2,*)VEC1(I), SUM
      ENDDO
c     
      END

      FUNCTION G(WI, WJ, GAMMA)
      IMPLICIT REAL*8(A-H,O-Z)
c
      G = DEXP(-((WJ - WI)/GAMMA)**2)
c
      RETURN
      END
