      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8  VEC(10000), MTR(5,10000)
c
      WRITE(*,*)'DIGITE O NUMERO DE PONTOS E O NUMERO DE DADOS'
      READ(*,*)N, M
      
      
      DO I=1,M,1
         DO J=1,N,1
            READ(1,*)VEC(J), MTR(I,J)
         ENDDO
          READ(1,*)
      ENDDO
c
      DO J=1,N,1
         SUM = 0.0D0
         DO I=1,M,1
            SUM = SUM + MTR(I,J)
         ENDDO
         WRITE(2,*)VEC(J), SUM
      ENDDO
c
      END

      
