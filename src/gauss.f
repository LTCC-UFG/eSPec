      REAL*8 FUNCTION GET_GAUSS(rA, COORD, rnorm)

      IMPLICIT NONE

      !Purpose:
      !creating a Gaussian wave packet

      !Constants
      REAL*8 LN_2, COEF
      PARAMETER (LN_2 = 0.693147D0, COEF = 0.68536D0)

      !Scalar arguments
      REAL*8 rA, COORD, rNORM, TEMP

      !Array arguments
d      REAL*8 rA(*)

      GET_GAUSS = rNORM * DEXP( - 0.5D0*LN_2*(COORD / rA)**2)
c      GET_GAUSS = rNORM * DEXP( - 0.5D0*LN_2*COORD**2)
      
!      TEMP = (1 / rA(2))
!      GET_GAUSS = COEF*DSQRT(TEMP)*
!     &DEXP( - 0.5D0*LN_2*(TEMP * COORD)**2)
      
      END FUNCTION GET_GAUSS
