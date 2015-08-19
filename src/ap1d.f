      SUBROUTINE AP1D(NP,AP,SHM,VPOT)
      IMPLICIT NONE
c     ** Array arguments
      INTEGER NP(*)
      REAL*8 AP(*),VPOT(*),SHM(*)
c     ** Parameters 
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1)
c     ** Local scalars
      INTEGER I,J,N
      REAL*8 SHT,APAUX
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*     
*     ..
*     Historic
*     ========
*     (03/02/2003) First version AU_2D written by Freddy 

      N = NP(1)

c     .. Generate half simetrical 1-dimensional hamiltonian matrix 
      DO J=1,N,1
         DO I=1,J,1
            IF(I.EQ.J)THEN
               APAUX = THIRTY*SHM(1) + VPOT(I)
            ELSEIF(I.EQ.J-1)THEN
               APAUX = -SIXTEEN*SHM(1)
            ELSEIF(I.EQ.J-2)THEN
               APAUX = ONE*SHM(1)
            ELSE
               APAUX = ZERO
            ENDIF
            AP(I+(J-1)*J/2) = APAUX
         ENDDO
      ENDDO

      END
