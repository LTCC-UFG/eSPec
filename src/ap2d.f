      SUBROUTINE AP2D(NP,AP,SHM,VPOT)
      IMPLICIT NONE
c     ** Array arguments
      INTEGER NP(*)
      REAL*8 AP(*),VPOT(*),SHM(*)
c     ** Parameters 
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1)
c     ** Local scalars
      INTEGER N,I,J,L,M
      REAL*8 SHT,APAUX
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*     
*     ..
*     Historic
*     ========
*     (03/02/2003) First version AU_2D written by Freddy 
*     (18/09/2014) small bug fixed by Vinicius and Freddy: The matrix was
*     incorrect near the edges of the grid.
*
*

c The restrictions on L arise from the grid limits, in order not to include
c points in the derivative that are not related to a given point (L,M)
c the restrictions on M appear automatically when we generate only the half
c symmetric matrix

      
      N = NP(1)*NP(2)
      SHT = SHM(1) + SHM(2)
d      SHT = -ONE
d      SHM(1) = -ONE
d      SHM(2) = -ONE
      DO M=1,NP(2),1
         DO L=1,NP(1),1
            J = L + (M-1)*NP(1)
            DO I=1,J,1
c     .. Generate half simetrical 2-dimensional hamiltonian matrix
cccccccccccccccccccccccccccccccccccccccccc
c     dydy 
               IF(J.EQ.(I+2*NP(1)))THEN
                  APAUX = + ONE*SHM(2)
               ELSEIF(J.EQ.I+NP(1))THEN
                  APAUX = - SIXTEEN*SHM(2) 
c     dxdx
               ELSEIF( J.EQ.(I+2) .AND. L.GT.2)THEN
                  APAUX = + ONE*SHM(1)
               ELSEIF(J.EQ.(I+1).AND. L.GT.1)THEN
                  APAUX = - SIXTEEN*SHM(1)
c     dxdx + dydy 
               ELSEIF(I.EQ.J)THEN
                  APAUX = + THIRTY*SHT + VPOT(I)
               ELSE
                  APAUX = ZERO
               ENDIF
               AP(I+(J-1)*J/2) = APAUX
            ENDDO
!            write(2,1111)(int(ap(i+(j-1)*j/2)), i=1,j,1)
         ENDDO
      ENDDO 
 1111 FORMAT(100(I3,1X))
      END
