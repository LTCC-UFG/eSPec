      SUBROUTINE AP2DCT(NP,AP,SHM,VPOT)
      IMPLICIT NONE
c     ** Array arguments
      INTEGER NP(*)
      REAL*8 AP(*),VPOT(*),SHM(*)
c     ** Parameters 
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1)
c     ** Local scalars
      INTEGER I,J,L,M,N
      REAL*8 SHT,APAUX
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*     Vinicius Vaz da Cruz
*     ..
*     Historic
*     ========
*     (03/02/2003) First version AU_2D written by Freddy 
*     (18/09/2014) AU_2DCT written based on AU_2D by Vinicius and Freddy
* 
*     

c This routine generates a Simmetric Hamiltonian matrix for
c H = d2/dx2 + d2/dy2 + d2/dxdy + V
c
      
      N = NP(1)*NP(2)
      SHT = SHM(1) + SHM(2)
      !SHT = -ONE
      !SHM(1) = -ONE
      !SHM(2) = -ONE
      DO M=1,NP(2),1
         DO L=1,NP(1),1
            J = L + (M-1)*NP(1)
            DO I=1,J,1
c     .. Generate half simetrical 2-dimensional hamiltonian matrix
c     .. with a d2/dxdy kinetic energy cross term 
cccccccccccccccccccccccccccccccccccccccccc
cc     dxdy (2,2)
               IF(J.EQ.I+2*NP(1)+2 .AND. L.GT.2 
     &              )THEN
                  APAUX = + ONE*SHM(3)
cc     dxdy (-2,2)  
               ELSEIF(J.EQ.I+2*NP(1)-2 .AND. L.LT.NP(1)-1
     &                 )THEN 
                  APAUX = - ONE*SHM(3)
cc     dxdy (1,1)
               ELSEIF(J.EQ.I+NP(1)+1 .AND. L.GT.1 
     &                 )THEN  
                  APAUX = - SIXTEEN*SHM(3)
cc     dxdy (-1,1)
               ELSEIF(J.EQ.I+NP(1)-1 .AND. L.LT.NP(1) )THEN
                  APAUX = + SIXTEEN*SHM(3)
c-----------
c     dydy 
               ELSEIF(J.EQ.(I+2*NP(1)))THEN
                  APAUX = + ONE*SHM(2)
c      
               ELSEIF(J.EQ.I+NP(1))THEN
                  APAUX = - SIXTEEN*SHM(2) 
c-----------
c     dxdx
               ELSEIF( J.EQ.(I+2) .AND. L.GT.2)THEN
                  APAUX = + ONE*SHM(1)
c     
               ELSEIF(J.EQ.(I+1).AND. L.GT.1)THEN
                  APAUX = - SIXTEEN*SHM(1)
c-----------
c     dxdx + dydy 
               ELSEIF(I.EQ.J)THEN
                  APAUX = + THIRTY*SHT + VPOT(I)
               ELSE
                  APAUX = ZERO
               ENDIF
               AP(I+(J-1)*J/2) = APAUX
            ENDDO
            !write(2,1111)(int(ap(i+(j-1)*j/2)), i=1,j,1)
            !write(2,1112)(ap(i+(j-1)*j/2), i=1,j,1)
         ENDDO
      ENDDO

 1111 FORMAT(100(I3,1X))
 1112 FORMAT(100(ES26.18,1X))
      END
