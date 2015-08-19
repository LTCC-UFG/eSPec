      PROGRAM TESTMTRX
      implicit none
      INTEGER I,J,K,NP(3),N
      REAL*8 APAUX,AP(60000),SHT,SHM(3)
      REAL*8        ZERO, ONE, SIXTEEN, THIRTY
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = 3.0D+1) 

      NP(1)=10
      NP(2)=10
      N=NP(1)*NP(2)
      SHM(1)=-ONE
      SHM(2)=-ONE
      SHM(3)=-ONE

      !SHT = SHM(1) + SHM(2)
      SHT = ONE
      DO J=1,N,1
         DO I=1,J,1
c     .. Generate half simetrical 2-dimensional hamiltonian matrix
c     .. with a d2/dxdy kinetic energy cross term 
d     write(*,*)'This option does not work'
d     stop
cccccccccccccccccccccccccccccccccccccccccc
c dxdy -2,-2
            IF(J.EQ.I+2*NP(1)+2)THEN
               APAUX = + ONE*SHM(3)
c dydy 
            ELSEIF(J.EQ.I+2*NP(1) )THEN
               APAUX = + ONE*SHM(2)
c dxdy 2,-2  
            ELSEIF(J.EQ.I+2*NP(1)-2)THEN 
               APAUX = - ONE*SHM(3)
c dxdy -1,-1
            ELSEIF(J.EQ.I+NP(1)+1 .AND. )THEN  
               APAUX = - SIXTEEN*SHM(3)
c dydy 
            ELSEIF(J.EQ.I+NP(1))THEN
               APAUX = - SIXTEEN*SHM(2) 
c dxdy 1,-1
            ELSEIF(J.EQ.I+NP(1)-1)THEN
               APAUX = + SIXTEEN*SHM(3)
c dxdx
            ELSEIF(J.EQ.I+2)THEN
               APAUX = + ONE*SHM(1)
c dxdx
            ELSEIF(J.EQ.I+1)THEN
               APAUX = - SIXTEEN*SHM(1) 
c dxdx + dydy 
            ELSEIF(I.EQ.J)THEN
               APAUX = + THIRTY*SHT !+ VPOT(I)
            ELSE
               APAUX = ZERO
            ENDIF
            AP(I+(J-1)*J/2) = APAUX 
         ENDDO
         write(1,1111)(int(ap(i+(j-1)*j/2)), i=1,j,1)
      ENDDO

c      DO K = 1,N*N/2
         

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 1111 FORMAT(100(I3,1X))
      STOP
      END
