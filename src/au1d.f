      SUBROUTINE AU_1D(N, SHM, VPOT, US, AU)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      INTEGER       N
cdel      REAL*8       
c     **
c     ** Array arguments
      REAL*8        SHM(*), VPOT(*), US(*), AU(*)    
c     **
*     ..
*     Purpose
*     =======
*
*     ..
*     Arguments
*     =========
*
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (03/02/2003) First version AU_1D written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8          ONE, SIXTEEN, THIRTY
      PARAMETER       (ONE = +1.0D+0, SIXTEEN = +1.6D+1, 
     &     THIRTY = +3.0D+1) 
c     **
c     ** local scalars **
      INTEGER        I
      REAL*8         THIRTYS, SIXTEENS
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*
cdel      INTEGER       
cdel      REAL*8
c     **
c     ** External subroutines 
cdel      EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel      INTRINSIC     
c     .. Start program 
c      write(*,*)'***********',n,us(1)
      THIRTYS  = THIRTY*SHM(1)
      SIXTEENS = SIXTEEN*SHM(1)
c
      AU(1) = (THIRTYS + VPOT(1))*US(1)
     &     - SIXTEENS*US(2) + SHM(1)*US(3)
      AU(2) = - SIXTEENS*US(1)
     &     + (THIRTYS + VPOT(2))*US(2)
     &     - SIXTEENS*US(3) + SHM(1)*US(4)
c paralelizar (estudar numero de processadores para otimizar a paral.)
c     os sinais estao corretos a autoenergia torna-se negativa quando os 
c     sinais sao invertidos.
      DO I=3,N-2,1
          AU(I) = SHM(1)*US(I-2) - SIXTEENS*US(I-1)
     &        + (THIRTYS + VPOT(I))*US(I)
     &        - SIXTEENS*US(I+1) + SHM(1)*US(I+2)
      ENDDO
c
      AU(N-1) = SHM(1)*US(N-3) - SIXTEENS*US(N-2)
     &     + (THIRTYS + VPOT(N-1))*US(N-1)
     &     - SIXTEENS*US(N)
      AU(N) = SHM(1)*US(N-2) - SIXTEENS*US(N-1)
     &     + (THIRTYS + VPOT(N))*US(N)
c     ..
      RETURN
      END

