c     .. Read potential from file called INFILE
c     =========================================
      SUBROUTINE RDPTE(INFILE, NT, ND, XI, XF, XMIN, V)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
cdel      LOGICAL
      CHARACTER*(*) INFILE
      INTEGER       INFLENGTH, NT, ND
      REAL*8        XMIN
c     **
c     ** Array arguments
cdel       LOGICAL
cdel       CHARACTER      
cdel       INTEGER  
      REAL*8        XI(*), XF(*), V(*) 
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
*     Historic
*     ========
*     (14/03/2003) First version RDPT written by Freddy. 
*
c     **
c     ** Parameters 
cdel      REAL*8        ZERO, ONE
cdel      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel       LOGICAL
c      CHARACTER*1   
      CHARACTER*72  TITLE
      INTEGER       I, J, IOTEST
      REAL*8        VAUX

c     **
c     ** External functions 
cdel       LOGICAL
cdel       CHARACTER
      INTEGER       ICHLENGTH
c      REAL*8        ECNORM  
c     **
c     ** External subroutines 
cdel       EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel       INTRINSIC   ECNORM    
c     .. Start program
c     ..
      INFLENGTH = ICHLENGTH(INFILE, 0)
c
      OPEN(69, STATUS='OLD', FILE=INFILE(1:INFLENGTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.0)THEN
         WRITE(*,*)'<<<>>> Potetial input file "', INFILE(1:INFLENGTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(69,'(A)',END=999,ERR=999)TITLE
c
      READ(69,*,END=999,ERR=999)(XI(I), I=ND,1,-1), V(1)
      VAUX = V(1)
      DO J=2,NT-1,1
D         write(*,*)j,nt
         READ(69,*,END=999,ERR=999)(XF(I), I=ND,1,-1), V(J) 
         IF(V(J).LT.VAUX)THEN
            VAUX = V(J)
            XMIN = XF(J)
         ENDIF 
      ENDDO 
      READ(69,*,END=999,ERR=999)(XF(I), I=ND,1,-1),V(NT)
      IF(V(NT).LT.VAUX)THEN
         VAUX = V(NT)
         XMIN = XF(NT)
      ENDIF
c
c      ORM = ECNORM(NT, V)
D      DO I=1,NT,1
D         V(I) = V(I) !/ORM
D      ENDDO
c
 992  FORMAT(1X,'End of file, reading finish!')
c      WRITE(*,992)
c     .. 
      RETURN
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
      END


	
	
