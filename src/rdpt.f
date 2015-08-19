c     .. Read potential from file called INFILE
c     =========================================
      SUBROUTINE RDPT(SHIFTP, INFILE, NT, ND, XI, XF, XMIN, V)
      IMPLICIT NONE
C 
C vinicius: this one reads potentials and does a shift (V)
C
c     **
c     ** Scalar arguments 
      LOGICAL       SHIFTP
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
cdel       REAL*8
c     **
c     ** External subroutines 
cdel       EXTERNAL      
c     **
c     ** Intrinsic functions 
cdel       INTRINSIC     
c     .. Start program
c     ..
      INFLENGTH = ICHLENGTH(INFILE, 0)
c
      OPEN(69, STATUS='OLD', FILE=INFILE(1:INFLENGTH), IOSTAT=IOTEST)
      IF(IOTEST.NE.0)THEN
         WRITE(*,*)'<<<>>> Input file "', INFILE(1:INFLENGTH), 
     &        '" cannot be opened <<<>>>'
         STOP
      ENDIF
      READ(69,'(A)',END=999,ERR=999)TITLE
c
      READ(69,*,END=999,ERR=999)(XI(I), I=ND,1,-1), V(1)
      IF(SHIFTP)VAUX = V(1)
c*123      DO J=2,NT-1,1
      DO J=2,NT,1
c         write(*,*)'j=',j
         READ(69,*,END=999,ERR=999)(XF(I), I=ND,1,-1), V(J)
c         write(*,*)'j,v(j)',j,v(j)
         IF(SHIFTP .AND. V(J).LT.VAUX)THEN
            VAUX = V(J)
c            write(*,*)'j1,v1(j)',j,v(j)
c array xf cannot have possition bigger than 3
c            XMIN = XF(J)
            xmin = 0
c            write(*,*)'j2,XF(j),xmin',j,xf(1),xmin
         ENDIF 
      ENDDO
c*123      READ(69,*,END=999,ERR=999)(XF(I), I=ND,1,-1),V(NT)
      IF(SHIFTP .AND. V(NT).LT.VAUX)THEN
         VAUX = V(NT)
c array xf cannot have possition bigger than 3
c     XMIN = XF(NT)
         xmin = 0
      ENDIF
c
      IF(SHIFTP)THEN
         DO I=1,NT,1
            V(I) = V(I) - VAUX
         ENDDO
      ENDIF
c     
 992  FORMAT(1X,'End of file, reading finish!')
      WRITE(*,992)
c     .. 
      RETURN
 999  WRITE(*,*)'<<<>>> Reading file error, check and resubmit <<<>>>'
      WRITE(*,*)INFILE
      STOP
      END


	
	
