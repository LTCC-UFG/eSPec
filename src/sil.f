c     subroutine PLNZAUX 
c     ==================
      SUBROUTINE PLNZAUX(ABSORB, CHANGE, DIM, TPABSOR, INFO, LMTREORT, 
     &     MP, MXDCT, N, ABSTOL, CST, NP, IWORK, SHM, VPOT, VABC, U, V, 
     &     WU, WV, WP, WORK, VAR, LNZVC)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      LOGICAL       ABSORB, CHANGE
      CHARACTER*(*) DIM, TPABSOR           
      INTEGER       INFO, LMTREORT, MP, MXDCT, N
      REAL*8        ABSTOL, CST
c     **
c     ** Array arguments
cdel      LOGICAL
      INTEGER       NP(*), IWORK(*)
      REAL*8        SHM(*), VPOT(*), VABC(*), U(*), V(*), WU(*) 
      REAL*8        WV(*), WP(*), WORK(*), VAR(*), LNZVC(MXDCT,*)
c     **
*     ..
*     Purpose
*     =======
*     Solve:
*     e^(2*i*|H*dt/hbar)*|p =  |L'*|Q'*e^(-i*|D*dt/hbar)*|Q*|L*|p  
*     ..
*     Arguments
*     =========
*     MP number of lanczos vector to be used.
*     ..
*     Authors
*     =======
*     Freddy Fernandes Guimaraes
*
*     ..
*     Historic
*     ========
*     (19/03/2003) First version PLNZAUX written by Freddy 
*
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
cdel      LOGICAL
cdel      CHARACTER*1
      INTEGER        I, M, IU, IL
      REAL*8         VU, VL, EIAUX
c     **
c     ** Local arrays
cdel      LOGICAL
cdel      CHARACTER*1
cdel      INTEGER  
      REAL*8         AL(MP), BT(MP), EIGVL(MP)
      REAL*8         EIGVC(MP,MP)
c     **
c     ** External functions 
cdel      LOGICAL
cdel      CHARACTER*(*)
cdel      INTEGER      
cdel      REAL*8
c     **
c     ** External subroutines 
      EXTERNAL       LNZ, DSTEVX, DGEMV       
c     **
c     ** Intrinsic functions 
      INTRINSIC      DCOS, DSIN
c     .. Start program
c     .. initiate values
      DO I=1,N,1
         WU(I) = U(I)
         WP(I) = V(I)
      ENDDO
c     ..
c     .. real
c     .. calcula |L
c      call lnz(dim, ndim, n, MP, MP, AC, ni, SH, vpot, wv, wu,
c     &     work, AL, BT, lnzvc)
      CALL LNZ(CHANGE, DIM, N, MP, MP, LMTREORT, MXDCT, NP, SHM,
     &      VPOT, WV, WU, WORK, AL, BT, VAR, LNZVC)
c     .. calcula |Q
      CALL DSTEVX('V', 'A', MP, AL, BT, VL, VU, IL, IU, ABSTOL,
     &     M, EIGVL, EIGVC, MP, WORK, IWORK, INFO)
c     .. testa matriz |Q
      IF(INFO.EQ.ZERO)THEN
c     .. |u1 = |L*|p
         CALL DGEMV( 'T', N, MP, ONE, LNZVC, MXDCT, U, 1,
     &        ZERO, WORK, 1 )
c     .. |u2 = |q*|u1
         CALL DGEMV( 'T', MP, MP, ONE, EIGVC, MP, WORK, 1,
     &        ZERO, U, 1 )
c     ..
c     .. real and complex
c     .. |u3 = e^(-i*|D*dt/hbar)*|u2
         DO I=1,MP,1
            EIAUX = EIGVL(I)*CST
d            write(15,*)i,EIGVL(I)
            U(I) = DCOS(EIAUX)*U(I)
            WV(I) = -DSIN(EIAUX)*U(I)
         ENDDO
d         write(15,*)
c     ..
c     .. real 
c     .. |u4 = |Q'*|u3
         CALL DGEMV( 'N', MP, MP, ONE, EIGVC, MP, U, 1,
     &        ZERO, WORK, 1 )
c     .. |u = |Q'*|u4
         CALL DGEMV( 'N', N, MP, ONE, LNZVC, MXDCT, WORK, 1,
     &        ZERO, U, 1 )  
c     ..
c     .. complex 
c     .. |u4 = |Q'*|u3
         CALL DGEMV( 'N', MP, MP, ONE, EIGVC, MP, WV, 1,
     &        ZERO, WORK, 1 )
c     .. |u = |Q'*|u4
         CALL DGEMV( 'N', N, MP, ONE, LNZVC, MXDCT, WORK, 1,
     &        ZERO, WV, 1 )   
      ELSE
         WRITE(*,1011)
         STOP
      ENDIF
c     ..
c     .. complex 
c     .. calcula |L
c      call lnz(dim, ndim, n, MP, MP, AC, ni, SH, vpot, wu, wp,
c     &     work, AL, BT, lnzvc)
      CALL LNZ(CHANGE, DIM, N, MP, MP, LMTREORT, MXDCT, NP, SHM,
     &      VPOT, WU, WP, WORK, AL, BT, VAR, LNZVC)
c     .. calcula |Q
      CALL DSTEVX('V', 'A', MP, AL, BT, VL, VU, IL, IU, ABSTOL,
     &     M, EIGVL, EIGVC, MP, WORK, IWORK, INFO)
c     ..testa matriz |Q
      IF(INFO.EQ.ZERO)THEN
c     .. |u1 = |L*|p
         CALL DGEMV( 'T', N, MP, ONE, LNZVC, MXDCT, V, 1,
     &        ZERO, WORK, 1 )
c     .. |u2 = |Q*|u1
         CALL DGEMV( 'T', MP, MP, ONE, EIGVC, MP, WORK, 1,
     &        ZERO, V, 1 )
c     ..
c     .. real and complex
c     .. |u3 = e^(-i*|D*dt/hbar)*|u2
         DO I=1,MP,1
            EIAUX = EIGVL(I)*CST
d            write(16,*)i,EIGVL(I)
            WU(I) = DSIN(EIAUX)*V(I)
            V(I) = DCOS(EIAUX)*V(I)
         ENDDO
d         write(16,*)
c     ..
c     .. real 
c     .. |u4 = |Q'*|u3
         CALL DGEMV( 'N', MP, MP, ONE, EIGVC, MP, V, 1,
     &        ZERO, WORK, 1 )
c     .. |u = |Q'*|u4
         CALL DGEMV( 'N', N, MP, ONE, LNZVC, MXDCT, WORK, 1,
     &        ZERO, V, 1 )
c     ..
c     .. complex 
c     .. |u4 = |Q'*|u3
         CALL DGEMV( 'N', MP, MP, ONE, EIGVC, MP, WU, 1,
     &        ZERO, WORK, 1 )
c     .. |u = |Q'*|u4
         CALL DGEMV( 'N', N, MP, ONE, LNZVC, MXDCT, WORK, 1,
     &        ZERO, WU, 1 )  
      ELSE
         WRITE(*,1011)
         STOP
      ENDIF
c
      IF(ABSORB .AND. TPABSOR(1:8).EQ.'.SMOOTHW')THEN
         DO I=1,N,1
            U(I) = VABC(I)*(U(I) + WU(I))
            V(I) = VABC(I)*(WV(I) + V(I))
         ENDDO
      ELSE
         DO I=1,N,1
            U(I) = U(I) + WU(I) 
            V(I) = WV(I) + V(I)
         ENDDO
      ENDIF
c     ..
      RETURN
 1011 FORMAT('<<<>>> Diagonalization error <<<>>>')
      END
