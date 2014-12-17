      SUBROUTINE NWLET(A, NN, NDIM, WKSP1, WKSP2, ISIGN)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      INTEGER       ISIGN, NDIM
c     del      REAL*8        wtstep
c     **
c     ** Array arguments
      INTEGER       NN(*)
      REAL*8        A(*), WKSP1(*), WKSP2(*)
c     
c     **
*     ..
*     Purpose
*     =======
*     Perform the discrete wavelet transform, 
*     
*     ..
*     Arguments
*     =========
*     A is a real array of length equal to the product of these lengths, in 
*     which the data are stored as in a multidimensional real FORTRAN array. 
*  
*     NN is an integer array of length ndim, containing the lengths of each 
*     dimension (number of real values), which MUST all be powers of 2
*
*     NDIM is the number of dimensions.
*
*     ISIGN if isign is input as 1 A is replaced by its wavelet transform.
*     If isign is input as -1 A is replaced by its inverse wavelet transform
*
*     The subroutine wtstep, whose actual name must be supplied in calling 
*     this routine, is the underlying wavelet filter. Examples of wtstep are 
*     daub4 and (preceded by pwtset) pwt
*
*     ..
*     Authors
*     =======
*     
*     ..
*     Historic
*     ========
*     
c     **
c     ** Parameters 
      REAL*8        ZERO, ONE
      PARAMETER     (ZERO = +0.0D+0, ONE = +1.0D+0)
c     **
c     ** Local scalars 
      INTEGER       I1, I2, I3, IDIM, K, N, NNEW, NPREV, NT, NTOT
c     **
c     ** Local arrays
c     **
c     ** External functions 
c     **
c     ** External subroutines 
      EXTERNAL      WTSTEP
c     **
c     ** Intrinsic functions 
c     .. Start program      
      NTOT=1
      DO IDIM=1,NDIM
         NTOT=NTOT*NN(IDIM)
      ENDDO
      NPREV=1
      DO IDIM=1,NDIM            
         N=NN(IDIM)
         NNEW=N*NPREV
         IF (N.GT.4) THEN
            DO I2=0,NTOT-1,NNEW
               DO I1=1,NPREV
                  I3=I1+I2
                  DO K=1,N   
                     WKSP1(K)=A(I3) 
                     I3=I3+NPREV
                  ENDDO  
                  IF (ISIGN.GE.0) THEN            
                     NT=N
 1                   IF (NT.GE.4) THEN
                        CALL WTSTEP(WKSP1, NT, WKSP2, ISIGN)
                        NT=NT/2
                        GOTO 1
                     ENDIF
                     
                  ELSE                           
                     NT=4
 2                   IF (NT.LE.N) THEN
                        CALL WTSTEP(WKSP1, NT, WKSP2, ISIGN)
                        NT=NT*2
                        GOTO 2
                     ENDIF
                  ENDIF
                  I3=I1+I2
                  DO K=1,N                    
                     A(I3)=WKSP1(K)
                     I3=I3+NPREV
                  ENDDO
               ENDDO 
            ENDDO 
         ENDIF
         NPREV=NNEW
      ENDDO
C     ..   
      RETURN
      END


      SUBROUTINE WTSTEP(A, N, WKSP, ISIGN)
      IMPLICIT NONE
c     **
c     ** Scalar arguments 
      INTEGER       N, ISIGN
c     **
c     ** Array arguments
      REAL*8        A(*), WKSP(*)     
c     **
c
*     ..
*     Purpose
*     =======
*     Applies the Daubechies 4-coefficient wavelet filter to data vector 
*     a(1:n) (for isign=1) or applies its transpose (for isign=-1). Used 
*     hierarchically by routines wt1 and wtn.
*     
*     ..
*     Arguments
*     =========
*     
*     ..
*     Authors
*     =======
*     
*     ..
*     Historic
*     ========
*     
c     **
c     ** Parameters 
      REAL*8        C0, C1, C2, C3
      PARAMETER     (C0=0.4829629131445341, C1=0.8365163037378079,
     *     C2=0.2241438680420134, C3=-0.1294095225512604)
c     **
c     ** Local scalars 
      INTEGER        NH, NH1, I, J
c     **
c     ** Local arrays
c     **
c     ** External functions 
c     **
c     ** External subroutines 
c     **
c     ** Intrinsic functions 
c     .. Start program 
      IF(N.LT.4)RETURN
      NH=N/2
      NH1=NH+1
      IF (ISIGN.GE.0) THEN      !apply filter.
      I=1
      DO J=1,N-3,2
         WKSP(I)=C0*A(J)+C1*A(J+1)+C2*A(J+2)+C3*A(J+3)
         WKSP(I+NH)=C3*A(J)-C2*A(J+1)+C1*A(J+2)-C0*A(J+3)
         I=I+1
      ENDDO 
      WKSP(I)=C0*A(N-1)+C1*A(N)+C2*A(1)+C3*A(2)
      WKSP(I+NH)=C3*A(N-1)-C2*A(N)+C1*A(1)-C0*A(2)
      ELSE                      !apply transpose filter.
         WKSP(1)=C2*A(NH)+C1*A(N)+C0*A(1)+C3*A(NH1)
         WKSP(2)=C3*A(NH)-C0*A(N)+C1*A(1)-C2*A(NH1)
         J=3
         DO I=1,NH-1
            WKSP(J)=C2*A(I)+C1*A(I+NH)+C0*A(I+1)+C3*A(I+NH1)
            WKSP(J+1)=C3*A(I)-C0*A(I+NH)+C1*A(I+1)-C2*A(I+NH1)
            J=J+2
         ENDDO
      ENDIF
      DO I=1,N
         A(I)=WKSP(I)
      ENDDO 
c     ..
      RETURN
      END
