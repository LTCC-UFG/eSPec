      SUBROUTINE FFTNT( X, Y, N, M, ITYPE )
      IMPLICIT NONE
c
c ... Scalar arguments ...
      INTEGER  N, M, ITYPE
c ... Array arguments ...
      REAL*8     X(*), Y(*)
c
c  This is a very simple Cooley-Tukey Radix-2 DIF Complex FFT.
c 
c  The data sequence is of size N = 2**M.
c  X and Y contain the real and imaginary parts of the data.
c
c  ITYPE .ne. -1 for forward transform
c  ITYPE .eq. -1 for backward transform
c
c  The forward transform computes
c     Z(k) = sum_{j=0}^{N-1} z(j)*exp(-2ijk*pi/N)
c
c  The backward transform computes
c     z(j) = (1/N) * sum_{k=0}^{N-1} Z(k)*exp(2ijk*pi/N)
c
c
c  Steve Kifowit, 31 October 1998
c       
c ... Local scalars ...
      INTEGER    I, J, K, L, N1, N2, IE, IA
      REAL*8     C, S, XT, YT, P, TWOPI, A
c ... Parameters ...
      PARAMETER  (TWOPI = 6.283185307179586476925287)
c ... Intrinsic functions
      INTRINSIC  SIN, COS
c
c ... Exe. statements ... 
c
c ... Quick return ...
      IF ( N .EQ. 1 ) RETURN
c
c ... Conjugate if necessary ...
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 1, I = 1, N
	    Y(I) = - Y(I)
 1       CONTINUE
      ENDIF
c
c ... Main loop ...
      P = TWOPI / N
      N2 = N
      DO 10, K = 1, M
         N1 = N2
         N2 = N2 / 2
         IE = N / N1
         IA = 1
         DO 20, J = 1, N2
	    A = ( IA - 1 ) * P
            C = COS( A ) 
            S = SIN( A )
            IA = IA + IE
            DO 30, I = J, N, N1
               L = I + N2
               XT = X(I) - X(L)
               X(I) = X(I) + X(L)
               YT = Y(I) - Y(L)
               Y(I) = Y(I) + Y(L)
               X(L) = C * XT + S * YT
               Y(L) = C * YT - S * XT
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
c
c ... Bit reversal permutation ...
 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
         IF ( I .GE. J ) GOTO 101
         XT = X(J)
         X(J) = X(I)
         X(I) = XT
         YT = Y(J)
         Y(J) = Y(I)
         Y(I) = YT
 101     K = N / 2
 102     IF ( K .GE. J ) GOTO 103
         J = J - K
         K = K / 2
         GOTO 102
 103     J = J + K
 104  CONTINUE
c
c ... Conjugate and normalize if necessary ...
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 3, I = 1, N
	    Y(I) = - Y(I) / N
            X(I) = X(I) / N
 3       CONTINUE
      ENDIF
      RETURN
c
c ... End of subroutine CFFT ...
c
      END
