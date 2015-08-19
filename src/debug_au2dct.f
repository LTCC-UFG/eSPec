      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,SHM(3),US(K-NP1AUX-2),
     &     SHM(3)*US(K-NP1AUX-2)  
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,SHM(2),US(K-NP1AUX),
     &     SHM(2)*US(K-NP1AUX)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,-SHM(3),US(K-NP1AUX+2),
     &     -SHM(3)*US(K-NP1AUX+2)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,CF4,US(K-NP(1)-1),
     &     CF4*US(K-NP(1)-1) 
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,CF3,US(K-NP(1)),
     &     CF3*US(K-NP(1))
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K, - CF4,US(K-NP(1)+1),
     &     -CF4*US(K-NP(1)+1)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,SHM(1),US(K-2),
     &     SHM(1)*US(K-2)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,CF2,US(K-1),
     &     CF2*US(K-1)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,(CF1 + VPOT(K)),US(K),
     &     (CF1 + VPOT(K))*US(K) 
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,CF2,US(K+1),
     &     CF2*US(K+1)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,SHM(1),US(K+2),
     &     SHM(1)*US(K+2)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,- CF4,US(K+NP(1)-1),
     &     - CF4*US(K+NP(1)-1)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,CF3,US(K+NP(1)),
     &     CF3*US(K+NP(1))
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,CF4,US(K+NP(1)+1),
     &     CF4*US(K+NP(1)+1)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,- SHM(3),US(K+NP1AUX-2),
     &     - SHM(3)*US(K+NP1AUX-2)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,SHM(2),US(K+NP1AUX),
     &     SHM(2)*US(K+NP1AUX)
      write(4,'(I6,ES26.18,ES26.18,ES26.18)')K,SHM(3),US(K+NP1AUX+2),
     &     SHM(3)*US(K+NP1AUX+2)
      write(4,'(I6,ES26.18)')K,AU(K)
