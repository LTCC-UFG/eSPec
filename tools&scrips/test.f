      implicit real*8(a-h,o-z)
c
      xi = -5
      xf =  10
      dx = 0.0001
      QC = 1.15625/5.291772083D-1
      DELQ = 0.052918/5.291772083D-1
c      QC = 0d0
      DELQ = 2d0
      do x=xi,xf,dx
         xd=x
         t1 = EXP(-((XD - QC)/DELQ)**2)
         t2 = (EXP(-((XD + dx - QC)/DELQ)**2) 
     &           - EXP(-((XD - dx - QC)/DELQ)**2))/(2*dx)
         t3 = -2*(XD - QC)/DELQ**2*t1
         write(1,*)x,t1,t2,t3
      enddo
c
      end
