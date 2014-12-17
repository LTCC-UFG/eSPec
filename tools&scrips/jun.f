      implicit real*8(a-h,o-z)
c
      character*20 lx1
c
      do i=1,10000,1
         write(*,*)'i',i
         read(1,*,end=10)lx1, lx1, lx1, dt1, dt2, dt3
          write(*,*)'i',i
         write(2,*)dt1, dt3
      enddo
 10   continue
c
      end
