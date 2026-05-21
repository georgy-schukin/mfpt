      program brbz003c
c
c     тестовая программа
c     двумерное уравнение Пуассона. Гран.условия 2-го рода.
c     Преобразование Фурье и прогонка.
c     Предварительно вычисленные синусы
c
      implicit none
      integer imp,kmp,im,km,i,k,n,k1,km2,j,k2
      include 'brz.par'
      real*8 br(imp,kmp),bf(imp,kmp),bz(imp,kmp)
      real*8 sb(2*imp),jb(2*imp),al(2*imp),be(2*imp)
      real*8 aa(2*imp,kmp),jf(2*imp,kmp),aa1(2*imp,kmp)
      real*8 gg(2*imp,kmp),bb(2*imp,kmp),ff(2*imp,kmp),phi(2*imp,kmp)
      real*8 dd(imp,kmp),phi1(2*imp,kmp),ff1(2*imp,kmp),ds(2*kmp)
      real*8 hr,hz,s,s1,s2,s3,eps,b0,pi,c,s4
c      real*8 B_sol,minim,maxim
c      real*8 r_coil, z_coil,coil_length,po
c      real*8 z_coil_left, z_coil_right
      real*8 alpha, delta,r,f,rm,zm,a,d,z,a0

      common/a/hr,hz
c      common/c/im,km
c      common/g/br,bf,bz
      common/jf/jf

      open(25,file='brbz003c.lst',form='formatted')
      open(16,file='aa11.txt',form='formatted') 
c      open(17,file='conv11.dat',form='formatted') 
c      open(18,file='ds.txt',form='formatted') 

      im=20
      km=60
      pi=3.14159265358979d0
      c=pi/km
c      z_coil=0.d0
c      coil_length=1.d0
c      r_coil=4.3d0
c      eps=1.d-10
      hr=0.2d0
      hz=0.2d0
      rm=im*hr
      zm=km*hz
c      po=2.d0
c      b0=20.d0

c     тестовое решение
c
      a0=-0.1d0
      a=1.d0
      d=1.d0
      do k=1,km+2
         z=hz*(k-1.5d0)
         s=a*z**2*(z-1.5d0*zm)+d
         s=a0*z**2*(z**2-2.d0*zm**2)+a*z**2*(z-1.5d0*zm)+d
c         s=dcos(pi*z/zm)
         do i=2,2*im+2
            aa1(i,k)=s*(hr*(i-1.5d0)*(2.d0*rm-hr*(i-2.d0)))
         enddo
            aa1(1,k)=-aa1(2,k)
      enddo

c     тестовые токи
c
      do k=2,km+1
         s=(1.5d0*aa1(3,k)-4.5d0*aa1(2,k))/hr**2+
     =        (aa1(2,k+1)-2.d0*aa1(2,k)+aa1(2,k-1))/hz**2
         jf(2,k)=-s
         do i=3,2*im+1
            s=(((i-0.5d0)*aa1(i+1,k)-(i-1.5d0)*aa1(i,k))/(i-1.d0)-
     =        ((i-1.5d0)*aa1(i,k)-(i-2.5d0)*aa1(i-1,k))/(i-2.d0))/hr**2+
     =        (aa1(i,k+1)-2.d0*aa1(i,k)+aa1(i,k-1))/hz**2
            jf(i,k)=-s
         enddo
      enddo

      write(25,*)
      write(25,*) 'aa1 aa1'
      call pr21(aa1,1,7,1,6)
      write(25,*)
      write(25,*) 'jf jf'
      call pr21(jf,1,7,1,6)
c      call pr21(jf,1,7,1,km+1)

c=======================================================================

c     вычисление разностей
c
      do i=2,2*im+1
         do k=2,km
            gg(i,k)=jf(i,k+1)-jf(i,k)
            phi1(i,k)=aa1(i,k+1)-aa1(i,k)
         enddo
         gg(i,1)=0.d0               !?
         gg(i,km+1)=0.d0            !?
         phi1(i,1)=0.d0               !?
         phi1(i,km+1)=0.d0            !?
      enddo

      write(25,*)
      write(25,*) 'gg gg'
      call pr21(gg,1,7,1,6)

c (0)
c      c=pi/km
c
c     вычисление синусов
c
      do k=1,2*km
         ds(k)=dsin(c*k)
      enddo
c
c
c     преобразование Фурье для правых частей и для решения
c
      km2=km/2
      do i=2,2*im+1
         do j=2,km
            s1=0.d0
            s2=0.d0
            k1=0
            do k=2,km
               k1=k1+j-1
               if(k1.gt.2*km) k1=k1-2*km
               s1=s1+gg(i,k)*ds(k1)
               s2=s2+phi1(i,k)*ds(k1)
c               s1=s1+gg(i,k)*dsin(c*(j-1.d0)*(k-1.d0))
c               s2=s2+phi1(i,k)*dsin(c*(j-1.d0)*(k-1.d0))
c               k2=(j-1.d0)*(k-1.d0)
c               s4=dsin(c*(j-1.d0)*(k-1.d0))
c               write(18,100) j,k,k1,k2,ds(k1),s4
c  100         format('j,k,k1,k2-',3i4,i6,2e12.4)
            enddo
            bb(i,j)=s1/km2
            ff1(i,j)=s2/km2
         enddo
         bb(i,1)=0.d0
         bb(i,km+1)=0.d0
         ff1(i,1)=0.d0
         ff1(i,km+1)=0.d0


      enddo    ! i

      write(25,*)
      write(25,*) 'bb bb'
      call pr21(bb,1,7,1,6)

c
c     прогонка по радиусу
c

      do j=2,km
         s=9.d0/(2.d0*hr**2)+(4.d0/hz**2)*(dsin(c*(j-1.d0)/2.d0))**2
         al(2)=3.d0/(2.d0*hr**2*s)
         be(2)=bb(2,j)/s

         do i=3,2*im+1
            s=(2.d0*((i-1.5d0)/hr)**2)/((i-1.d0)*(i-2.d0))+
     =         (4.d0/hz**2)*(dsin(c*(j-1.d0)/2.d0))**2-
     =         al(i-1)*(i-2.5d0)/((i-2.d0)*hr**2)
            al(i)=(i-0.5d0)/(s*(i-1.d0)*hr**2)
            be(i)=(be(i-1)*(i-2.5d0)/((i-2.d0)*hr**2)+bb(i,j))/s
         enddo

         ff(2*im+2,j)=0.d0
         do i=2*im+1,2,-1
            ff(i,j)=al(i)*ff(i+1,j)+be(i)
         enddo

      enddo     !   j

      write(25,*)
      write(25,*) 'ff ff'
      call pr21(ff,1,7,1,6)
      write(25,*)
      write(25,*) 'ff1 ff1'
      call pr21(ff1,1,7,1,6)

c--------------------------proverka1 решения dd dd
c      do k=2,km
c         dd(2,k)=bb(2,k)+ff(3,k)*3.d0/(2.d0*hr**2)-
c     =ff(2,k)*(9.d0/(2.d0*hr**2)+(dsin(c*(k-1.d0)/2.d0)*2.d0/hz)**2)
c         do i=3,2*im+1
c            dd(i,k)=bb(i,k)+ff(i+1,k)*(i-0.5d0)/((i-1.d0)*hr**2)+
c     =         ff(i-1,k)*(i-2.5d0)/((i-2.d0)*hr**2)-
c     =         ff(i,k)*((dsin(c*(k-1.d0)/2.d0)*2.d0/hz)**2+
c     =         2.d0*(((i-1.5d0)/hr)**2)/((i-1.d0)*(i-2.d0)))
c         enddo
c      enddo     !   k
c      write(25,*)
c      write(25,*) 'proverka1 dd dd'
c      call pr2(dd,1,7,1,6)
c
c     обратное преобразование Фурье 
c

      do i=2,2*im+1     !    im+1 ?

         do k=2,km
            s1=0.d0
            k1=0
            do j=2,km
               k1=k1+k-1
               if(k1.gt.2*km) k1=k1-2*km
               s1=s1+ff(i,j)*ds(k1)
c               s1=s1+ff(i,j)*dsin(c*(j-1.d0)*(k-1.d0))
            enddo
            phi(i,k)=s1
         enddo
         phi(i,1)=0.d0
         phi(i,km+1)=0.d0
         
      enddo     !   i

      write(25,*)
      write(25,*) 'phi phi'
      call pr21(phi,1,7,1,6)
      write(25,*)
      write(25,*) 'phi1 phi1'
      call pr21(phi1,1,7,1,6)

c--------------------------proverka2 решения dd dd

      do i=3,im+1
         do k=2,km
           dd(i,k)=(((i-0.5d0)*phi(i+1,k)-(i-1.5d0)*phi(i,k))/(i-1.d0)-
     =       ((i-1.5d0)*phi(i,k)-(i-2.5d0)*phi(i-1,k))/(i-2.d0))/hr**2+
     =       (phi(i,k+1)-2.d0*phi(i,k)+phi(i,k-1))/hz**2+gg(i,k)
         enddo
      enddo

      write(25,*)
      write(25,*) 'proverka2 dd dd'
      call pr2(dd,1,7,1,6)
c--------
c------------ решение при к=2


      al(1)=1.d0/3.d0
      be(1)=(2.d0*hr**2/9.d0)*(jf(2,2)+phi(2,2)/hz**2)

      do i=2,2*im
         s=2.d0*(i-0.5d0)**2/(i*(i-1.d0))-al(i-1)*(i-1.5d0)/(i-1.d0)
         al(i)=(i+0.5d0)/(i*s)
         be(i)=(be(i-1)*(i-1.5d0)/(i-1.d0)+
     =      (hr**2)*(jf(i+1,2)+phi(i+1,2)/hz**2))/s
      enddo

      aa(2*im+2,2)=0.d0
      do i=2*im,1,-1
         aa(i+1,2)=al(i)*aa(i+2,2)+be(i)
      enddo

c------------ решение во всей области

      do i=2,im+2
         aa(i,1)=aa(i,2)
         do k=2,km+1
            aa(i,k+1)=aa(i,k)+phi(i,k)
         enddo
      enddo

c--------
c-------------------- proverka3 решения dd dd

      write(25,*)
      write(25,*) 'aa aa'
      call pr21(aa,1,7,1,6)

      do i=3,im+1
         do k=2,km+1
            dd(i,k)=(((i-0.5d0)*aa(i+1,k)-(i-1.5d0)*aa(i,k))/(i-1.d0)-
     =        ((i-1.5d0)*aa(i,k)-(i-2.5d0)*aa(i-1,k))/(i-2.d0))/hr**2+
     =        (aa(i,k+1)-2.d0*aa(i,k)+aa(i,k-1))/hz**2+jf(i,k)
         enddo
      enddo

      write(25,*)
      write(25,*) 'dd dd'
      call pr2(dd,1,7,1,6)

c-------- вычисление магнитных полей

      do k=1,km+2
         bz(1,k)=4.d0*aa(2,k)/hr
         do i=2,im+2
         bz(i,k)=((i-0.5d0)*aa(i+1,k)-(i-1.5d0)*aa(i,k))/(hr*(i-1.d0))
         enddo
      enddo

      do k=1,km+1
         do i=1,im+2
            br(i,k)=-(aa(i,k+1)-aa(i,k))/hz
         enddo
      enddo

      end
c===================================================
      subroutine pr2(b,l1,m1,l2,m2)
      include 'brz.par'
c      parameter(im2=42,km2=122)
      real*8 b(imp,kmp),v(8)
      integer mv(8)
      n1=l1
      n2=l1+7
    5 if(n2.gt.m1) n2=m1
      n3=n2-n1+1
      do 2 i=n1,n2
    2 mv(i+1-n1)=i
      write(25,902) (mv(i),i=1,n3)
c  902 format(7x,8(i3,9x))
  902 format(7x,8(i3,6x))
      do 4 l0=l2,m2
      l=m2+l2-l0
      l02=l
      do 3 i=n1,n2
    3 v(i+1-n1)=b(i,l)
      write(25,903) l02,(v(i),i=1,n3)
c  903 format(i3,1x,8e12.4)
  903 format(i3,1x,8f9.3)
    4 continue
      if(n2.eq.m1) goto 1
      n1=n1+8
      n2=n2+8
      goto 5
    1 continue
      return
      end
c======================================================
      subroutine pr21(b,l1,m1,l2,m2)
      include 'brz.par'
c      parameter(im2=42,km2=122)
      real*8 b(2*imp,kmp),v(8)
      integer mv(8)
      n1=l1
      n2=l1+7
    5 if(n2.gt.m1) n2=m1
      n3=n2-n1+1
      do 2 i=n1,n2
    2 mv(i+1-n1)=i
      write(25,902) (mv(i),i=1,n3)
c  902 format(7x,8(i3,9x))
  902 format(7x,8(i3,6x))
      do 4 l0=l2,m2
      l=m2+l2-l0
      l02=l
      do 3 i=n1,n2
    3 v(i+1-n1)=b(i,l)
      write(25,903) l02,(v(i),i=1,n3)
c  903 format(i3,1x,8e12.4)
  903 format(i3,1x,8f9.3)
    4 continue
      if(n2.eq.m1) goto 1
      n1=n1+8
      n2=n2+8
      goto 5
    1 continue
      return
      end