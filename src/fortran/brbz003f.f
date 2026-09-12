      program brbz003f
c
c     тестовая программа
c     двумерное уравнение Пуассона. Гран.условия 2-го рода.
c     Быстрое преобразование Фурье и прогонка.
c     
c
      implicit none
      integer imp,kmp,im,km,i,k,n,k1,km2,j,k2,i1,m,nj,j1
c      parameter(im=20,n=6,km=2**n,imp=im+2,kmp=km+2)
      include 'brz2.par'
      real*8 br(imp,kmp),bf(imp,kmp),bz(imp,kmp)
      real*8 al(2*imp),be(2*imp)
      real*8 aa(2*imp,kmp),jf(2*imp,kmp),aa1(2*imp,kmp)
      real*8 gg(2*imp,kmp),phi(2*imp,kmp),phi2(2*imp,kmp)
      real*8 dd(imp,kmp),phi1(2*imp,kmp)
      real*8 hr,hz,s,s1,s2,s3,eps,b0,pi,c,s4,t1,t2
c      real*8 B_sol,minim,maxim
c      real*8 r_coil, z_coil,coil_length,po
c      real*8 z_coil_left, z_coil_right
      real*8 alpha, delta,r,f,rm,zm,a,d,z,a0
      complex*16 dan(2*km)
      real*8 bb(2*imp,4*km),ff(2*imp,4*km),ff1(2*imp,4*km)

      common/a/hr,hz
c      common/c/im,km
c      common/g/br,bf,bz
      common/jf/jf

      open(25,file='brbz003f.lst',form='formatted')
c      open(16,file='aa11.txt',form='formatted') 
c      open(17,file='conv11.dat',form='formatted') 
c      open(18,file='ds.txt',form='formatted') 

c      im=20
c      km=60
      pi=3.14159265358979d0
      c=0.5d0*pi/km
c      z_coil=0.d0
c      coil_length=1.d0
c      r_coil=4.3d0
c      eps=1.d-10
      rm=4.d0
      zm=12.d0
      hr=rm/im
      hz=zm/km
c      po=2.d0
c      b0=20.d0
      write(25,*) 'im,km=',im,km
      write(25,*)

      call cpu_time(t1)

      do m=1,1000
c      do m=1,1
c     тестовое решение
c
      a0=-0.1d0
      a=1.d-3
      d=1.d0+1.d-5*m

      do k=1,km+2
         z=hz*(k-1.5d0)
         s=a*z**2*(z-1.5d0-zm)+d
         do i=2,2*im+2
            aa1(i,k)=s*(i-2*im-2.d0)*(i-1.5d0)*hr**2
         enddo
            aa1(1,k)=-aa1(2,k)
            aa1(2*im+2,k)=0.d0
      enddo
         do i=1,2*im+2
            aa1(i,1)=aa1(i,2)
            aa1(i,km+2)=aa1(i,km+1)
         enddo
c      write(25,*)
c      write(25,*) 'aa1 aa1'
c      call pr21(aa1,1,7,1,km+2)

c     тестовые токи
c
      do k=2,km+1
         s=(1.5d0*aa1(3,k)-4.5d0*aa1(2,k))/hr**2+
     =        (aa1(2,k+1)-2.d0*aa1(2,k)+aa1(2,k-1))/hz**2
         jf(2,k)=-s
         do i=3,2*im+1
            s=(((i-0.5d0)*aa1(i+1,k)-(i-1.5d0)*aa1(i,k))/(i-1.d0)-
     =       ((i-1.5d0)*aa1(i,k)-(i-2.5d0)*aa1(i-1,k))/(i-2.d0))/hr**2+
     =        (aa1(i,k+1)-2.d0*aa1(i,k)+aa1(i,k-1))/hz**2
            jf(i,k)=-s
         enddo
      enddo

c      write(25,*)
c      write(25,*) 'jf jf'
c      call pr21(jf,1,7,1,6)

c=======================================================================
c
c     преобразование Фурье для правых частей 
c
      do i=2,2*im+1

         do k=1,km
            dan(k)=dcmplx(jf(i,k+1),0.d0)
            dan(k+km)=dcmplx(jf(i,km+2-k),0.d0)
         enddo

         call fftc(dan,2*n,1,2*km)

         do k=1,2*km
            bb(i,k)=dreal(dan(k))
            bb(i,2*km+k)=dimag(dan(k))
         enddo
      enddo    ! i

c      write(25,*)
c      write(25,*) 'bb bb'
c      call pr22(bb,3,3,1,4*km)
c
c     прогонка по радиусу ===============================
c
c            write(16,*)
c            write(16,*) 'progonka'
c            write(16,*)

      do j=1,4*km
         j1=j-1
         if(j1.ge.2*km) j1=j1-2*km

c      s1=dsin(c*j1)
c      write(16,*) j,j1,s1

         s=9.d0/(2.d0*hr**2)+(4.d0/hz**2)*(dsin(c*j1))**2
         al(2)=3.d0/(2.d0*hr**2*s)
         be(2)=bb(2,j)/s

         do i=3,2*im+1
            s=(2.d0*((i-1.5d0)/hr)**2)/((i-1.d0)*(i-2.d0))+
     =         (4.d0/hz**2)*(dsin(c*j1))**2-
     =         al(i-1)*(i-2.5d0)/((i-2.d0)*hr**2)
            al(i)=(i-0.5d0)/(s*(i-1.d0)*hr**2)
            be(i)=(be(i-1)*(i-2.5d0)/((i-2.d0)*hr**2)+bb(i,j))/s
         enddo

         ff(2*im+2,j)=0.d0

         do i=2*im+1,2,-1
            ff(i,j)=al(i)*ff(i+1,j)+be(i)
         enddo

      enddo     !   j

c      write(25,*)
c      write(25,*) 'ff ff'
c      call pr22(ff,1,7,1,4*km)



c-------------------------- konets proverka 1 решения dd dd
c
c     обратное преобразование Фурье 
c
      do i=2,2*im+1

         do k=1,2*km
            dan(k)=dcmplx(ff(i,k),ff(i,k+2*km))
         enddo

         call fftc(dan,2*n,-1,2*km)

         do k=1,km
            phi(i,k+1)=dreal(dan(k))
         enddo
         phi(i,1)=phi(i,2)
         phi(i,km+2)=phi(i,km+1)

c         do k=1,2*km
c            ff1(i,k)=dreal(dan(k))
c            ff1(i,k+2*km)=dimag(dan(k))
c         enddo

      enddo    ! i

      enddo    ! m

      call cpu_time(t2)

      write(25,*) 't1:',t1 
      write(25,*) 't2:',t2 
      write(25,*) 't2-t1=',t2-t1


c      write(25,*)
c      write(25,*) 'ff1 ff1'
c      call pr22(ff1,1,7,1,4*km)

c      write(25,*)
c      write(25,*) 'aa1 aa1'
c      call pr21(aa1,1,7,1,6)
      write(25,*)
      write(25,*) 'phi phi'
      call pr21(phi,1,7,1,km+2)

      do i=3,im+1
         do k=2,km
c           dd(i,k)=(((i-0.5d0)*phi(i+1,k)-(i-1.5d0)*phi(i,k))/(i-1.d0)-
c     =       ((i-1.5d0)*phi(i,k)-(i-2.5d0)*phi(i-1,k))/(i-2.d0))/hr**2+
c     =       (phi(i,k+1)-2.d0*phi(i,k)+phi(i,k-1))/hz**2+jf(i,k)
            dd(i,k)=phi(i,k)-aa1(i,k)
         enddo
      enddo

      write(25,*)
      write(25,*) 'proverka 2 dd dd'
      call pr2(dd,1,7,1,6)

c-------------------- proverka 3 решения dd dd

c      write(25,*)
c      write(25,*) 'proverka 3 dd dd'
c      s=0.d0
c      do k=2,km+1
c         do i=2,2*im+1
c            s1=dabs(phi(i,k)-aa1(i,k))
c            if(s1.gt.s) then
c               s=s1
c               i1=i
c               k1=k
c            endif
c         enddo
c      enddo
c      write(25,*)
c      write(25,*) 'max|phi-aa1|=',s,i1,k1

c      write(25,*)
c      write(25,*) 'aa aa'
c      call pr21(phi,1,7,1,6)

c      do i=3,im+1
c         do k=2,km+1
c          dd(i,k)=(((i-0.5d0)*phi(i+1,k)-(i-1.5d0)*phi(i,k))/(i-1.d0)-
c     =      ((i-1.5d0)*phi(i,k)-(i-2.5d0)*phi(i-1,k))/(i-2.d0))/hr**2+
c     =        (phi(i,k+1)-2.d0*phi(i,k)+phi(i,k-1))/hz**2+jf(i,k)
c         enddo
c      enddo

c      write(25,*)
c      write(25,*) 'dd dd'
c      call pr2(dd,1,7,1,6)

      end
c=====================================================
      subroutine fftc(a,n,isi,np)
      implicit real*8(a-h,o-z)
      complex*16 a(np),t,w,w1
      do 7 i=1,np
      if(isi.gt.0) a(i)=a(i)/np
    7 continue
      pi2=8.d-0*datan(1.d-0)
      nn=np
      j=1
      do 3 i=1,nn
      if(i.ge.j) go to 1
      t=a(j)
      a(j)=a(i)
      a(i)=t
    1 m=nn/2
    2 if(j.le.m) go to 3
      j=j-m
      m=m/2
      if(m.ge.1) go to 2
    3 j=j+m
      mm=1
    4 if(mm.ge.nn) return
      ii=2*mm
c      th= pi2/dsign(ii,n*isi)
      th= pi2/isign(ii,n*isi)
      w1=dcmplx(-2.0 d-0*dsin(th/2)**2,dsin(th))
      w=1
      do 6 m=1,mm
      do 5 i=m,nn,ii
      t=w*a(i+mm)
      a(i+mm)=a(i)-t
      a(i)=a(i)+t
    5 continue
      w=w1*w+w
    6 continue
      mm=ii
      go to 4
      end
c===================================================
      subroutine pr2(b,l1,m1,l2,m2)
      include 'brz2.par'
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
      include 'brz2.par'
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
c======================================================
      subroutine pr22(b,l1,m1,l2,m2)
      include 'brz2.par'
c      parameter(im2=42,km2=122)
      real*8 b(2*imp,4*km),v(8)
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