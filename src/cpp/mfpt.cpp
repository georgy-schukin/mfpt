#include "array2d.h"

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

using DArray2 = Array2D<double>;
using DArray1 = std::vector<double>;

void pr2(DArray2 &b, int l1, int m1, int l2, int m2, ofstream &out) {
    vector<double> v(8);
    vector<int> mv(8);

    int n1 = l1;
    int n2 = l1 + 7;

    while (true) {
        if (n2 > m1) {
            n2 = m1;
        }
        int n3 = n2 - n1 + 1;
        for (int i = n1; i <= n2; i++) {
            mv[i + 1 - n1 - 1] = i;
        }
        //write(25,902) (mv(i),i=1,n3)
        //902 format(7x,8(i3,6x))
        out << setw(7) << "";
        for (int i = 0; i < n3; i++) {
            out << setw(3) << mv[i] << setw(6) << "";
        }
        out << endl;
        for (int l0 = l2; l0 <= m2; l0++) {
            int l = m2 + l2 - l0;
            int l02 = l;
            for (int i = n1; i <= n2; i++) {
                v[i - n1] = b(i - 1, l - 1);
            }
            //write(25,903) l02,(v(i),i=1,n3)
            //903 format(i3,1x,8f9.3)
            out << setw(3) << l02 << setw(1) << "";
            for (int i = 0; i < n3; i++) {
                out << setw(9) << std::fixed << setprecision(3) << v[i];
            }
            out << endl;
        }
        if (n2 == m1) {
            return;
        }
        n1 = n1 + 8;
        n2 = n2 + 8;
    }
}

void output(const string &header, DArray2 &b, int l1, int m1, int l2, int m2, ofstream &out) {
    out << "\n";
    out << header << "\n";
    pr2(b, l1, m1, l2, m2, out);
}

void output(const string &header, DArray2 &array, const std::array<int, 4> &range, ofstream &out) {
    out << "\n";
    out << header << "\n";
    pr2(array, range[0], range[1], range[3], range[4], out);
}

int main() {
/*
    program brbz003c
    Тестовая программа.
    Двумерное уравнение Пуассона. Гран.условия 2-го рода.
    Преобразование Фурье и прогонка.
    Предварительно вычисленные синусы.
*/

    const size_t imp = 42;
    const size_t kmp = 122;

/*
    real*8 br(imp,kmp),bf(imp,kmp),bz(imp,kmp)
    real*8 sb(2*imp),jb(2*imp),al(2*imp),be(2*imp)
    real*8 aa(2*imp,kmp),jf(2*imp,kmp),aa1(2*imp,kmp)
    real*8 gg(2*imp,kmp),bb(2*imp,kmp),ff(2*imp,kmp),phi(2*imp,kmp)
    real*8 dd(imp,kmp),phi1(2*imp,kmp),ff1(2*imp,kmp),ds(2*kmp)
*/

    DArray2 br(imp, kmp), bf(imp, kmp), bz(imp, kmp);
    DArray1 sb(2 * imp), jb(2 * imp), al(2 * imp), be(2 * imp);
    DArray2 aa(2 * imp, kmp), jf(2 * imp, kmp), aa1(2 * imp, kmp);
    DArray2 gg(2 * imp, kmp), bb(2 * imp, kmp), ff(2 * imp, kmp), phi(2 * imp, kmp);
    DArray2 dd(imp, kmp), phi1(2 * imp, kmp), ff1(2 * imp, kmp);
    DArray1 ds(2 * kmp);

/*
      open(25,file='brbz003c.lst',form='formatted')
      open(16,file='aa11.txt',form='formatted')
c      open(17,file='conv11.dat',form='formatted')
c      open(18,file='ds.txt',form='formatted')
*/

    ofstream out_25("brbz003c.lst");
    ofstream out_16("aa11.txt");

    const int im = 20;
    const int km = 60;
    const double pi = 3.14159265358979;
    double c = pi / km;
    const double hr = 0.2;
    const double hz = 0.2;
    const double rm = im * hr;
    const double zm = km * hz;
    const double hr2 = hr * hr;
    const double hz2 = hz * hz;

    const std::array<int, 4> output_range = {1, 7, 1, 6};

/*
    тестовое решение

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
*/

    double a0 = -0.1;
    double a = 1.0;
    double d = 1.0;
    for (int k = 1; k <= km + 2; k++) {
        double z = hz * (k - 1.5);
        double s = a * z * z * (z - 1.5 * zm) + d;
        s = a0 * z * z * (z * z - 2.0 * zm * zm) + a * z * z * (z - 1.5 * zm) + d;
        for (int i = 2; i <= 2 * im + 2; i++) {
            aa1(i - 1, k - 1) = s * (hr * (i - 1.5) * (2.0 * rm - hr * (i - 2.0)));
        }
        aa1(0, k - 1) = -aa1(1, k - 1);
    }

/*
    тестовые токи

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
*/

    for (int k = 2; k <= km + 1; k++) {
        double s = (1.5 * aa1(2, k - 1) - 4.5 * aa1(1, k - 1)) / hr2 +
                   (aa1(1, k) - 2.0 * aa1(1, k - 1) + aa1(1, k - 2)) / hz2;
        jf(1, k - 1) = -s;
        for (int i = 3; i < 2 * im + 1; i++) {
            s = (((i - 0.5) * aa1(i + 1, k) - (i - 1.5) * aa1(i, k)) / (i - 1.0) -
                 ((i - 1.5) * aa1(i, k) - (i - 2.5) * aa1(i - 1, k)) / (i - 2.0)) / hr2 +
                (aa1(i, k + 1) - 2.0 * aa1(i, k) + aa1(i, k - 1)) / hz2;
            jf(i, k) = -s;
        }
    }

    output("aa1 aa1", aa1, output_range, out_25);
    output("jf jf", jf, output_range, out_25);

/*
    вычисление разностей

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
*/
    for (int i = 2; i < 2 * im + 1; i++) {
        for (int k = 2; k < km; k++) {
            gg(i, k) = jf(i, k + 1) - jf(i, k);
            phi1(i, k) = aa1(i, k + 1) - aa1(i, k);
        }
        gg(i, 1) = 0.0;
        gg(i, km + 1) = 0.0;
        phi1(i, 1) = 0.0;
        phi1(i, km + 1) = 0.0;
    }

    output("gg gg", gg, output_range, out_25);

/*
    вычисление синусов

      do k=1,2*km
         ds(k)=dsin(c*k)
      enddo
*/
    for (int k = 1; k < 2 * km; k++) {
        ds[k] = sin(c * k);
    }

/*
    преобразование Фурье для правых частей и для решения

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
*/
    auto km2 = km / 2;
    for (int i = 2; i < 2 * im + 1; i++) {
        for (int j = 2; j < km; j++) {
            double s1 = 0.0, s2 = 0.0;
            int k1 = 0;
            for (int k = 2; k < km; k++) {
                k1 = k1 + j - 1;
                if (k1 > 2 * km) {
                    k1 = k1 - 2 * km;
                }
                s1 = s1 + gg(i, k) * ds[k1];
                s2 = s2 + phi1(i, k) * ds[k1];
            }
            bb(i, j) = s1 / km2;
            ff1(i, j) = s2 / km2;
        }
        bb(i, 1) = 0.0;
        bb(i, km + 1) = 0.0;
        ff1(i, 1) = 0.0;
        ff1(i, km + 1) = 0.0;
    }

    output("bb bb", bb, 1, 7, 1, 6, out_25);

/*
    прогонка по радиусу

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
*/
    for (int j = 2; j < km; j++) {
        double s = 9.0 / (2.0 * hr * hr) + (4.0 / hz2) * (sin(c * (j - 1.0) / 2.0)) * (sin(c * (j - 1.0) / 2.0));
        al[2] = 3.0 / (2.0 * hr * hr * s);
        be[2] = bb(2, j) / s;
        for (int i = 3; i < 2 * im + 1; i++) {
            const auto dsin = sin(c * (j - 1.0) / 2.0);
            s = (2.0 * ((i - 1.5) / hr) * ((i - 1.5) / hr)) / ((i - 1.0) * (i - 2.0)) +
                (4.0 / hz2) * dsin * dsin -
                al[i - 1] * (i - 2.5) / ((i - 2.0) * hr * hr);
            al[i] = (i - 0.5) / (s * (i - 1.0) * hr * hr);
            be[i] = (be[i - 1] * (i - 2.5) / ((i - 2.0) * hr * hr) + bb(i, j)) / s;
        }
        ff(2 * im + 2, j) = 0.0;
        for (int i = 2 * im + 1; i > 2; i--) {
            ff(i, j) = al[i] * ff(i + 1, j) + be[i];
        }
    }

    output("ff ff", ff, output_range, out_25);
    output("ff1 ff1", ff1, output_range, out_25);

/*
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
*/

/*
    обратное преобразование Фурье

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
*/
    for (int i = 2; i < 2 * im + 1; i++) {
        for (int k = 2; k < km; k++) {
            double s1 = 0.0;
            int k1 = 0;
            for (int j = 2; j < km; j++) {
                k1 = k1 + k - 1;
                if (k1 > 2 * km) {
                    k1 = k1 - 2 * km;
                }
                s1 = s1 + ff(i, j) * ds[k1];
            }
            phi(i, k) = s1;
        }
        phi(i, 1) = 0.0;
        phi(i, km + 1) = 0.0;
    }

    output("phi phi", phi, output_range, out_25);
    output("phi1 phi1", phi1, output_range, out_25);

/*
    proverka2 решения dd dd

      do i=3,im+1
         do k=2,km
           dd(i,k)=(((i-0.5d0)*phi(i+1,k)-(i-1.5d0)*phi(i,k))/(i-1.d0)-
     =       ((i-1.5d0)*phi(i,k)-(i-2.5d0)*phi(i-1,k))/(i-2.d0))/hr**2+
     =       (phi(i,k+1)-2.d0*phi(i,k)+phi(i,k-1))/hz**2+gg(i,k)
         enddo
      enddo
*/

    for (int i = 3; i < im + 1; i++) {
        for (int k = 2; k < km; k++) {
            dd(i, k) = (((i - 0.5) * phi(i + 1, k) - (i - 1.5) * phi(i, k)) / (i - 1.0) -
                        ((i - 1.5) * phi(i, k) - (i - 2.5) * phi(i - 1, k)) / (i - 2.0)) / hr2 +
                       (phi(i, k + 1) - 2.0 * phi(i, k) + phi(i, k - 1)) / hz2 + gg(i, k);
        }
    }

    output("proverka2 dd dd", dd, output_range, out_25);

/*
    решение при к=2

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
*/
    al[1] = 1.0 / 3.0;
    be[1] = (2.0 * hr * hr / 9.0) * (jf(2, 2) + phi(2,2) / hz2);

    for (int i = 2; i < 2 * im; i++) {
        double s = 2.0 * (i - 0.5) * (i - 0.5) / (i * (i - 1.0)) - al[i - 1] * (i - 1.5) / (i - 1.0);
        al[i] = (i + 0.5) / (i * s);
        be[i] = (be[i - 1] * (i - 1.5) / (i - 1.0) + hr * hr * (jf(i + 1, 2) + phi(i + 1, 2) / hz2)) / s;
    }

    aa(2 * im + 2, 2) = 0.0;
    for (int i = 2 * im; i > 1; i--) {
        aa(i + 1, 2) = al[i] * aa(i + 2, 2) + be[i];
    }

/*
    решение во всей области

      do i=2,im+2
         aa(i,1)=aa(i,2)
         do k=2,km+1
            aa(i,k+1)=aa(i,k)+phi(i,k)
         enddo
      enddo
*/

    for (int i = 2; i < im + 2; i++) {
        aa(i, 1) = aa(i, 2);
        for (int k = 2; k < km + 1; k++) {
            aa(i, k + 1) = aa(i, k) + phi(i, k);
        }
    }

/*
    proverka3 решения dd dd

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
*/

    for (int i = 3; i < im + 1; i++) {
        for (int k = 2; k < km + 1; k++) {
            dd(i, k) = (((i - 0.5) * aa(i + 1, k) - (i - 1.5) * aa(i, k)) / (i - 1.0) -
                        ((i - 1.5) * aa(i, k) - (i - 2.5) * aa(i - 1, k)) / (i - 2.0)) / hr2 +
                       (aa(i, k + 1) - 2.0 * aa(i, k) + aa(i, k - 1)) / hz2 + jf(i, k);

        }
    }

    output("dd dd", dd, output_range, out_25);

/*
    вычисление магнитных полей

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
*/
    for (int k = 1; k < km + 2; k++) {
        bz(1, k) = 4.0 * aa(2, k) / hr;
        for (int i = 2; i < im + 2; i++) {
            bz(i, k) = ((i - 0.5) * aa(i + 1, k) - (i - 1.5) * aa(i, k)) / (hr * (i - 1.0));
        }
    }

    for (int k = 1; k < km + 1; k++) {
        for (int i = 1; i < im + 2; i++) {
            br(i, k) = -(aa(i, k + 1) - aa(i, k)) / hz;
        }
    }

    return 0;
}
