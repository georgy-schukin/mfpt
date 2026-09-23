#include "defs.h"
#include "../common/output.h"
#include "../common/timer.h"

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <complex>

using namespace std;

/*
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
 */

using Complex = std::complex<double>;
using CArray1 = std::vector<Complex>;
using namespace std::complex_literals;

void fftc(CArray1 &a, int n, int isi, int np) {
    if (isi > 0) {
        for (int i = 0; i < np; i++) {
            a[i] /= np;
        }
    }
    const double pi2 = 8.0 * std::atan(1.0);
    int j = 0;
    for (int i = 0; i < np; i++) {
        if (i < j) {
            std::swap(a[i], a[j]);
        }
        int m = np / 2;
        do {
            if (j < m) {
                break;
            }
            j -= m;
            m /= 2;
        } while (m >= 1);
        j += m;
    }
    int mm = 1;
    while (mm < np) {
        const int ii = 2 * mm;
        const double th = pi2 / (n * isi >= 0 ? std::abs(ii) : -std::abs(ii));
        const auto sin_th = std::sin(th);
        const auto sin_th_d2 = std::sin(th * 0.5);
        Complex w1 {-2.0 * sin_th_d2 * sin_th_d2, sin_th};
        Complex w = 1.0;
        for (int m = 0; m < mm; m++) {
            for (int i = m; i < np; i += ii) {
                const Complex t = w * a[i + mm];
                a[i + mm] = a[i] - t;
                a[i] += t;
            }
            w += w1 * w;
        }
        mm = ii;
    }
}

int main(int argc, char **argv) {
/*
    Двумерное уравнение Пуассона. Гран.условия 2-го рода.
    Быстрое преобразование Фурье и прогонка.
*/

    const int IM_DEF = 80;
    const int N_DEF = 7;
    const int FOUT_DEF = 1;
    const int FULL_OUTPUT_DEF = 0;

    if (argc > 1) {
        const auto s = string(argv[1]);
        if (s == "-h" || s == "--help") {
            cout << argv[0] <<
                " [im=" << IM_DEF << "]" <<
                " [n=" << N_DEF << "]" <<
                " [file_output=" << FOUT_DEF << "]" <<
                " [full_output=" << FULL_OUTPUT_DEF << "]" <<
                endl;
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int n = (argc > 2) ? stoi(argv[2]) : N_DEF;
    const int km = std::pow(2, n);
    const bool file_output = (argc > 3) ? stoi(argv[3]) : FOUT_DEF;
    const bool full_output = (argc > 4) ? stoi(argv[4]) : FULL_OUTPUT_DEF;

/*
    integer imp,kmp,im,km,i,k,n,k1,km2,j,k2,i1,m,nj,j1
    parameter(im=20,n=6,km=2**n,imp=im+2,kmp=km+2)
    real*8 br(imp,kmp),bf(imp,kmp),bz(imp,kmp)
    real*8 al(2*imp),be(2*imp)
    real*8 aa(2*imp,kmp),jf(2*imp,kmp),aa1(2*imp,kmp)
    real*8 gg(2*imp,kmp),phi(2*imp,kmp),phi2(2*imp,kmp)
    real*8 dd(imp,kmp),phi1(2*imp,kmp)
    real*8 hr,hz,s,s1,s2,s3,eps,b0,pi,c,s4,t1,t2
    real*8 alpha, delta,r,f,rm,zm,a,d,z,a0
    complex*16 dan(2*km)
    real*8 bb(2*imp,4*km),ff(2*imp,4*km),ff1(2*imp,4*km)
*/

    const size_t imp = im + 2;
    const size_t imp2 = 2 * imp;
    const size_t kmp = km + 2;

    //DArray2 br(ims, kms), bf(ims, kms), bz(ims, kms);
    DArray2 aa(imp2, kmp), jf(imp2, kmp), aa1(imp2, kmp);
    DArray2 gg(imp2, kmp), bb(imp2, 4 * km), ff(imp2, 4 * km), phi(imp2, kmp);
    //DArray2 dd(ims, kms), phi1(ims2, kms), ff1(ims2, kms);

/*
    pi=3.14159265358979d0
    c=0.5d0*pi/km
    rm=4.d0
    zm=12.d0
    hr=rm/im
    hz=zm/km
*/

    const double pi = 3.14159265358979;
    const double c = 0.5 * pi / km;
    const double rm = 4.0;
    const double zm = 12.0;
    const double hr = rm / im;
    const double hz = zm / km;
    const double hr2 = hr * hr;
    const double hz2 = hz * hz;
    const double rhr2 = 1.0 / hr2;
    const double rhz2 = 1.0 / hz2;

    const std::array<int, 4> output_range = {0, 7, 0, 6};

    ofstream out_lst;
    if (file_output) {
        out_lst.open("output2.lst");
    }

    auto doOutput = [&](const std::string &header, const DArray2 &data) {
        if (file_output) {
            if (full_output) {
                output(header, data, out_lst);
            } else {
                output(header, data, output_range, out_lst);
            }
        }
    };

    auto doOutputRange = [&](const std::string &header, const DArray2 &data, const std::array<int, 4> &range) {
        if (file_output) {
            if (full_output) {
                output(header, data, out_lst);
            } else {
                output(header, data, range, out_lst);
            }
        }
    };

    double ft_time = 0;
    double prog_time = 0;
    //double sol_time = 0;
    //double mag_time = 0;

    // Init functions

/*
    тестовое решение

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
*/

    auto initTestSolution = [im, km, zm, hz, hr2](DArray2 &output, int m) {
        const double a0 = -0.1;
        const double a = 1e-3;
        const double d = 1.0 + 1e-5 * m;
        for (int k = 0; k < km + 2; k++) {
            const double z = hz * (k + 1 - 1.5);
            const double z2 = z * z;
            const double s = a * z2 * (z - 1.5 - zm) + d;
            for (int i = 1; i < 2 * im + 2; i++) {
                output(i, k) = s * (i + 1 - 2 * im - 2.0) * (i + 1 - 1.5) * hr2;
            }
            output(0, k) = -output(1, k);
            output(2 * im + 1, k) = 0.0;
        }
        for (int i = 0; i < 2 * im + 2; i++) {
            output(i, 0) = output(i, 1);
            output(i, km + 1) = output(i, km);
        }
    };

/*
    тестовые токи

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
*/

    auto getSolution = [rhr2, rhz2](const DArray2 &input, int i, int k) {
        return (((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (i) -
                ((i - 0.5) * input(i, k) - (i - 1.5) * input(i - 1, k)) / (i - 1.0)) * rhr2 +
               (input(i, k + 1) - 2.0 * input(i, k) + input(i, k - 1)) * rhz2;
    };

    auto initTestCurrent = [im, km, rhr2, rhz2, &getSolution](const DArray2 &input, DArray2 &output) {
        for (int k = 1; k < km + 1; k++) {
            double s = (1.5 * input(2, k) - 4.5 * input(1, k)) * rhr2 +
                       (input(1, k + 1) - 2.0 * input(1, k) + input(1, k - 1)) * rhz2;
            output(1, k) = -s;
            for (int i = 2; i < 2 * im + 1; i++) {
                output(i, k) = -getSolution(input, i, k);
            }
        }
    };

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

    auto computeDifference = [&](DArray2 &input, DArray2 &output) {
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = 1; k < km; k++) {
                output(i, k) = input(i, k + 1) - input(i, k);
            }
            output(i, 0) = 0.0;
            output(i, km) = 0.0;
        }
    };

/*
    преобразование Фурье для правых частей

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
*/

    auto computeFFT = [im, km, n, &ft_time](const DArray2 &input, DArray2 &output) {
        Timer tm;
        CArray1 dan(2 * km);
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = 0; k < km; k++) {
                dan[k] = Complex {input(i, k + 1), 0.0};
                dan[k + km] = Complex {input(i, km - k), 0.0};
            }
            fftc(dan, 2 * n, 1, 2 * km);
            for (int k = 0; k < 2 * km; k++) {
                output(i, k) = dan[k].real();
                output(i, 2 * km + k) = dan[k].imag();
            }
        }
        ft_time += tm.time();
    };

/*
    прогонка по радиусу

    do j=1,4*km
        j1=j-1
        if(j1.ge.2*km) j1=j1-2*km

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
*/

    auto computeProgonka = [im, km, imp2, hr, hr2, hz, hz2, c, &prog_time](const DArray2 &input, DArray2 &output) {
        Timer tm;
        DArray1 al(imp2), be(imp2);
        for (int k = 0; k < 4 * km; k++) {
            const int k1 = (k >= 2 * km ? k - 2 * km : k);
            const double dsin = sin(c * k1);
            double s = 9.0 / (2.0 * hr2) + (4.0 / hz2) * dsin * dsin;
            al[1] = 3.0 / (2.0 * hr2 * s);
            be[1] = input(1, k) / s;
            for (int i = 2; i < 2 * im + 1; i++) {
                s = (2.0 * ((i - 0.5) / hr) * ((i - 0.5) / hr)) / ((i) * (i - 1.0)) +
                    (4.0 / hz2) * dsin * dsin -
                    al[i - 1] * (i - 1.5) / ((i - 1.0) * hr2);
                al[i] = (i + 0.5) / (s * (i) * hr2);
                be[i] = (be[i - 1] * (i - 1.5) / ((i - 1.0) * hr2) + input(i, k)) / s;
            }
            output(2 * im + 1, k) = 0.0;
            for (int i = 2 * im; i >= 1; i--) {
                output(i, k) = al[i] * output(i + 1, k) + be[i];
            }
        }
        prog_time += tm.time();
    };

/*
    обратное преобразование Фурье

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
      enddo    ! i
*/

    auto computeFFTInverse = [&](const DArray2 &input, DArray2 &output) {
        Timer tm;
        CArray1 dan(2 * km);
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = 0; k < 2 * km; k++) {
                dan[k] = Complex {input(i, k), input(i, k + 2 * km)};
            }
            fftc(dan, 2 * n, -1, 2 * km);
            for (int k = 0; k < km; k++) {
                output(i, k + 1) = dan[k].real();
            }
            output(i, 0) = output(i, 1);
            output(i, km + 1) = output(i, km);
        }
        ft_time += tm.time();
    };

/*
    proverka2 решения dd dd

      do i=3,im+1
         do k=2,km
            dd(i,k)=phi(i,k)-aa1(i,k)
         enddo
      enddo
*/

    auto computeSolutionDifference = [&](DArray2 &first, const DArray2 &second) -> DArray2 {
        DArray2 output(im + 2, km + 2);
        for (int i = 2; i < im + 1; i++) {
            for (int k = 1; k < km; k++) {
                output(i, k) = first(i, k) - second(i, k);
            }
        }
        return output;
    };

    // The main program's body

    Timer full_time;

    for (int m = 1; m <= 1000; m++) {
        initTestSolution(aa1, m);
        initTestCurrent(aa1, jf);
        //doOutput("aa1 aa1", aa1);
        //doOutput("jf jf", jf);

        //computeDifference(jf, gg);
        //computeDifference(aa1, phi1);
        //doOutput("gg gg", gg);

        computeFFT(jf, bb);
        //doOutput("bb bb", bb);

        computeProgonka(bb, ff);
        //doOutput("ff ff", ff);
        //doOutput("ff1 ff1", ff1);

        computeFFTInverse(ff, phi);
        //doOutput("phi phi", phi);
        //doOutput("phi1 phi1", phi1);
    }

    doOutputRange("phi phi", phi, {0, 7, 0, km + 2});
    doOutput("proverka 2 dd dd", computeSolutionDifference(phi, aa1));

    if (file_output) {
        out_lst.close();
    }

    auto time = full_time.time();
    cout << "Im: " << im << ", Km: " << km << endl;
    cout << "TIME: " << time << endl;
    cout << "FT: " << ft_time <<
        ", Prog: " << prog_time << endl;

    return 0;
}
