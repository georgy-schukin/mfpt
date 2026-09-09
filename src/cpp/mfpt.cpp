#include "defs.h"
#include "output.h"
#include "timer.h"

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;

DArray1 computeSins(size_t size, double coeff) {
    DArray1 dsins(size);
    for (size_t k = 0; k < size; k++) {
        dsins[k] = sin(coeff * k);
    }
    return dsins;
}

int main(int argc, char **argv) {
/*
    program brbz003c
    Тестовая программа.
    Двумерное уравнение Пуассона. Гран.условия 2-го рода.
    Преобразование Фурье и прогонка.
    Предварительно вычисленные синусы.
*/

    const int IM_DEF = 20;
    const int KM_DEF = 60;
    const int FOUT_DEF = 1;
    const int FULL_OUTPUT_DEF = 0;

    if (argc > 1) {
        const auto s = string(argv[1]);
        if (s == "-h" || s == "--help") {
            cout << argv[0] <<
                " [im=" << IM_DEF << "]" <<
                " [km=" << KM_DEF << "]" <<
                " [file_output=" << FOUT_DEF << "]" <<
                " [full_output=" << FULL_OUTPUT_DEF << "]" <<
                endl;
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int km = (argc > 2) ? stoi(argv[2]) : KM_DEF;
    const bool file_output = (argc > 3) ? stoi(argv[3]) : FOUT_DEF;
    const bool full_output = (argc > 4) ? stoi(argv[4]) : FULL_OUTPUT_DEF;

/*
    real*8 br(imp,kmp),bf(imp,kmp),bz(imp,kmp)
    real*8 sb(2*imp),jb(2*imp),al(2*imp),be(2*imp)
    real*8 aa(2*imp,kmp),jf(2*imp,kmp),aa1(2*imp,kmp)
    real*8 gg(2*imp,kmp),bb(2*imp,kmp),ff(2*imp,kmp),phi(2*imp,kmp)
    real*8 dd(imp,kmp),phi1(2*imp,kmp),ff1(2*imp,kmp),ds(2*kmp)
*/

    // aa1(2 * im + 2, km + 2)
    // jf(2 * im + 2, km + 2)
    // gg(2 * im + 2, km + 1), phi1(2 * im + 2, km + 1)
    // ds(2 * km)
    // bb(2 * im + 1, km + 1), ff1(2 * im + 1, km + 1)
    // al(2 * im + 1), be(2 * im + 1)
    // ff(2 * im + 2, km)
    // phi(2 * im + 1, km + 1)
    // dd(im + 1, km + 1)
    // aa(2 * im + 2, km + 2)
    // bz(im + 2, km + 2)
    // br(im + 2, km + 1)

    const size_t ims = im + 2;
    const size_t ims2 = 2 * im + 2;
    const size_t kms = km + 2;

    DArray2 br(ims, kms), bf(ims, kms), bz(ims, kms);
    DArray1 sb(ims2), jb(ims2), al(ims2), be(ims2);
    DArray2 aa(ims2, kms), jf(ims2, kms), aa1(ims2, kms);
    DArray2 gg(ims2, kms), bb(ims2, kms), ff(ims2, kms), phi(ims2, kms);
    DArray2 dd(ims, kms), phi1(ims2, kms), ff1(ims2, kms);

    ofstream out_lst;
    if (file_output) {
        out_lst.open("output.lst");
    }

    const double pi = 3.14159265358979;
    const double c = pi / km;
    const double hr = 0.2;
    const double hz = 0.2;
    const double rm = im * hr;
    const double zm = km * hz;
    const double hr2 = hr * hr;
    const double hz2 = hz * hz;

    const std::array<int, 4> output_range = {0, 7, 0, 6};

    auto doOutput = [&](const std::string &header, const DArray2 &data) {
        if (file_output) {
            if (full_output) {
                output(header, data, out_lst);
            } else {
                output(header, data, output_range, out_lst);
            }
        }
    };

    Timer full_time;

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

    // k, i: aa1(i, k) <- expr
    double a0 = -0.1;
    double a = 1.0;
    double d = 1.0;
    for (int k = 0; k < km + 2; k++) {
        const double z = hz * (k + 1 - 1.5);
        const double z2 = z * z;
        //double s = a * z * z * (z - 1.5 * zm) + d;
        const double s = a0 * z2 * (z2 - 2.0 * zm * zm) + a * z2 * (z - 1.5 * zm) + d;
        for (int i = 1; i < 2 * im + 2; i++) {
            aa1(i, k) = s * (hr * (i + 1 - 1.5) * (2.0 * rm - hr * (i + 1 - 2.0)));
        }
        aa1(0, k) = -aa1(1, k);
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

    auto compSolution = [hr2, hz2](const DArray2 &phi, int i, int k) {
        return (((i + 0.5) * phi(i + 1, k) - (i - 0.5) * phi(i, k)) / (i) -
                ((i - 0.5) * phi(i, k) - (i - 1.5) * phi(i - 1, k)) / (i - 1.0)) / hr2 +
               (phi(i, k + 1) - 2.0 * phi(i, k) + phi(i, k - 1)) / hz2;
    };

    auto compCheck = [&compSolution](const DArray2 &phi, const DArray2 &gg, int i, int k) {
        return compSolution(phi, i, k) + gg(i, k);
    };

    // k, i: jf(i, k) <- aa1(i+-1, k+-1)
    for (int k = 1; k < km + 1; k++) {
        double s = (1.5 * aa1(2, k) - 4.5 * aa1(1, k)) / hr2 +
                   (aa1(1, k + 1) - 2.0 * aa1(1, k) + aa1(1, k - 1)) / hz2;
        jf(1, k) = -s;
        for (int i = 2; i < 2 * im + 1; i++) {
            jf(i, k) = -compSolution(aa1, i, k);
        }
    }

    doOutput("aa1 aa1", aa1);
    doOutput("jf jf", jf);

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

    // i, k: gg(i, k) <- jf(i, k), jf(i, k + 1)
    // i, k: phi(i, k) <- aa1(i, k), aa1(i, k + 1)
    for (int i = 1; i < 2 * im + 1; i++) {
        for (int k = 1; k < km; k++) {
            gg(i, k) = jf(i, k + 1) - jf(i, k);
            phi1(i, k) = aa1(i, k + 1) - aa1(i, k);
        }
        gg(i, 0) = 0.0;
        gg(i, km) = 0.0;
        phi1(i, 0) = 0.0;
        phi1(i, km) = 0.0;
    }

    doOutput("gg gg", gg);

/*
    вычисление синусов

      do k=1,2*km
         ds(k)=dsin(c*k)
      enddo
*/

    const auto dsins = computeSins(2 * km, c);

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

    auto computeFFTAux = [&](const DArray2 &input, DArray2 &output, double coeff) {
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int j = 1; j < km; j++) {
                double s = 0.0;
                int k1 = 0;
                for (int k = 1; k < km; k++) {
                    k1 = k1 + j;
                    if (k1 >= 2 * km) {
                        k1 = k1 - 2 * km;
                    }
                    s += input(i, k) * dsins[k1];
                }
                output(i, j) = s * coeff;
            }
            output(i, 0) = 0.0;
            output(i, km) = 0.0;
        }
    };

    auto computeFFT = [&computeFFTAux, km](const DArray2 &input, DArray2 &output) {
        computeFFTAux(input, output, 2.0 / km);
    };

    auto computeFFTInverse = [&computeFFTAux](const DArray2 &input, DArray2 &output) {
        computeFFTAux(input, output, 1.0);
    };

    // i, j: bb(i, j) <- k, gg(i, k)
    // i, j: ff1(i, j) <- k, phi1(i, k)
    computeFFT(gg, bb);
    computeFFT(phi1, ff1);

    doOutput("bb bb", bb);

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

    // j, i: al(i) <- al(i - 1)
    // j, i: be(i) <- be(i - 1)
    // j, i: ff(i, j) <- al(i), be(i), ff(i + 1, j)
    for (int j = 1; j < km; j++) {
        const double ss = sin(c * (j + 1 - 1.0) / 2.0);
        double s = 9.0 / (2.0 * hr2) + (4.0 / hz2) * ss * ss;
        al[1] = 3.0 / (2.0 * hr2 * s);
        be[1] = bb(1, j) / s;
        for (int i = 2; i < 2 * im + 1; i++) {
            const auto dsin = sin(c * (j + 1 - 1.0) / 2.0);
            s = (2.0 * ((i + 1 - 1.5) / hr) * ((i + 1 - 1.5) / hr)) / ((i + 1 - 1.0) * (i + 1 - 2.0)) +
                (4.0 / hz2) * dsin * dsin -
                al[i - 1] * (i + 1 - 2.5) / ((i + 1 - 2.0) * hr2);
            al[i] = (i + 1 - 0.5) / (s * (i + 1 - 1.0) * hr2);
            be[i] = (be[i - 1] * (i + 1 - 2.5) / ((i + 1 - 2.0) * hr2) + bb(i, j)) / s;
        }
        ff(2 * im + 1, j) = 0.0;
        for (int i = 2 * im; i >= 1; i--) {
            ff(i, j) = al[i] * ff(i + 1, j) + be[i];
        }
    }

    doOutput("ff ff", ff);
    doOutput("ff1 ff1", ff1);

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
            enddo
            phi(i,k)=s1
         enddo
         phi(i,1)=0.d0
         phi(i,km+1)=0.d0
      enddo     !   i
*/

    // i, k: phi(i, k) <- j, ff(i, j)
    computeFFTInverse(ff, phi);

    doOutput("phi phi", phi);
    doOutput("phi1 phi1", phi1);

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

    auto compSolutionDifference = [&](const DArray2 &phi, const DArray2 &gg) -> DArray2 {
        DArray2 diff(im + 2, km + 2);
        for (int i = 2; i < im + 1; i++) {
            for (int k = 1; k < km + 1; k++) {
                diff(i, k) = std::abs(compCheck(phi, gg, i, k));
            }
        }
        return diff;
    };

    // i, k: dd(i, k) <- phi(i+-1, k+-1)
    doOutput("proverka2 dd dd", compSolutionDifference(phi, gg));

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

    // i: al(i) <- al(i - 1)
    // i: be(i) <- be(i - 1), phi(i + 1, 1), jf(i + 1, 1)
    // i: aa(i + 1, 1) <- al(i), be(i), aa(i + 2, 1)
    al[0] = 1.0 / 3.0;
    be[0] = (2.0 * hr2 / 9.0) * (jf(1, 1) + phi(1,1) / hz2);

    for (int i = 1; i < 2 * im; i++) {
        double s = 2.0 * (i + 1 - 0.5) * (i + 1 - 0.5) / ((i + 1) * (i + 1 - 1.0)) - al[i - 1] * (i + 1 - 1.5) / (i + 1 - 1.0);
        al[i] = (i + 1 + 0.5) / ((i + 1) * s);
        be[i] = (be[i - 1] * (i + 1 - 1.5) / (i + 1 - 1.0) + hr2 * (jf(i + 1, 1) + phi(i + 1, 1) / hz2)) / s;
    }

    aa(2 * im + 1, 1) = 0.0;
    for (int i = 2 * im - 1; i >= 0; i--) {
        aa(i + 1, 1) = al[i] * aa(i + 2, 1) + be[i];
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

    // i, k: aa(i, k + 1) <- aa(i, k), phi(i, k)
    for (int i = 1; i < im + 2; i++) {
        aa(i, 0) = aa(i, 1);
        for (int k = 1; k < km + 1; k++) {
            aa(i, k + 1) = aa(i, k) + phi(i, k);
        }
    }

    doOutput("aa aa", aa);

/*
    proverka3 решения dd dd

      do i=3,im+1
         do k=2,km+1
            dd(i,k)=(((i-0.5d0)*aa(i+1,k)-(i-1.5d0)*aa(i,k))/(i-1.d0)-
     =        ((i-1.5d0)*aa(i,k)-(i-2.5d0)*aa(i-1,k))/(i-2.d0))/hr**2+
     =        (aa(i,k+1)-2.d0*aa(i,k)+aa(i,k-1))/hz**2+jf(i,k)
         enddo
      enddo
*/

    // i, k: dd(i, k) <- aa(i+-1, k+-1), jf(i, k)
    doOutput("dd dd", compSolutionDifference(aa, jf));

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

    // k, i: bz(i, k) <- aa(i, k), aa(i + 1, k)
    for (int k = 0; k < km + 2; k++) {
        bz(0, k) = 4.0 * aa(1, k) / hr;
        for (int i = 1; i < im + 2; i++) {
            bz(i, k) = ((i + 1 - 0.5) * aa(i + 1, k) - (i + 1 - 1.5) * aa(i, k)) / (hr * (i + 1 - 1.0));
        }
    }

    // k, i: br(i, k) <- aa(i, k), aa(i, k + 1)
    for (int k = 0; k < km + 1; k++) {
        for (int i = 0; i < im + 2; i++) {
            br(i, k) = -(aa(i, k + 1) - aa(i, k)) / hz;
        }
    }

    if (file_output) {
        out_lst.close();
    }

    auto time = full_time.time();
    cout << "Time: " << time << endl;

    return 0;
}
