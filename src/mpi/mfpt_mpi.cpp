#include "common.h"
#include "block_decomp.h"
#include "comm_util.h"
#include "output.h"
#include "timer.h"

#include <mpi.h>

#include <cmath>
#include <string>
#include <chrono>
#include <functional>
#include <sstream>

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

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    // 1d decomposition by k.
    BlockDecomposition kms_decomp(kms, size);
    const size_t km_bsize = kms_decomp.getBlockSize(rank);
    const auto my_km_range = kms_decomp.getRange(rank);

    DArray2 br(ims, km_bsize), bf(ims, km_bsize), bz(ims, km_bsize);
    //DArray1 sb(ims2), jb(ims2), al(ims2), be(ims2);
    DArray2 aa(ims2, km_bsize, 0, 1), jf(ims2, km_bsize, 0, 1), aa1(ims2, km_bsize, 0, 1), phi(ims2, km_bsize, 0, 1);
    DArray2 gg(ims2, km_bsize), bb(ims2, km_bsize), ff(ims2, km_bsize);
    DArray2 dd(ims, km_bsize), phi1(ims2, km_bsize), ff1(ims2, km_bsize);
    //DArray1 dsins(2 * kms);

    auto col_type = makeColType(aa);

    const std::array<int, 4> output_range = {0, 7, 0, 6};

    ofstream out_lst;
    if (file_output && rank == 0) {
        out_lst.open("output.lst");
    }

    auto outputArray = [&](const std::string &header, const DArray2 &data) {
        if (rank == 0) {
            if (full_output) {
                output(header, data, {0, (int)data.size(0), 0, (int)data.size(1)}, out_lst);
            } else {
                output(header, data, output_range, out_lst);
            }
        }
    };

    auto gatherAndOutput = [&](const std::string &header, const DArray2 &local_data) {
        if (file_output) {
            const auto arr = gatherArrayK(local_data, kms_decomp, rank, size);
            outputArray(header, arr);
        }
    };

    const double pi = 3.14159265358979;
    const double c = pi / km;
    const double hr = 0.2;
    const double hz = 0.2;
    const double rm = im * hr;
    const double zm = km * hz;
    const double hr2 = hr * hr;
    const double hz2 = hz * hz;

    Timer work_timer;
    double comp_time = 0;
    double shadow_time = 0;
    double reduce_time = 0;

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

    auto initTestSolution = [&](DArray2 &output) {
        Timer tm;
        const double a0 = -0.1;
        const double a = 1.0;
        const double d = 1.0;
        for (int k = my_km_range.localStart(0); k < my_km_range.localEnd(km + 2); k++) {
            const auto kk = my_km_range.toGlobal(k);
            const double z = hz * (kk - 0.5);
            const double z2 = z * z;
            const double s = a0 * z2 * (z2 - 2.0 * zm * zm) + a * z2 * (z - 1.5 * zm) + d;
            for (int i = 1; i < 2 * im + 2; i++) {
                output(i, k) = s * (hr * (i - 0.5) * (2.0 * rm - hr * (i - 1.0)));
            }
            output(0, k) = -output(1, k);
        }
        comp_time += tm.time();
    };

    // k: 0..km+2, i: 1..2*im+2: aa1(i, k) <- expr
    initTestSolution(aa1);

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

    auto getSolution = [hr2, hz2](const DArray2 &phi, int i, int k) {
        return (((i + 0.5) * phi(i + 1, k) - (i - 0.5) * phi(i, k)) / (i) -
                ((i - 0.5) * phi(i, k) - (i - 1.5) * phi(i - 1, k)) / (i - 1.0)) / hr2 +
               (phi(i, k + 1) - 2.0 * phi(i, k) + phi(i, k - 1)) / hz2;
    };

    auto initTestCurrent = [&](DArray2 &input, DArray2 &output) {
        Timer tm;
        syncShadowsK(input, col_type, rank, size);
        shadow_time += tm.time();
        tm.reset();
        for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km + 1); k++) {
            const double s = (1.5 * input(2, k) - 4.5 * input(1, k)) / hr2 +
                             (input(1, k + 1) - 2.0 * input(1, k) + input(1, k - 1)) / hz2;
            output(1, k) = -s;
            for (int i = 2; i < 2 * im + 1; i++) {
                output(i, k) = -getSolution(input, i, k);
            }
        }
        comp_time += tm.time();
    };

    // k: 1..km+1, i: 2..2*im+1: jf(i, k) <- aa1(i+-1, k+-1)
    initTestCurrent(aa1, jf);

    gatherAndOutput("aa1 aa1", aa1);
    gatherAndOutput("jf jf", jf);

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

    auto setOnK = [&my_km_range](DArray2 &arr, int i, int k, double value) {
        if (my_km_range.hasIndex(k)) {
            arr(i, my_km_range.toLocal(k)) = value;
        }
    };

    auto computeDifference = [&](DArray2 &input, DArray2 &output) {
        Timer tm;
        syncShadowsK(input, col_type, rank, size);
        shadow_time += tm.time();
        tm.reset();
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km); k++) {
                output(i, k) = input(i, k + 1) - input(i, k);
            }
            setOnK(output, i, 0, 0.0);
            setOnK(output, i, km, 0.0);
        }
        comp_time += tm.time();
    };

    // i: 1..2*im+1, k: 1..km: gg(i, k) <- jf(i, k), jf(i, k + 1)
    // i: 1..2*im+1, k: 1..km: phi1(i, k) <- aa1(i, k), aa1(i, k + 1)
    computeDifference(jf, gg);
    computeDifference(aa1, phi1);

    gatherAndOutput("gg gg", gg);

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

    auto vector_sum_op = makeVectorSumOp();

    auto computeFFTAux = [&](const DArray2 &input, DArray2 &output, double coeff) {
        for (int r = 0; r < size; r++) {
            Timer tm;
            const auto curr_range = kms_decomp.getRange(r);
            const auto j_start = curr_range.localStart(1);
            const auto j_end = curr_range.localEnd(km);
            DArray2 tmp(output.size(0), curr_range.size(), 0.0, output.shadowSize(0), output.shadowSize(1));
            for (int i = 1; i < 2 * im + 1; i++) {
                for (int j = j_start; j < j_end; j++) {
                    double s = 0.0;
                    const auto jj = curr_range.toGlobal(j);
                    for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km); k++) {
                        const auto kk = my_km_range.toGlobal(k);
                        s += input(i, k) * dsins[(kk * jj) % dsins.size()];
                    }
                    tmp(i, j) = s * coeff;
                }
            }
            comp_time += tm.time();
            auto type = makeDataVectorType(tmp);
            tm.reset();
            MPI_Reduce(&tmp(0, 0), &output(0, 0), 1, type, vector_sum_op, r, MPI_COMM_WORLD);
            reduce_time += tm.time();
            MPI_Type_free(&type);
        }
    };

    auto computeFFT = [&computeFFTAux, km](const DArray2 &input, DArray2 &output) {
        computeFFTAux(input, output, 2.0 / km);
    };

    auto computeFFTInverse = [&computeFFTAux](const DArray2 &input, DArray2 &output) {
        computeFFTAux(input, output, 1.0);
    };

    // i: 1..2*im+1, j: 1..km: bb(i, j) <- k: 1..km, gg(i, k)
    // i: 1..2*im+1, j: 1..km: ff1(i, j) <- k: 1..km, phi1(i, k)
    computeFFT(gg, bb);
    computeFFT(phi1, ff1);

    gatherAndOutput("bb bb", bb);

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

    auto computeProgonka = [&](const DArray2 &input, DArray2 &output) {
        Timer tm;
        DArray1 al(ims2), be(ims2);
        for (int k = my_km_range.localEnd(1); k < my_km_range.localEnd(km); k++) {
            const auto kk = my_km_range.toGlobal(k);
            const double dsin = sin(c * kk / 2.0);
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
        comp_time += tm.time();
    };

    // k: 1..km, i: 2..2*im+1: al(i) <- al(i - 1)
    // k: 1..km, i: 2..2*im+1: be(i) <- be(i - 1), bb(i, k)
    // k: 1..km, i: 2*im..1: ff(i, k) <- al(i), be(i), ff(i + 1, k)
    // i: seq, k: par
    computeProgonka(bb, ff);

    gatherAndOutput("ff ff", ff);
    gatherAndOutput("ff1 ff1", ff1);

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

    // i: 1..2*im+1, k: 1..km: phi(i, k) <- j: 1..km, ff(i, j)
    computeFFTInverse(ff, phi);

    gatherAndOutput("phi phi", phi);
    gatherAndOutput("phi1 phi1", phi1);

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

    auto computeSolutionDifference = [&](DArray2 &first, const DArray2 &second) -> DArray2 {
        Timer tm;
        syncShadowsK(first, col_type, rank, size);
        shadow_time += tm.time();
        tm.reset();
        DArray2 output(im + 2, my_km_range.size());
        for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km + 1); k++) {
            for (int i = 2; i < im + 1; i++) {
                output(i, k) = std::abs(getSolution(first, i, k) + second(i, k));
            }
        }
        comp_time += tm.time();
        return output;
    };

    // i: 2..im+1, k: 1..km: dd(i, k) <- phi(i+-1, k+-1), gg(i, k)
    gatherAndOutput("proverka2 dd dd", computeSolutionDifference(phi, gg));

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

    решение во всей области

      do i=2,im+2
         aa(i,1)=aa(i,2)
         do k=2,km+1
            aa(i,k+1)=aa(i,k)+phi(i,k)
         enddo
      enddo
*/

    // TODO: parallelize
    auto compSolution = [&](const DArray2 &phi, const DArray2 &jf, DArray2 &aa) {
        Timer tm;
        DArray1 al(ims2), be(ims2);
        al[0] = 1.0 / 3.0;
        be[0] = (2.0 * hr2 / 9.0) * (jf(1, 1) + phi(1,1) / hz2);

        // i: 1..2*im: al(i) <- al(i - 1)
        // i: 1..2*im: be(i) <- be(i - 1), phi(i + 1, 1), jf(i + 1, 1)
        // i: seq
        for (int i = 1; i < 2 * im; i++) {
            double s = 2.0 * (i + 0.5) * (i + 0.5) / ((i + 1) * (i)) - al[i - 1] * (i - 0.5) / (i);
            al[i] = (i + 1.5) / ((i + 1) * s);
            be[i] = (be[i - 1] * (i - 0.5) / (i) + hr2 * (jf(i + 1, 1) + phi(i + 1, 1) / hz2)) / s;
        }

        // i: 2*im-1..0: aa(i + 1, 1) <- al(i), be(i), aa(i + 2, 1)
        // i: seq
        aa(2 * im + 1, 1) = 0.0;
        for (int i = 2 * im - 1; i >= 0; i--) {
            aa(i + 1, 1) = al[i] * aa(i + 2, 1) + be[i];
        }

        // i: 1..im+2, k: 1..km+1: aa(i, k + 1) <- aa(i, k), phi(i, k)
        // i: par, k: seq
        for (int i = 1; i < im + 2; i++) {
            aa(i, 0) = aa(i, 1);
            for (int k = 1; k < km + 1; k++) {
                aa(i, k + 1) = aa(i, k) + phi(i, k);
            }
        }
        comp_time += tm.time();
    };

    // Temporary fix: compute solution on a root(0) node
    auto phi_global = gatherArrayK(phi, kms_decomp, rank, size);
    auto jf_global = gatherArrayK(jf, kms_decomp, rank, size);
    DArray2 aa_global;
    if (rank == 0) {
        aa_global = DArray2(ims2, kms);
        compSolution(phi_global, jf_global, aa_global);
    }
    aa = scatterArrayK(aa_global, ims2, kms_decomp, 0, 1, rank, size);

    outputArray("aa aa", aa_global);

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

    // i: 2..im+1, k: 1..km+1: dd(i, k) <- aa(i+-1, k+-1), jf(i, k)
    gatherAndOutput("dd dd", computeSolutionDifference(aa, jf));

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
    auto computeMagnetics = [&](DArray2 &input, DArray2 &bz, DArray2 &br) {
        Timer tm;
        for (int k = my_km_range.localStart(0); k < my_km_range.localEnd(km + 2); k++) {
            bz(0, k) = 4.0 * input(1, k) / hr;
            for (int i = 1; i < im + 2; i++) {
                bz(i, k) = ((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (hr * (i));
            }
        }
        comp_time += tm.time();

        tm.reset();
        syncShadowsKNext(input, col_type, rank, size);
        shadow_time += tm.time();
        tm.reset();
        for (int k = my_km_range.localStart(0); k < my_km_range.localEnd(km + 1); k++) {
            for (int i = 0; i < im + 2; i++) {
                br(i, k) = -(input(i, k + 1) - input(i, k)) / hz;
            }
        }
        comp_time += tm.time();
    };

    // k: 0..km+2, i: 1..im+2: bz(i, k) <- aa(i, k), aa(i + 1, k)
    // k: 0..km+1, i: 0..im+2: br(i, k) <- aa(i, k), aa(i, k + 1)
    // i: par, k: par
    computeMagnetics(aa, bz, br);

    if (file_output && rank == 0) {
        out_lst.close();
    }

    auto work_time = work_timer.time();

    double time = 0;
    MPI_Reduce(&work_time, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    std::ostringstream out;

    if (rank == 0) {
        out << "Im: " << im << ", Km: " << km << ", Nodes: " << size << endl;
        out << "TIME: " << time << endl;
    }

    out << rank << ": Work time: " << work_time <<
        ", Comp time: " << comp_time <<
        ", Shadow time: " << shadow_time <<
        ", Reduce time: " << reduce_time << endl;

    std::cout << out.str();

    MPI_Finalize();

    return 0;
}
