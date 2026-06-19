#include "shadowed_array2d.h"
#include "block_decomp.h"

#include <mpi.h>

#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <chrono>
#include <functional>
#include <exception>

using namespace std;

using DArray2 = ShadowedArray2D<double>;
using DArray1 = std::vector<double>;

const int TAG_PREV = 0;
const int TAG_NEXT = 1;
const int GATHER_TAG = 2;

void print(const DArray2 &data, int i_start, int i_end, int j_start, int j_end, ofstream &out) {
    out << setw(7) << "";
    for (int i = i_start; i < i_end; i++) {
        out << setw(3) << i + 1;
        if (i < i_end - 1) {
            out << setw(6) << "";
        }
    }
    out << std::endl;
    for (int j = j_end - 1; j >= j_start; j--) {
        out << setw(3) << j + 1 << setw(1) << "";
        for (int i = i_start; i < i_end; i++) {
            out << setw(9) << std::fixed << setprecision(3) << data(i, j);
        }
        out << std::endl;
    }
}

void output(const string &header, const DArray2 &data, int i_start, int i_end, int j_start, int j_end, ofstream &out) {
    out << "\n";
    out << " " << header << "\n";
    print(data, i_start, i_end, j_start, j_end, out);
}

void output(const string &header, const DArray2 &data, const std::array<int, 4> &range, ofstream &out) {
    output(header, data, range[0], range[1], range[2], range[3], out);
}

MPI_Datatype makeColType(const DArray2 &array) {
    MPI_Datatype col_type;
    if (array.isRowMajorOrder()) {
        MPI_Datatype temp_type;
        MPI_Type_vector(array.size(0), 1, array.fullSize(1), MPI_DOUBLE, &temp_type);
        MPI_Type_create_resized(temp_type, 0, sizeof(double), &col_type);
    } else {
        MPI_Type_contiguous(array.size(0), MPI_DOUBLE, &col_type);
    }
    MPI_Type_commit(&col_type);
    return col_type;
}

MPI_Datatype makeDataSendType(const DArray2 &array) {
    MPI_Datatype send_type;
    MPI_Type_vector(array.size(0), array.size(1), array.size(1) + 2 * array.shadowSize(1), MPI_DOUBLE, &send_type);
    MPI_Type_commit(&send_type);
    return send_type;
}

void syncKPrev(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    if (rank > 0) {
        MPI_Request req[2];
        // Send data column.
        MPI_Isend(&arr(0, 0), 1, col_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &req[0]);
        // Receive in shadow.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), size_t(0)), 1, col_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    }
}

void syncKNext(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    if (rank < size - 1) {
        MPI_Request req[2];
        // Send data column.
        MPI_Isend(&arr(size_t(0), arr.size(1) - 1), 1, col_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &req[0]);
        // Receive in shadow.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), arr.fullSize(1) - 1), 1, col_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    }
}

void syncK(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    syncKPrev(arr, col_type, rank, size);
    syncKNext(arr, col_type, rank, size);
}

DArray2 gatherArrayK(const DArray2 &local_data, const BlockDecomposition &k_decomp, int im_size, int km_size, int rank, int size) {
    std::vector<MPI_Request> reqs;
    if (rank == 0) {
        reqs.resize(size + 1);
    } else {
        reqs.resize(1);
    }

    auto send_type = makeDataSendType(local_data);

    MPI_Isend(&local_data(0, 0), 1, send_type, 0, GATHER_TAG, MPI_COMM_WORLD, &reqs[0]);

    if (rank == 0) {
        std::vector<DArray2> parts;
        for (int i = 0; i < size; i++) {
            parts.push_back(DArray2(im_size, k_decomp.getBlockSize(i)));
        }
        for (int i = 0; i < size; i++) {
            MPI_Irecv(parts[i].data(), parts[i].size(), MPI_DOUBLE, i, GATHER_TAG, MPI_COMM_WORLD, &reqs[i + 1]);
        }
        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

        DArray2 arr(im_size, km_size);
        for (int r = 0; r < size; r++) {
            const auto &part = parts[r];
            const auto k_shift = k_decomp.getBlockShift(r);
            for (int i = 0; i < (int)part.size(0); i++) {
                for (int k = 0; k < (int)part.size(1); k++) {
                    arr(i, k + k_shift) = part(i, k);
                }
            }
        }
        MPI_Type_free(&send_type);
        return arr;
    } else {
        MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
        MPI_Type_free(&send_type);
    }
    return DArray2();
}

DArray1 computeSins(size_t size, double coeff) {
    DArray1 dsins(size);
    for (int k = 0; k < size; k++) {
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
                " [file_output=" << FOUT_DEF << "]" << endl;
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int km = (argc > 2) ? stoi(argv[2]) : KM_DEF;
    const bool file_output = (argc > 3) ? stoi(argv[3]) : FOUT_DEF;

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
    DArray1 sb(ims2), jb(ims2), al(ims2), be(ims2);
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

    auto gatherAndOutput = [&](const std::string &header, const DArray2 &local_data) {
        const auto arr = gatherArrayK(local_data, kms_decomp, local_data.size(0), kms, rank, size);
        if (rank == 0) {
            output(header, arr, output_range, out_lst);
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

    auto ts = chrono::steady_clock::now();

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

    auto initSolution = [&](DArray2 &output) {
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
    };

    // k, i: aa1(i, k) <- expr
    initSolution(aa1);

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

    auto initCurrent = [&](DArray2 &input, DArray2 &output) {
        syncK(input, col_type, rank, size);
        for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km + 1); k++) {
            const double s = (1.5 * input(2, k) - 4.5 * input(1, k)) / hr2 +
                             (input(1, k + 1) - 2.0 * input(1, k) + input(1, k - 1)) / hz2;
            output(1, k) = -s;
            for (int i = 2; i < 2 * im + 1; i++) {
                output(i, k) = -compSolution(input, i, k);
            }
        }
    };

    // k, i: jf(i, k) <- aa1(i+-1, k+-1)
    initCurrent(aa1, jf);

    if (file_output) {
        gatherAndOutput("aa1 aa1", aa1);
        gatherAndOutput("jf jf", jf);
    }

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
    auto doOnK = [&my_km_range](int k, std::function<void(int)> f) {
        if (my_km_range.hasIndex(k)) {
            const auto kk = my_km_range.toLocal(k);
            f(kk);
        }
    };

    auto computeDifference = [&](DArray2 &input, DArray2 &output) {
        syncKNext(input, col_type, rank, size);
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km); k++) {
                output(i, k) = input(i, k + 1) - input(i, k);
            }
            setOnK(output, i, 0, 0.0);
            setOnK(output, i, km, 0.0);
        }
    };

    // i, k: gg(i, k) <- jf(i, k), jf(i, k + 1)
    // i, k: phi1(i, k) <- aa1(i, k), aa1(i, k + 1)
    computeDifference(jf, gg);
    computeDifference(aa1, phi1);

    if (file_output) {
        gatherAndOutput("gg gg", gg);
    }

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
        for (int r = 0; r < size; r++) {
            DArray2 tmp(output.size(0), kms_decomp.getBlockSize(r));
            for (int i = 1; i < 2 * im + 1; i++) {
                const auto j_start = kms_decomp.getRange(r).localStart(1);
                const auto j_end = kms_decomp.getRange(r).localEnd(km);
                for (int j = j_start; j < j_end; j++) {
                    double s = 0.0;
                    for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km); k++) {
                        const auto kk = my_km_range.toGlobal(k);
                        const auto jj = my_km_range.toGlobal(j);
                        s += input(i, k) * dsins[(kk * jj) % dsins.size()];
                    }
                    tmp(i, j) = s * coeff;
                }
            }
            // Create tmp array with no shadows to use MPI Reduce.
            DArray2 output_tmp;
            if (rank == r) {
                output_tmp = DArray2(output.size(0), output.size(1));
            }
            MPI_Reduce(tmp.data(), output_tmp.data(), tmp.size(), MPI_DOUBLE, MPI_SUM, r, MPI_COMM_WORLD);
            copy(output_tmp, output);
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

    if (file_output) {
        gatherAndOutput("bb bb", bb);
    }

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

    // k, i: al(i) <- al(i - 1)
    // k, i: be(i) <- be(i - 1)
    // k, i: ff(i, k) <- al(i), be(i), ff(i + 1, k)
    for (int k = my_km_range.localEnd(1); k < my_km_range.localEnd(km); k++) {
        const auto kk = my_km_range.toGlobal(k);
        const double dsin = sin(c * kk / 2.0);
        double s = 9.0 / (2.0 * hr2) + (4.0 / hz2) * dsin * dsin;
        al[1] = 3.0 / (2.0 * hr2 * s);
        be[1] = bb(1, k) / s;
        for (int i = 2; i < 2 * im + 1; i++) {
            s = (2.0 * ((i - 0.5) / hr) * ((i - 0.5) / hr)) / ((i) * (i - 1.0)) +
                (4.0 / hz2) * dsin * dsin -
                al[i - 1] * (i - 1.5) / ((i - 1.0) * hr2);
            al[i] = (i + 0.5) / (s * (i) * hr2);
            be[i] = (be[i - 1] * (i - 1.5) / ((i - 1.0) * hr2) + bb(i, k)) / s;
        }
        ff(2 * im + 1, k) = 0.0;
        for (int i = 2 * im; i >= 1; i--) {
            ff(i, k) = al[i] * ff(i + 1, k) + be[i];
        }
    }

    if (file_output) {
        gatherAndOutput("ff ff", ff);
        gatherAndOutput("ff1 ff1", ff1);
    }

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

    if (file_output) {
        gatherAndOutput("phi phi", phi);
        gatherAndOutput("phi1 phi1", phi1);
    }

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

    auto computeSolutionDifference = [&](DArray2 &first, const DArray2 &second, DArray2 &output) {
        syncK(first, col_type, rank, size);
        for (int k = my_km_range.localStart(1); k < my_km_range.localEnd(km + 1); k++) {
            for (int i = 2; i < im + 1; i++) {
                output(i, k) = compSolution(first, i, k) + second(i, k);
            }
        }
    };

    // i, k: dd(i, k) <- phi(i+-1, k+-1), gg(i, k)
    computeSolutionDifference(phi, gg, dd);

    if (file_output) {
        gatherAndOutput("proverka2 dd dd", dd);
    }

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
    // TODO: parallelize
    al[0] = 1.0 / 3.0;
    be[0] = (2.0 * hr2 / 9.0) * (jf(1, 1) + phi(1,1) / hz2);

    for (int i = 1; i < 2 * im; i++) {
        double s = 2.0 * (i + 0.5) * (i + 0.5) / ((i + 1) * (i)) - al[i - 1] * (i - 0.5) / (i);
        al[i] = (i + 1.5) / ((i + 1) * s);
        be[i] = (be[i - 1] * (i - 0.5) / (i) + hr2 * (jf(i + 1, 1) + phi(i + 1, 1) / hz2)) / s;
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
    // TODO: parallelize
    for (int i = 1; i < im + 2; i++) {
        aa(i, 0) = aa(i, 1);
        for (int k = 1; k < km + 1; k++) {
            aa(i, k + 1) = aa(i, k) + phi(i, k);
        }
    }

    if (file_output) {
        gatherAndOutput("aa aa", aa);
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

    // i, k: dd(i, k) <- aa(i+-1, k+-1), jf(i, k)
    computeSolutionDifference(aa, jf, dd);

    if (file_output) {
        gatherAndOutput("dd dd", dd);
    }

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
    for (int k = my_km_range.localStart(0); k < my_km_range.localEnd(km + 2); k++) {
        bz(0, k) = 4.0 * aa(1, k) / hr;
        for (int i = 1; i < im + 2; i++) {
            bz(i, k) = ((i + 0.5) * aa(i + 1, k) - (i - 0.5) * aa(i, k)) / (hr * (i));
        }
    }

    // k, i: br(i, k) <- aa(i, k), aa(i, k + 1)
    syncKNext(aa, col_type, rank, size);
    for (int k = my_km_range.localStart(0); k < my_km_range.localEnd(km + 1); k++) {
        for (int i = 0; i < im + 2; i++) {
            br(i, k) = -(aa(i, k + 1) - aa(i, k)) / hz;
        }
    }

    if (file_output && rank == 0) {
        out_lst.close();
    }

    auto te = chrono::steady_clock::now();
    auto time = chrono::duration<double>(te - ts).count();
    cout << "Time: " << time << endl;

    MPI_Finalize();

    return 0;
}
