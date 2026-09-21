#include "defs.h"
#include "block_decomp.h"
#include "comm_util.h"
#include "../common/output.h"
#include "../common/timer.h"

#include <mpi.h>

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
    int nn = np;
    int j = 1;
    for (int i = 1; i <= nn; i++) {
        if (i < j) {
            std::swap(a[i - 1], a[j - 1]);
        }
        int m = nn / 2;
        do {
            if (j <= m) {
                break;
            }
            j -= m;
            m /= 2;
        } while (m >= 1);
        j += m;
    }
    int mm = 1;
    while (mm < nn) {
        const int ii = 2 * mm;
        const double th = pi2 / (n * isi >= 0 ? std::abs(ii) : -std::abs(ii));
        const auto sin_th = std::sin(th);
        const auto sin_th_d2 = std::sin(th * 0.5);
        Complex w1 {-2.0 * sin_th_d2 * sin_th_d2, sin_th};
        Complex w = 1.0;
        for (int m = 1; m <= mm; m++) {
            for (int i = m; i <= nn; i += ii) {
                const Complex t = w * a[i + mm - 1];
                a[i + mm - 1] = a[i - 1] - t;
                a[i - 1] += t;
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

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    const size_t imp = im + 2;
    const size_t imp2 = 2 * imp;
    const size_t kmp = km + 2;

    BlockDecomposition km_decomp(kmp, size);
    const size_t km_bsize = km_decomp.getBlockSize(rank);
    const auto my_km_range = km_decomp.getRange(rank);

    BlockDecomposition km4_decomp(4 * km, size);
    const size_t km4_bsize = km4_decomp.getBlockSize(rank);
    const auto my_km4_range = km4_decomp.getRange(rank);

    DArray2 jf(imp2, km_bsize), aa1(imp2, km_bsize);
    DArray2 bb(imp2, km4_bsize), ff(imp2, km4_bsize), phi(imp2, km_bsize);

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
    if (file_output && rank == 0) {
        out_lst.open("output2.lst");
    }

    auto outputArray = [rank, full_output, &out_lst](const std::string &header, const DArray2 &data, const std::array<int, 4> &range) {
        if (rank == 0) {
            if (full_output) {
                output(header, data, out_lst);
            } else {
                output(header, data, range, out_lst);
            }
        }
    };

    auto gatherAndOutput = [rank, size, file_output, &outputArray](const std::string &header, const DArray2 &local_data, const BlockDecomposition &decomp, const std::array<int, 4> &range) {
        if (file_output) {
            const auto arr = gatherArrayK(local_data, decomp, rank, size);
            outputArray(header, arr, range);
        }
    };

    double ft_time = 0;
    double prog_time = 0;
    double comp_time = 0;
    double shadow_time = 0;

    auto syncShadows = [rank, size, &shadow_time](DArray2 &array) {
        Timer tm;
        auto col_type = makeColType(array);
        syncShadowsK(array, col_type, rank, size);
        MPI_Type_free(&col_type);
        shadow_time += tm.time();
    };

    // Init functions

    auto initTestSolution = [im, km, zm, hz, hr2, &my_km_range, &comp_time](DArray2 &output, int m) {
        Timer tm;
        const double a0 = -0.1;
        const double a = 1e-3;
        const double d = 1.0 + 1e-5 * m;
        const auto k_start = my_km_range.localStart(0);
        const auto k_end = my_km_range.localEnd(km + 2);
        for (int k = k_start; k < k_end; k++) {
            const auto kk = my_km_range.toGlobal(k);
            const double z = hz * (kk + 1 - 1.5);
            const double z2 = z * z;
            const double s = a * z2 * (z - 1.5 - zm) + d;
            for (int i = 1; i < 2 * im + 2; i++) {
                output(i, k) = s * (i + 1 - 2 * im - 2.0) * (i + 1 - 1.5) * hr2;
            }
            output(0, k) = -output(1, k);
            output(2 * im + 1, k) = 0.0;
        }
        if (my_km_range.hasIndex(0)) {
            const auto dst_loc = my_km_range.toLocal(0);
            const auto src_loc = my_km_range.toLocal(1);
            for (int i = 0; i < 2 * im + 2; i++) {
                output(i, dst_loc) = output(i, src_loc);
            }
        }
        if (my_km_range.hasIndex(km + 1)) {
            const auto dst_loc = my_km_range.toLocal(km + 1);
            const auto src_loc = my_km_range.toLocal(km);
            for (int i = 0; i < 2 * im + 2; i++) {
                output(i, dst_loc) = output(i, src_loc);
            }
        }
        comp_time += tm.time();
    };

    auto getSolution = [rhr2, rhz2](const DArray2 &input, int i, int k) {
        return (((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (i) -
                ((i - 0.5) * input(i, k) - (i - 1.5) * input(i - 1, k)) / (i - 1.0)) * rhr2 +
               (input(i, k + 1) - 2.0 * input(i, k) + input(i, k - 1)) * rhz2;
    };

    auto initTestCurrent = [im, km, rhr2, rhz2, &getSolution, &syncShadows, &my_km_range, &comp_time](DArray2 &input, DArray2 &output) {
        syncShadows(input);
        Timer tm;
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km + 1);
        for (int k = k_start; k < k_end; k++) {
            double s = (1.5 * input(2, k) - 4.5 * input(1, k)) * rhr2 +
                       (input(1, k + 1) - 2.0 * input(1, k) + input(1, k - 1)) * rhz2;
            output(1, k) = -s;
            for (int i = 2; i < 2 * im + 1; i++) {
                output(i, k) = -getSolution(input, i, k);
            }
        }
        comp_time += tm.time();
    };

    /*auto computeDifference = [&](DArray2 &input, DArray2 &output) {
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = 1; k < km; k++) {
                output(i, k) = input(i, k + 1) - input(i, k);
            }
            output(i, 0) = 0.0;
            output(i, km) = 0.0;
        }
    };*/

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

    auto computeProgonka = [im, km, imp2, hr, hr2, hz, hz2, c, &my_km4_range, &prog_time, &comp_time](const DArray2 &input, DArray2 &output) {
        Timer tm;
        DArray1 al(imp2), be(imp2);
        const auto k_start = my_km4_range.localStart(0);
        const auto k_end = my_km4_range.localEnd(4 * km);
        for (int k = k_start; k < k_end; k++) {
            const int kk = my_km4_range.toGlobal(k);
            const int k1 = (kk >= 2 * km ? kk - 2 * km : kk);
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
        comp_time += tm.time();
    };

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

    auto computeSolutionDifference = [im, km, &my_km_range, &comp_time](DArray2 &first, const DArray2 &second) -> DArray2 {
        DArray2 output(im + 2, my_km_range.size());
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km);
        Timer tm;
        for (int i = 2; i < im + 1; i++) {
            for (int k = k_start; k < k_end; k++) {
                output(i, k) = first(i, k) - second(i, k);
            }
        }
        comp_time += tm.time();
        return output;
    };

    // The main program's body

    Timer work_timer;

    for (int m = 1; m <= 1000; m++) {
        initTestSolution(aa1, m);
        initTestCurrent(aa1, jf);

        computeFFT(jf, bb);
        computeProgonka(bb, ff);
        computeFFTInverse(ff, phi);
    }

    gatherAndOutput("phi phi", phi, km_decomp, {0, 7, 0, km + 2});
    gatherAndOutput("proverka 2 dd dd", computeSolutionDifference(phi, aa1), km_decomp, {0, 7, 0, 6});

    auto work_time = work_timer.time();

    double time = 0;
    MPI_Reduce(&work_time, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (file_output && rank == 0) {
        out_lst.close();
    }

    std::ostringstream out;

    if (rank == 0) {
        out << "Im: " << im << ", Km: " << km << ", Nodes: " << size << endl;
        out << "TIME: " << time << endl;
    }

    out << rank << ": Work time: " << work_time <<
        ", Comp time: " << comp_time <<
        ", Shadow time: " << shadow_time <<
        //", Reduce time: " << reduce_time <<
        ", FT: " << ft_time <<
        ", Prog: " << prog_time << endl;

    std::cout << out.str();

    MPI_Finalize();

    return 0;
}
