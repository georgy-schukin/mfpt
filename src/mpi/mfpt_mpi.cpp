#include "defs.h"
#include "block_decomp.h"
#include "comm_util.h"
#include "../common/output.h"
#include "../common/timer.h"

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

    auto col_type = makeColType(aa);
    auto vector_sum_op = makeVectorSumOp();

    const std::array<int, 4> output_range = {0, 7, 0, 6};

    ofstream out_lst;
    if (file_output && rank == 0) {
        out_lst.open("output.lst");
    }

    auto outputArray = [&](const std::string &header, const DArray2 &data) {
        if (rank == 0) {
            if (full_output) {
                output(header, data, out_lst);
            } else {
                output(header, data, output_range, out_lst);
            }
        }
    };

    auto gatherAndOutput = [&](const std::string &header, const DArray2 &local_data) {
        if (file_output) {
            const auto arr = gatherArrayCols(local_data, kms_decomp, rank, size);
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
    const double rhr2 = 1.0 / hr2;
    const double rhz2 = 1.0 / hz2;

    double comp_time = 0;
    double shadow_time = 0;
    double reduce_time = 0;
    double ft_time = 0;
    double prog_time = 0;
    double sol_time = 0;
    double mag_time = 0;

    // Init functions

    auto syncShadows = [col_type, rank, size, &shadow_time](DArray2 &array) {
        Timer tm;
        syncShadowsCols(array, col_type, rank, size);
        shadow_time += tm.time();
    };

    auto syncShadowsNext = [col_type, rank, size, &shadow_time](DArray2 &array) {
        Timer tm;
        syncShadowsColsNext(array, col_type, rank, size);
        shadow_time += tm.time();
    };

    auto initTestSolution = [im, km, hr, hz, rm, zm, &my_km_range, &comp_time](DArray2 &output) {
        Timer tm;
        const double a0 = -0.1;
        const double a = 1.0;
        const double d = 1.0;
        const auto k_start = my_km_range.localStart(0);
        const auto k_end = my_km_range.localEnd(km + 2);
        for (int k = k_start; k < k_end; k++) {
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

    auto getSolution = [hr2, hz2](const DArray2 &input, int i, int k) {
        return (((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (i) -
                ((i - 0.5) * input(i, k) - (i - 1.5) * input(i - 1, k)) / (i - 1.0)) / hr2 +
               (input(i, k + 1) - 2.0 * input(i, k) + input(i, k - 1)) / hz2;
    };

    auto initTestCurrent = [im, km, rhr2, rhz2, &my_km_range, &syncShadows, &getSolution, &comp_time](DArray2 &input, DArray2 &output) {
        syncShadows(input);
        Timer tm;
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km + 1);
        for (int k = k_start; k < k_end; k++) {
            const double s = (1.5 * input(2, k) - 4.5 * input(1, k)) * rhr2 +
                             (input(1, k + 1) - 2.0 * input(1, k) + input(1, k - 1)) * rhz2;
            output(1, k) = -s;
            for (int i = 2; i < 2 * im + 1; i++) {
                output(i, k) = -getSolution(input, i, k);
            }
        }
        comp_time += tm.time();
    };

    auto setOnK = [&my_km_range](DArray2 &arr, int i, int k, double value) {
        if (my_km_range.hasIndex(k)) {
            arr(i, my_km_range.toLocal(k)) = value;
        }
    };

    auto computeDifference = [im, km, &my_km_range, &syncShadows, &setOnK, &comp_time](DArray2 &input, DArray2 &output) {
        syncShadows(input);
        Timer tm;
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km);
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = k_start; k < k_end; k++) {
                output(i, k) = input(i, k + 1) - input(i, k);
            }
            setOnK(output, i, 0, 0.0);
            setOnK(output, i, km, 0.0);
        }
        comp_time += tm.time();
    };

    auto computeFFTAux = [im, km, size, &vector_sum_op, &kms_decomp, &my_km_range, &ft_time, &reduce_time, &comp_time](const DArray2 &input, DArray2 &output, const DArray1 &dsins, double coeff) {
        Timer ftm;
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km);
        const auto k_global_start = my_km_range.toGlobal(k_start);
        const auto dsins_size = dsins.size();
        for (int r = 0; r < size; r++) {
            Timer tm;
            const auto curr_j_range = kms_decomp.getRange(r);
            const auto j_start = curr_j_range.localStart(1);
            const auto j_end = curr_j_range.localEnd(km);
            DArray2 tmp(output.size(0), curr_j_range.size(), 0.0, output.shadowSize(0), output.shadowSize(1));
            for (int i = 1; i < 2 * im + 1; i++) {
                size_t jj = curr_j_range.toGlobal(j_start);
                for (int j = j_start; j < j_end; j++, jj++) {
                    double s = 0.0;
                    size_t k1 = (k_global_start * jj) % dsins_size;
                    for (int k = k_start; k < k_end; k++) {
                        if (k1 >= dsins_size) {
                            k1 -= dsins_size;
                        }
                        s += input(i, k) * dsins[k1];
                        k1 += jj;
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
        ft_time += ftm.time();
    };

    auto computeFFT = [&computeFFTAux, km](const DArray2 &input, DArray2 &output, const DArray1 &dsins) {
        computeFFTAux(input, output, dsins, 2.0 / km);
    };

    auto computeFFTInverse = [&computeFFTAux](const DArray2 &input, DArray2 &output, const DArray1 &dsins) {
        computeFFTAux(input, output, dsins, 1.0);
    };

    auto computeProgonka = [im, km, ims2, hr, hz, hr2, hz2, c, &my_km_range, &comp_time, &prog_time](const DArray2 &input, DArray2 &output) {
        Timer tm;
        DArray1 al(ims2), be(ims2);
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km);
        for (int k = k_start; k < k_end; k++) {
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
        prog_time += tm.time();
    };

    auto computeSolutionDifference = [im, km, &my_km_range, &syncShadows, &getSolution, &comp_time](DArray2 &first, const DArray2 &second) -> DArray2 {
        Timer tm;
        DArray2 output(im + 2, my_km_range.size());
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km + 1);
        for (int k = k_start; k < k_end; k++) {
            for (int i = 2; i < im + 1; i++) {
                output(i, k) = std::abs(getSolution(first, i, k) + second(i, k));
            }
        }
        comp_time += tm.time();
        return output;
    };

    // TODO: parallelize
    auto compSolution = [im, km, ims2, hr, hz, hr2, hz2, &comp_time, &sol_time](const DArray2 &phi, const DArray2 &jf, DArray2 &aa) {
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
        sol_time += tm.time();
    };

    auto computeMagnetics = [im, km, hr, hz, &my_km_range, &syncShadowsNext, &comp_time, &mag_time](DArray2 &input, DArray2 &bz, DArray2 &br) {
        Timer mtm;
        Timer tm;
        const int k_start = my_km_range.localStart(0);
        const int k_end = my_km_range.localEnd(km + 2);
        for (int k = k_start; k < k_end; k++) {
            bz(0, k) = 4.0 * input(1, k) / hr;
            for (int i = 1; i < im + 2; i++) {
                bz(i, k) = ((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (hr * (i));
            }
        }
        comp_time += tm.time();

        syncShadowsNext(input);
        tm.reset();
        const int k_end2 = my_km_range.localEnd(km + 1);
        for (int k = k_start; k < k_end2; k++) {
            for (int i = 0; i < im + 2; i++) {
                br(i, k) = -(input(i, k + 1) - input(i, k)) / hz;
            }
        }
        comp_time += tm.time();
        mag_time += mtm.time();
    };

    // The main program

    Timer work_timer;

    // k: 0..km+2, i: 1..2*im+2: aa1(i, k) <- expr
    initTestSolution(aa1);

    // k: 1..km+1, i: 2..2*im+1: jf(i, k) <- aa1(i+-1, k+-1)
    initTestCurrent(aa1, jf);

    gatherAndOutput("aa1 aa1", aa1);
    gatherAndOutput("jf jf", jf);

    // i: 1..2*im+1, k: 1..km: gg(i, k) <- jf(i, k), jf(i, k + 1)
    // i: 1..2*im+1, k: 1..km: phi1(i, k) <- aa1(i, k), aa1(i, k + 1)
    computeDifference(jf, gg);
    computeDifference(aa1, phi1);

    gatherAndOutput("gg gg", gg);

    const auto dsins = computeSins(2 * km, c);

    // i: 1..2*im+1, j: 1..km: bb(i, j) <- k: 1..km, gg(i, k)
    // i: 1..2*im+1, j: 1..km: ff1(i, j) <- k: 1..km, phi1(i, k)
    computeFFT(gg, bb, dsins);
    computeFFT(phi1, ff1, dsins);

    gatherAndOutput("bb bb", bb);

    // k: 1..km, i: 2..2*im+1: al(i) <- al(i - 1)
    // k: 1..km, i: 2..2*im+1: be(i) <- be(i - 1), bb(i, k)
    // k: 1..km, i: 2*im..1: ff(i, k) <- al(i), be(i), ff(i + 1, k)
    // i: seq, k: par
    computeProgonka(bb, ff);

    gatherAndOutput("ff ff", ff);
    gatherAndOutput("ff1 ff1", ff1);

    // i: 1..2*im+1, k: 1..km: phi(i, k) <- j: 1..km, ff(i, j)
    computeFFTInverse(ff, phi, dsins);

    gatherAndOutput("phi phi", phi);
    gatherAndOutput("phi1 phi1", phi1);

    // i: 2..im+1, k: 1..km: dd(i, k) <- phi(i+-1, k+-1), gg(i, k)
    gatherAndOutput("proverka2 dd dd", computeSolutionDifference(phi, gg));

    // Temporary fix: compute solution on a root(0) node
    auto phi_global = gatherArrayCols(phi, kms_decomp, rank, size);
    auto jf_global = gatherArrayCols(jf, kms_decomp, rank, size);
    DArray2 aa_global;
    if (rank == 0) {
        aa_global = DArray2(ims2, kms);
        compSolution(phi_global, jf_global, aa_global);
    }
    aa = scatterArrayCols(aa_global, ims2, kms_decomp, 0, 1, rank, size);

    outputArray("aa aa", aa_global);

    // i: 2..im+1, k: 1..km+1: dd(i, k) <- aa(i+-1, k+-1), jf(i, k)
    gatherAndOutput("dd dd", computeSolutionDifference(aa, jf));

    // k: 0..km+2, i: 1..im+2: bz(i, k) <- aa(i, k), aa(i + 1, k)
    // k: 0..km+1, i: 0..im+2: br(i, k) <- aa(i, k), aa(i, k + 1)
    // i: par, k: par
    computeMagnetics(aa, bz, br);

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
        ", Reduce time: " << reduce_time <<
        ", FT: " << ft_time <<
        ", Prog: " << prog_time <<
        ", Sol: " << sol_time <<
        ", Mag: " << mag_time << endl;

    std::cout << out.str();

    MPI_Finalize();

    return 0;
}
