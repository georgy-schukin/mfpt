#include "defs.h"
#include "block_decomp.h"
#include "comm_util.h"
#include "../common/output.h"
#include "../common/timer.h"
#include "distributed_array2d.h"

#include <mpi.h>

#include <cmath>
#include <string>
#include <chrono>
#include <functional>
#include <sstream>

using namespace std;
using DDArray2 = DistributedArray2D;

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
            MPI_Finalize();
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int km = (argc > 2) ? stoi(argv[2]) : KM_DEF;
    const bool file_output = (argc > 3) ? stoi(argv[3]) : FOUT_DEF;
    const bool full_output = (argc > 4) ? stoi(argv[4]) : FULL_OUTPUT_DEF;

    const size_t ims = im + 2;
    const size_t ims2 = 2 * im + 2;

    // 1d decomposition by k.
    const auto kms_decomp = BlockDecomposition(km, size).grown(1); // km + 2
    const auto im_km_decomp = BlockDecomposition2D(ims, kms_decomp);
    const auto im2_km_decomp = BlockDecomposition2D(ims2, kms_decomp);

    DDArray2 br(im_km_decomp, rank), bf(im_km_decomp, rank), bz(im_km_decomp, rank);
    DDArray2 aa(im2_km_decomp, rank, 1), jf(im2_km_decomp, rank, 1), aa1(im2_km_decomp, rank, 1), phi(im2_km_decomp, rank, 1);
    DDArray2 gg(im2_km_decomp, rank), bb(im2_km_decomp, rank), ff(im2_km_decomp, rank);
    DDArray2 phi1(im2_km_decomp, rank), ff1(im2_km_decomp, rank);

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

    auto gatherAndOutput = [&](const std::string &header, const DDArray2 &local_data) {
        if (file_output) {
            const auto arr = local_data.gather(0);
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
    double sweep_time = 0;
    double sol_time = 0;
    double mag_time = 0;

    // Init functions

    auto syncShadows = [&shadow_time](DDArray2 &array) {
        Timer tm;
        array.syncShadows();
        shadow_time += tm.time();
    };

    auto initTestSolution = [im, km, hr, hz, rm, zm, &comp_time](DDArray2 &output) {
        Timer tm;
        const double a0 = -0.1;
        const double a = 1.0;
        const double d = 1.0;
        const auto k_range = output.range();
        const auto k_start = k_range.localStart(0);
        const auto k_end = k_range.localEnd(km + 2);
        for (int k = k_start; k < k_end; k++) {
            const auto kk = k_range.toGlobal(k);
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

    auto getSolution = [hr2, hz2](const DDArray2 &input, int i, int k) {
        return (((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (i) -
                ((i - 0.5) * input(i, k) - (i - 1.5) * input(i - 1, k)) / (i - 1.0)) / hr2 +
               (input(i, k + 1) - 2.0 * input(i, k) + input(i, k - 1)) / hz2;
    };

    auto initTestCurrent = [im, km, rhr2, rhz2, &syncShadows, &getSolution, &comp_time](DDArray2 &input, DDArray2 &output) {
        syncShadows(input);
        Timer tm;
        const auto k_range = input.range();
        const auto k_start = k_range.localStart(1);
        const auto k_end = k_range.localEnd(km + 1);
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

    auto setOnK = [](DDArray2 &arr, int i, int k, double value) {
        if (arr.range().hasIndex(k)) {
            arr(i, arr.range().toLocal(k)) = value;
        }
    };

    auto computeDifference = [im, km, &syncShadows, &setOnK, &comp_time](DDArray2 &input, DDArray2 &output) {
        syncShadows(input);
        Timer tm;
        const auto k_range = input.range();
        const auto k_start = k_range.localStart(1);
        const auto k_end = k_range.localEnd(km);
        for (int i = 1; i < 2 * im + 1; i++) {
            for (int k = k_start; k < k_end; k++) {
                output(i, k) = input(i, k + 1) - input(i, k);
            }
            setOnK(output, i, 0, 0.0);
            setOnK(output, i, km, 0.0);
        }
        comp_time += tm.time();
    };

    auto computeFFTAux = [im, km, size, &vector_sum_op, &ft_time, &reduce_time, &comp_time](const DDArray2 &input, DDArray2 &output, const DArray1 &dsins, double coeff) {
        Timer ftm;
        const auto k_decomp = input.decomp(1);
        const auto k_range = input.range();
        const auto k_start = k_range.localStart(1);
        const auto k_end = k_range.localEnd(km);
        const auto k_global_start = k_range.toGlobal(k_start);
        const auto dsins_size = dsins.size();
        for (int r = 0; r < size; r++) {
            Timer tm;
            const auto curr_j_range = k_decomp.range(r);
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

    auto computeFFT = [&computeFFTAux, km](const DDArray2 &input, DDArray2 &output, const DArray1 &dsins) {
        computeFFTAux(input, output, dsins, 2.0 / km);
    };

    auto computeFFTInverse = [&computeFFTAux](const DDArray2 &input, DDArray2 &output, const DArray1 &dsins) {
        computeFFTAux(input, output, dsins, 1.0);
    };

    auto computeSweep = [im, km, ims2, hr, hz, hr2, hz2, c, &comp_time, &sweep_time](const DDArray2 &input, DDArray2 &output) {
        Timer tm;
        DArray1 al(ims2), be(ims2);
        const auto k_range = input.range();
        const auto k_start = k_range.localStart(1);
        const auto k_end = k_range.localEnd(km);
        for (int k = k_start; k < k_end; k++) {
            const auto kk = k_range.toGlobal(k);
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
        sweep_time += tm.time();
    };

    auto computeSolutionDifference = [im, km, &syncShadows, &getSolution, &comp_time](DDArray2 &first, const DDArray2 &second) -> DDArray2 {
        DDArray2 output(first.decomp2d(), first.rank());
        const auto k_range = first.range();
        const auto k_start = k_range.localStart(1);
        const auto k_end = k_range.localEnd(km + 1);
        Timer tm;
        for (int i = 2; i < im + 1; i++) {
            for (int k = k_start; k < k_end; k++) {
                output(i, k) = std::abs(getSolution(first, i, k) + second(i, k));
            }
        }
        comp_time += tm.time();
        return output;
    };

    auto computeSolutionFirstColumn = [im, km, ims2, hr, hz, hr2, hz2, &comp_time](const DDArray2 &phi, const DDArray2 &jf, DDArray2 &aa) {
        Timer tm;
        DArray1 al(ims2), be(ims2);
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
        comp_time += tm.time();
    };

    auto computeSolutionRows = [im, km, ims2, hr, hz, hr2, hz2, &comp_time](const DDArray2 &phi, DDArray2 &aa) {
        const auto i_range = aa.range(0);
        const auto i_start = i_range.localStart(1);
        const auto i_end = i_range.localEnd(im + 2);
        Timer tm;
        for (int i = i_start; i < i_end; i++) {
            aa(i, 0) = aa(i, 1);
            for (int k = 1; k < km + 1; k++) {
                aa(i, k + 1) = aa(i, k) + phi(i, k);
            }
        }
        comp_time += tm.time();
    };

    auto computeSolution = [&sol_time, &computeSolutionFirstColumn, &computeSolutionRows](const DDArray2 &phi, const DDArray2 &jf, DDArray2 &aa) {
        Timer tm;
        computeSolutionFirstColumn(phi, jf, aa);
        const BlockDecomposition2D by_rows{BlockDecomposition(phi.decomp(0).fullSize(), phi.numOfNodes()), phi.decomp(1).fullSize()};
        DDArray2 phi_rows(by_rows, phi.rank()), aa_rows(by_rows, aa.rank());
        phi_rows.combineFrom(phi);
        aa_rows.combineFrom(aa);
        computeSolutionRows(phi_rows, aa_rows);
        aa.combineFrom(aa_rows);
        sol_time += tm.time();
    };

    auto computeMagnetics = [im, km, hr, hz, &syncShadows, &comp_time, &mag_time](DDArray2 &input, DDArray2 &bz, DDArray2 &br) {
        Timer mtm;
        Timer tm;
        const auto k_range = input.range();
        const int k_start = k_range.localStart(0);
        const int k_end = k_range.localEnd(km + 2);
        for (int k = k_start; k < k_end; k++) {
            bz(0, k) = 4.0 * input(1, k) / hr;
            for (int i = 1; i < im + 2; i++) {
                bz(i, k) = ((i + 0.5) * input(i + 1, k) - (i - 0.5) * input(i, k)) / (hr * (i));
            }
        }
        comp_time += tm.time();

        syncShadows(input);
        tm.reset();
        const int k_end2 = k_range.localEnd(km + 1);
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

    initTestSolution(aa1);
    initTestCurrent(aa1, jf);

    computeDifference(jf, gg);
    computeDifference(aa1, phi1);

    const auto dsins = computeSins(2 * km, c);

    computeFFT(gg, bb, dsins);
    computeFFT(phi1, ff1, dsins);

    computeSweep(bb, ff);

    computeFFTInverse(ff, phi, dsins);

    computeSolution(phi, jf, aa);

    computeMagnetics(aa, bz, br);

    gatherAndOutput("aa1 aa1", aa1);
    gatherAndOutput("jf jf", jf);
    gatherAndOutput("gg gg", gg);
    gatherAndOutput("bb bb", bb);
    gatherAndOutput("ff ff", ff);
    gatherAndOutput("ff1 ff1", ff1);
    gatherAndOutput("phi phi", phi);
    gatherAndOutput("phi1 phi1", phi1);
    gatherAndOutput("proverka2 dd dd", computeSolutionDifference(phi, gg));
    gatherAndOutput("aa aa", aa);
    gatherAndOutput("dd dd", computeSolutionDifference(aa, jf));

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
        ", Sweep: " << sweep_time <<
        ", Sol: " << sol_time <<
        ", Mag: " << mag_time << endl;

    std::cout << out.str();

    MPI_Finalize();

    return 0;
}
