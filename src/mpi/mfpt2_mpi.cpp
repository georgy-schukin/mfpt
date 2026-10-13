#include "defs.h"
#include "block_decomp.h"
#include "comm_util.h"
#include "distributed_array2d.h"
#include "../common/output.h"
#include "../common/timer.h"
#include "async_op.h"

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

using DDArray2 = DistributedArray2D;
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
    const int REPEATS_DEF = 1000;
    const int FOUT_DEF = 1;
    const int FULL_OUTPUT_DEF = 0;
    const int DECOMP_BY_ROWS_DEF = 0;
    const int TASK_SIZE_DEF = 16;

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
                " [repeats=" << REPEATS_DEF << "]" <<
                " [file_output=" << FOUT_DEF << "]" <<
                " [full_output=" << FULL_OUTPUT_DEF << "]" <<
                " [decomp_by_rows=" << DECOMP_BY_ROWS_DEF << "]" <<
                " [tas_size=" << TASK_SIZE_DEF << "]" <<
                endl;
            MPI_Finalize();
            return 0;
        }
    }

    const int im = (argc > 1) ? stoi(argv[1]) : IM_DEF;
    const int n = (argc > 2) ? stoi(argv[2]) : N_DEF;
    const int km = std::pow(2, n);
    const int repeats = (argc > 3) ? stoi(argv[3]) : REPEATS_DEF;
    const bool file_output = (argc > 4) ? stoi(argv[4]) : FOUT_DEF;
    const bool full_output = (argc > 5) ? stoi(argv[5]) : FULL_OUTPUT_DEF;
    const bool decomp_by_rows = (argc > 6) ? stoi(argv[6]) : DECOMP_BY_ROWS_DEF;
    const int task_size = (argc > 7) ? stoi(argv[7]) : TASK_SIZE_DEF;

    //const size_t imp = im + 2;
    const size_t imp2 = 2 * im + 2;    

    const auto km_decomp = BlockDecomposition(km, size).grown(1); // km + 2
    const auto km4_decomp = BlockDecomposition(4 * km, size); // 4 * km
    const auto im2_decomp = BlockDecomposition(2 * im, size).grown(1); // 2 * im + 2

    const auto decomp_im2_km = decomp_by_rows ? BlockDecomposition2D(im2_decomp, km_decomp.fullSize()) : BlockDecomposition2D(im2_decomp.fullSize(), km_decomp);
    const auto decomp_im2_km4 = decomp_by_rows ? BlockDecomposition2D(im2_decomp, km4_decomp.fullSize()) : BlockDecomposition2D(im2_decomp.fullSize(), km4_decomp);

    DDArray2 jf(decomp_im2_km, rank);
    DDArray2 aa1(decomp_im2_km, rank, 1);
    DDArray2 bb(decomp_im2_km4, rank), al(decomp_im2_km4, rank, 1), be(decomp_im2_km4, rank, 1);
    DDArray2 ff(decomp_im2_km4, rank, 1);
    DDArray2 phi(decomp_im2_km, rank);

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

    auto gatherAndOutput = [rank, size, file_output, &outputArray](const std::string &header, const DDArray2 &array, const std::array<int, 4> &range) {
        if (file_output) {
            const auto arr = array.gather(0);
            outputArray(header, arr, range);
        }
    };

    double ft_time = 0;
    double sweep_time = 0;
    double comp_time = 0;
    double shadow_time = 0;
    double redistr_time = 0;

    auto syncShadows = [rank, size, &shadow_time](DDArray2 &array) {
        Timer tm;
        array.syncShadows();
        shadow_time += tm.time();
    };

    const int TAG_SHADOW_PIECE = 1000;

    auto sendShadowPieceAsync = [TAG_SHADOW_PIECE, &shadow_time](DDArray2 &src, int i, int k, int num, int dest_rank) -> AsyncOp {
        Timer tm;
        MPI_Request req;
        MPI_Isend(&src(i, k), num, MPI_DOUBLE, dest_rank, TAG_SHADOW_PIECE, MPI_COMM_WORLD, &req);
        shadow_time += tm.time();
        return AsyncOp(req);
    };

    auto recvShadowPiece = [TAG_SHADOW_PIECE, &shadow_time](DDArray2 &dest, int i, int k, int num, int src_rank) {
        Timer tm;
        MPI_Recv(&dest(i, k), num, MPI_DOUBLE, src_rank, TAG_SHADOW_PIECE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        shadow_time += tm.time();
    };

    auto redistrArray = [&redistr_time](const DDArray2 &src, DDArray2 &dst) {
        Timer tm;
        dst.combineFrom(src);
        redistr_time += tm.time();
    };

    // Init functions

    auto initTestSolution = [im, km, zm, hz, hr2, &comp_time](DDArray2 &output, int m) {
        Timer tm;
        const double a0 = -0.1;
        const double a = 1e-3;
        const double d = 1.0 + 1e-5 * m;
        const auto my_im_range = output.range(0);
        const auto i_start = my_im_range.localStart(1);
        const auto i_end = my_im_range.localEnd(2 * im + 2);
        const auto my_km_range = output.range(1);
        const auto k_start = my_km_range.localStart(0);
        const auto k_end = my_km_range.localEnd(km + 2);
        for (int k = k_start; k < k_end; k++) {
            const auto kk = my_km_range.toGlobal(k);
            const double z = hz * (kk + 1 - 1.5);
            const double z2 = z * z;
            const double s = a * z2 * (z - 1.5 - zm) + d;
            for (int i = i_start; i < i_end; i++) {
                const auto ii = my_im_range.toGlobal(i);
                output(i, k) = s * (ii - 2 * im - 1.0) * (ii - 0.5) * hr2;
            }
            //output(0, k) = -output(1, k);
            //output(2 * im + 1, k) = 0.0;
        }
        if (my_im_range.hasIndex(0)) {
            for (int k = k_start; k < k_end; k++) {
                output(0, k) = -output(1, k);
            }
        }
        if (my_im_range.hasIndex(2 * im + 1)) {
            const auto dst_loc = my_im_range.toLocal(2 * im + 1);
            for (int k = k_start; k < k_end; k++) {
                output(dst_loc, k) = 0.0;
            }
        }
        if (my_km_range.hasIndex(0)) {
            const auto dst_loc = my_km_range.toLocal(0);
            const auto src_loc = my_km_range.toLocal(1);
            for (int i = 0; i < i_end; i++) {
                output(i, dst_loc) = output(i, src_loc);
            }
        }
        if (my_km_range.hasIndex(km + 1)) {
            const auto dst_loc = my_km_range.toLocal(km + 1);
            const auto src_loc = my_km_range.toLocal(km);
            for (int i = 0; i < i_end; i++) {
                output(i, dst_loc) = output(i, src_loc);
            }
        }
        comp_time += tm.time();
    };

    auto getSolution = [rhr2, rhz2](const DDArray2 &input, int i, int k) {
        const auto ii = input.range(0).toGlobal(i);
        return (((ii + 0.5) * input(i + 1, k) - (ii - 0.5) * input(i, k)) / (ii) -
                ((ii - 0.5) * input(i, k) - (ii - 1.5) * input(i - 1, k)) / (ii - 1.0)) * rhr2 +
               (input(i, k + 1) - 2.0 * input(i, k) + input(i, k - 1)) * rhz2;
    };

    auto initTestCurrent = [im, km, rhr2, rhz2, &getSolution, &syncShadows, &comp_time](DDArray2 &input, DDArray2 &output) {
        syncShadows(input);
        const auto my_im_range = output.range(0);
        const auto i_start = my_im_range.localStart(2);
        const auto i_end = my_im_range.localEnd(2 * im + 1);
        const auto my_km_range = output.range(1);
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km + 1);
        const auto has_one = my_im_range.hasIndex(1);
        Timer tm;
        for (int k = k_start; k < k_end; k++) {
            if (has_one) {
                double s = (1.5 * input(2, k) - 4.5 * input(1, k)) * rhr2 +
                        (input(1, k + 1) - 2.0 * input(1, k) + input(1, k - 1)) * rhz2;
                output(1, k) = -s;
            }
            for (int i = i_start; i < i_end; i++) {
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

    auto computeFFTByRows = [im, km, n, &redistrArray, &comp_time](const DDArray2 &input, DDArray2 &output) {
        CArray1 dan(2 * km);
        Timer ctm;
        const auto my_im_range = input.range();
        const auto i_start = my_im_range.localStart(1);
        const auto i_end = my_im_range.localEnd(2 * im + 1);
        for (int i = i_start; i < i_end; i++) {
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
        comp_time += ctm.time();
    };

    auto computeFFTByCols = [&im2_decomp, &redistrArray, &computeFFTByRows](const DDArray2 &input, DDArray2 &output) {
        DDArray2 input_i(BlockDecomposition2D(im2_decomp, input.decomp(1).fullSize()), input.rank());
        DDArray2 output_i(BlockDecomposition2D(im2_decomp, output.decomp(1).fullSize()), output.rank());
        redistrArray(input, input_i);
        computeFFTByRows(input_i, output_i);
        redistrArray(output_i, output);
    };

    auto computeFFT = [&ft_time, &computeFFTByRows, &computeFFTByCols](const DDArray2 &input, DDArray2 &output) {
        Timer tm;
        if (input.isDistributedByRows()) {
            computeFFTByRows(input, output);
        } else {
            computeFFTByCols(input, output);
        }
        ft_time += tm.time();
    };

    auto computeSweepByCols = [im, km, imp2, hr, hr2, hz, hz2, c, &comp_time](const DDArray2 &input, DDArray2 &output) {
        DArray1 al(imp2), be(imp2);
        const auto my_km4_range = input.range();
        const auto k_start = my_km4_range.localStart(0);
        const auto k_end = my_km4_range.localEnd(4 * km);
        Timer ctm;
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
        comp_time += ctm.time();
    };

    auto computeSweepByRows = [im, km, imp2, hr, hr2, hz, hz2, c, task_size, &al, &be, &comp_time,
                                &sendShadowPieceAsync, &recvShadowPiece](const DDArray2 &input, DDArray2 &output) {
        const auto num_of_tasks = input.size(1) / task_size;
        BlockDecomposition task_decomp(input.size(1), num_of_tasks);
        //DArray2 al(input.size(0), input.size(1), 1), be(input.size(0), input.size(1), 1);
        // Front sweep
        const auto my_im_range = input.range();
        const auto i_start1 = my_im_range.localStart(2);
        const auto i_end1 = my_im_range.localEnd(2 * im + 1);
        const auto has_one = my_im_range.hasIndex(1);
        const int my_rank = input.rank();
        const int num_of_nodes = input.numOfNodes();
        AsyncOps send_ops;
        DArray1 sins(4 * km);
        for (int k = 0; k < sins.size(); k++) {
            sins[k] = sin(c * (k >= 2 * km ? k - 2 * km : k));
        }
        for (int task = 0; task < task_decomp.numOfBlocks(); task++) {
            const auto k_start = task_decomp.blockShift(task);
            const auto k_size = task_decomp.blockSize(task);
            const auto k_end = k_start + k_size;
            // Recv al, be from prev rank
            if (my_rank > 0) {
                recvShadowPiece(al, -1, k_start, k_size, my_rank - 1);
                recvShadowPiece(be, -1, k_start, k_size, my_rank - 1);
            }
            Timer ctm;
            if (has_one) {
                for (int k = k_start; k < k_end; k++) {
                    const auto dsin = sins[k];
                    const double s = 9.0 / (2.0 * hr2) + (4.0 / hz2) * dsin * dsin;
                    al(1, k) = 3.0 / (2.0 * hr2 * s);
                    be(1, k) = input(1, k) / s;
                }
            }
            for (int i = i_start1; i < i_end1; i++) {
                const auto ii = my_im_range.toGlobal(i);
                for (int k = k_start; k < k_end; k++) {
                    const auto dsin = sins[k];
                    const double s = (2.0 * ((ii - 0.5) / hr) * ((ii - 0.5) / hr)) / ((ii) * (ii - 1.0)) +
                        (4.0 / hz2) * dsin * dsin -
                        al(i - 1, k) * (ii - 1.5) / ((ii - 1.0) * hr2);
                    al(i, k) = (ii + 0.5) / (s * (ii) * hr2);
                    be(i, k) = (be(i - 1, k) * (ii - 1.5) / ((ii - 1.0) * hr2) + input(i, k)) / s;
                }
            }
            comp_time += ctm.time();
            // Send al, be to next rank
            if (my_rank < num_of_nodes - 1) {
                send_ops.add(sendShadowPieceAsync(al, al.size(0) - 1, k_start, k_size, my_rank + 1));
                send_ops.add(sendShadowPieceAsync(be, be.size(0) - 1, k_start, k_size, my_rank + 1));
            }
        }
        send_ops.wait();
        // Back sweep
        const auto i_start2 = my_im_range.localStart(1);
        const auto i_end2 = my_im_range.localEnd(2 * im + 1);
        const auto has_last = my_im_range.hasIndex(2 * im + 1);
        const auto last_ind = my_im_range.toLocal(2 * im + 1);
        for (int task = 0; task < task_decomp.numOfBlocks(); task++) {
            const auto k_start = task_decomp.blockShift(task);
            const auto k_size = task_decomp.blockSize(task);
            const auto k_end = k_start + k_size;
            // Recv output from next rank
            if (my_rank < num_of_nodes - 1) {
                recvShadowPiece(output, output.size(0), k_start, k_size, my_rank + 1);
            }
            Timer ctm;
            if (has_last) {
                for (int k = k_start; k < k_end; k++) {
                    output(last_ind, k) = 0.0;
                }
            }
            for (int i = i_end2 - 1; i >= i_start2; i--) {
                for (int k = k_start; k < k_end; k++) {
                    output(i, k) = al(i, k) * output(i + 1, k) + be(i, k);
                }
            }
            comp_time += ctm.time();
            // Send output to prev rank
            if (my_rank > 0) {
                send_ops.add(sendShadowPieceAsync(output, 0, k_start, k_size, my_rank - 1));
            }
        }
        send_ops.wait();
    };

    auto computeSweep = [&sweep_time, &computeSweepByCols, &computeSweepByRows](const DDArray2 &input, DDArray2 &output) {
        Timer tm;
        if (input.isDistributedByRows()) {
            computeSweepByRows(input, output);
        } else {
            computeSweepByCols(input, output);
        }
        sweep_time += tm.time();
    };

    auto computeFFTInverseByRows = [im, km, n, &comp_time](const DDArray2 &input, DDArray2 &output) {
        CArray1 dan(2 * km);
        const auto my_im_range = input.range(0);
        const auto i_start = my_im_range.localStart(1);
        const auto i_end = my_im_range.localEnd(2 * im + 1);
        Timer ctm;
        for (int i = i_start; i < i_end; i++) {
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
        comp_time += ctm.time();
    };

    auto computeFFTInverseByCols = [&redistrArray, &im2_decomp, &computeFFTInverseByRows](const DDArray2 &input, DDArray2 &output) {
        DDArray2 input_i(BlockDecomposition2D(im2_decomp, input.decomp(1).fullSize()), input.rank());
        DDArray2 output_i(BlockDecomposition2D(im2_decomp, output.decomp(1).fullSize()), output.rank());
        redistrArray(input, input_i);
        computeFFTInverseByRows(input_i, output_i);
        redistrArray(output_i, output);
    };

    auto computeFFTInverse = [&ft_time, &computeFFTInverseByRows, &computeFFTInverseByCols](const DDArray2 &input, DDArray2 &output) {
        Timer tm;
        if (input.isDistributedByRows()) {
            computeFFTInverseByRows(input, output);
        } else {
            computeFFTInverseByCols(input, output);
        }
        ft_time += tm.time();
    };

    auto computeSolutionDifference = [im, km, &comp_time](const DDArray2 &first, const DDArray2 &second) -> DDArray2 {
        DDArray2 output(first.decomp2d(), first.rank());
        const auto my_im_range = output.range(0);
        const auto i_start = my_im_range.localStart(2);
        const auto i_end = my_im_range.localEnd(im + 1);
        const auto my_km_range = output.range(1);
        const auto k_start = my_km_range.localStart(1);
        const auto k_end = my_km_range.localEnd(km);
        Timer tm;
        for (int i = i_start; i < i_end; i++) {
            for (int k = k_start; k < k_end; k++) {
                output(i, k) = first(i, k) - second(i, k);
            }
        }
        comp_time += tm.time();
        return output;
    };

    // The main program's body

    Timer work_timer;

    for (int m = 0; m < repeats; m++) {
        initTestSolution(aa1, m + 1);
        initTestCurrent(aa1, jf);

        computeFFT(jf, bb);
        computeSweep(bb, ff);
        computeFFTInverse(ff, phi);
    }

    gatherAndOutput("phi phi", phi, {0, 7, 0, km + 2});
    gatherAndOutput("proverka 2 dd dd", computeSolutionDifference(phi, aa1), {0, 7, 0, 6});

    auto work_time = work_timer.time();

    double time = 0;
    MPI_Reduce(&work_time, &time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (file_output && rank == 0) {
        out_lst.close();
    }

    std::ostringstream out;

    if (rank == 0) {
        out << "Im: " << im << ", Km: " << km << ", Repeats: " << repeats <<
            ", Nodes: " << size <<
            ", Decomp by rows: " << decomp_by_rows <<
            ", Task size: " << task_size << endl;
        out << "TIME: " << time << endl;
    }

    out << rank << ": Work time: " << work_time <<
        ", Comp time: " << comp_time <<
        ", Shadow time: " << shadow_time <<
        ", Redistr time: " << redistr_time <<
        ", FT: " << ft_time <<
        ", Sweep: " << sweep_time << endl;

    std::cout << out.str();

    MPI_Finalize();

    return 0;
}
