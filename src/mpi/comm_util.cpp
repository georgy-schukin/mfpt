#include "comm_util.h"

#include <exception>
#include <stdexcept>

namespace {

enum Tags: int {
    TAG_PREV = 0,
    TAG_NEXT = 1,
    TAG_GATHER = 2,
    TAG_SCATTER = 3
};

}

MPI_Datatype makeColType(const DArray2 &array) {
    return makeDataVectorType(array.size(0), 1, array.fullSize(1), 1);
}

MPI_Datatype makeRowType(const DArray2 &array) {
    return makeRowType(array.size(1), array.fullSize(1));
}

MPI_Datatype makeRowType(int block_size, int extent) {
    MPI_Datatype tmp_type, row_type;
    MPI_Type_contiguous(block_size, MPI_DOUBLE, &tmp_type);
    MPI_Type_create_resized(tmp_type, 0, extent * sizeof(double), &row_type);
    MPI_Type_commit(&row_type);
    return row_type;
}

MPI_Datatype makeDataVectorType(int num_of_blocks, int block_size, int stride, int extent) {
    MPI_Datatype type, vtype;
    MPI_Type_vector(num_of_blocks, block_size, stride, MPI_DOUBLE, &type);
    if (extent > 0) {
        MPI_Type_create_resized(type, 0, extent * sizeof(double), &vtype);
    } else {
        MPI_Type_dup(type, &vtype);
    }
    MPI_Type_commit(&vtype);
    return vtype;
}

MPI_Datatype makeDataVectorType(const DArray2 &array) {
    return makeDataVectorType(array.size(0), array.size(1), array.fullSize(1));
}

void syncShadowsColsPrev(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank > 0) {
        // Recv in shadow from prev.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), size_t(0)), 1, col_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &r1);
    }
    if (rank < size - 1) {
        // Send data column to next.
        MPI_Isend(&arr(size_t(0), arr.size(1) - 1), 1, col_type, rank + 1, TAG_PREV, MPI_COMM_WORLD, &r2);
    }
    if (rank > 0) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsColsNext(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank < size - 1) {
        // Recv in shadow from next.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), arr.fullSize(1) - 1), 1, col_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &r1);
    }
    if (rank > 0) {
        // Send data column to prev.
        MPI_Isend(&arr(0, 0), 1, col_type, rank - 1, TAG_NEXT, MPI_COMM_WORLD, &r2);
    }
    if (rank < size - 1) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank > 0) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsRowsPrev(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank > 0) {
        // Recv in shadow from prev.
        MPI_Irecv(&arr.raw(size_t(0), arr.shadowSize(1)), 1, row_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &r1);
    }
    if (rank < size - 1) {
        // Send data row to next.
        MPI_Isend(&arr(arr.size(0) - 1, size_t(0)), 1, row_type, rank + 1, TAG_PREV, MPI_COMM_WORLD, &r2);
    }
    if (rank > 0) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsRowsNext(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank < size - 1) {
        // Recv in shadow from next.
        MPI_Irecv(&arr.raw(arr.fullSize(0) - 1, arr.shadowSize(1)), 1, row_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &r1);
    }
    if (rank > 0) {
        // Send data row to prev.
        MPI_Isend(&arr(0, 0), 1, row_type, rank - 1, TAG_NEXT, MPI_COMM_WORLD, &r2);
    }
    if (rank < size - 1) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank > 0) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsCols(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    syncShadowsColsPrev(arr, col_type, rank, size);
    syncShadowsColsNext(arr, col_type, rank, size);
}

void syncShadowsRows(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    syncShadowsRowsPrev(arr, row_type, rank, size);
    syncShadowsRowsNext(arr, row_type, rank, size);
}

DArray2 gatherArrayCols(const DArray2 &local_data, const BlockDecomposition &cols_decomp, int rank, int size, int root) {
    std::vector<MPI_Request> reqs;
    if (rank == root) {
        reqs.resize(size + 1);
    } else {
        reqs.resize(1);
    }

    auto send_type = makeDataVectorType(local_data);
    MPI_Isend(&local_data(0, 0), 1, send_type, 0, TAG_GATHER, MPI_COMM_WORLD, &reqs[0]);
    MPI_Type_free(&send_type);

    DArray2 data;
    if (rank == root) {
        data = DArray2(local_data.size(0), cols_decomp.fullSize());
        for (int r = 0; r < size; r++) {
            auto recv_type = makeDataVectorType(data.size(0), cols_decomp.getBlockSize(r), data.fullSize(1));
            MPI_Irecv(&data(0, cols_decomp.getBlockShift(r)), 1, recv_type, r, TAG_GATHER, MPI_COMM_WORLD, &reqs[r + 1]);
            MPI_Type_free(&recv_type);
        }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    return data;
}

DArray2 scatterArrayCols(const DArray2 &data, int size_x, const BlockDecomposition &cols_decomp, int shadow_x, int shadow_y, int rank, int size, int root) {
    std::vector<MPI_Request> reqs;
    if (rank == root) {
        reqs.resize(size + 1);
    } else {
        reqs.resize(1);
    }

    DArray2 local_data(size_x, cols_decomp.getBlockSize(rank), shadow_x, shadow_y);

    auto recv_type = makeDataVectorType(local_data);
    MPI_Irecv(&local_data(0, 0), 1, recv_type, 0, TAG_SCATTER, MPI_COMM_WORLD, &reqs[0]);
    MPI_Type_free(&recv_type);

    if (rank == root) {
        for (int r = 0; r < size; r++) {
            auto send_type = makeDataVectorType(data.size(0), cols_decomp.getBlockSize(r), data.size(1) + 2 * data.shadowSize(1));
            MPI_Isend(&data(0, cols_decomp.getBlockShift(r)), 1, send_type, r, TAG_SCATTER, MPI_COMM_WORLD, &reqs[r + 1]);
            MPI_Type_free(&send_type);
        }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    return local_data;
}

void combineRowsFromCols(const DArray2 &src, DArray2 &dst, const BlockDecomposition &rows_decomp, const BlockDecomposition &cols_decomp, int rank, int size) {
    auto send_row_type = makeRowType(src);
    std::vector<int> send_displs(size), recv_displs(size);
    std::vector<int> send_counts(size), recv_counts(size);
    std::vector<MPI_Datatype> send_types(size, send_row_type);
    std::vector<MPI_Datatype> recv_types(size);

    for (int i = 0; i < size; i++) {
        send_counts[i] = rows_decomp.getBlockSize(i);
        send_displs[i] = (src.at(rows_decomp.getBlockShift(i), 0) - src.at(0, 0)) * sizeof(double);
        recv_counts[i] = rows_decomp.getBlockSize(rank);
        recv_displs[i] = (dst.at(0, cols_decomp.getBlockShift(i)) - dst.at(0, 0)) * sizeof(double);
        recv_types[i] = makeRowType(cols_decomp.getBlockSize(i), dst.fullSize(1));
        //std::cout << rank << " RC: send " << i << ": " << send_counts[i] << " " << send_displs[i] << ", recv " << recv_counts[i] << " " << recv_displs[i] << std::endl;
    }

    MPI_Alltoallw(&src(0, 0), send_counts.data(), send_displs.data(), send_types.data(),
                  &dst(0, 0), recv_counts.data(), recv_displs.data(), recv_types.data(), MPI_COMM_WORLD);
    MPI_Type_free(&send_row_type);
    for (auto &tp: recv_types) {
        MPI_Type_free(&tp);
    }
}

void combineColsFromRows(const DArray2 &src, DArray2 &dst, const BlockDecomposition &cols_decomp, const BlockDecomposition &rows_decomp, int rank, int size) {
    auto recv_row_type = makeRowType(dst);
    std::vector<int> send_displs(size), recv_displs(size);
    std::vector<int> send_counts(size), recv_counts(size);
    std::vector<MPI_Datatype> send_types(size);
    std::vector<MPI_Datatype> recv_types(size, recv_row_type);

    for (int i = 0; i < size; i++) {
        send_counts[i] = rows_decomp.getBlockSize(rank);
        send_displs[i] = (src.at(0, cols_decomp.getBlockShift(i)) - src.at(0, 0)) * sizeof(double);
        send_types[i] = makeRowType(cols_decomp.getBlockSize(i), src.fullSize(1));
        recv_counts[i] = rows_decomp.getBlockSize(i);
        recv_displs[i] = (dst.at(rows_decomp.getBlockShift(i), 0) - dst.at(0, 0)) * sizeof(double);
        //std::cout << rank << " CR: send " << i << ": " << send_counts[i] << " " << send_displs[i] << ", recv " << recv_counts[i] << " " << recv_displs[i] << std::endl;
    }

    MPI_Alltoallw(&src(0, 0), send_counts.data(), send_displs.data(), send_types.data(),
                  &dst(0, 0), recv_counts.data(), recv_displs.data(), recv_types.data(), MPI_COMM_WORLD);
    MPI_Type_free(&recv_row_type);
    for (auto &tp: send_types) {
        MPI_Type_free(&tp);
    }
}

void reduceSumOpVector(void *in, void *inout, int *len, MPI_Datatype *dtype) {
    if (*len != 1) {
        throw std::runtime_error("Not implemented for len != 1");
    }

    int nints, naddresses, ntypes, combiner;
    MPI_Type_get_envelope(*dtype, &nints, &naddresses, &ntypes, &combiner);
    if (combiner != MPI_COMBINER_VECTOR) {
        throw std::runtime_error("Non-vector datatype");
    }

    int vecargs [nints];
    MPI_Aint vecaddrs[naddresses];
    MPI_Datatype vectypes[ntypes];
    MPI_Type_get_contents(*dtype, nints, naddresses, ntypes, vecargs, vecaddrs, vectypes);

    if (vectypes[0] != MPI_DOUBLE) {
        throw std::runtime_error("Not a vector of doubles");
    }

    int count = vecargs[0];
    int blocklen = vecargs[1];
    int stride = vecargs[2];

    double *invec = (double*)in;
    double *inoutvec = (double*)inout;
    for (int i = 0; i < count; i++) {
        const auto shift = i * stride;
        for(int j = 0; j < blocklen; j++) {
            inoutvec[shift + j] += invec[shift + j];
        }
    }
}

MPI_Op makeVectorSumOp() {
    MPI_Op op;
    MPI_Op_create(&reduceSumOpVector, 1, &op);
    return op;
}
