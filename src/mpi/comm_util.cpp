#include "comm_util.h"

#include <exception>

namespace {

enum Tags: int {
    TAG_PREV = 0,
    TAG_NEXT = 1,
    TAG_GATHER = 2,
    TAG_SCATTER = 3
};

}

MPI_Datatype makeColType(const DArray2 &array) {
    MPI_Datatype col_type;
    MPI_Type_vector(array.size(0), 1, array.fullSize(1), MPI_DOUBLE, &col_type);
    MPI_Type_commit(&col_type);
    return col_type;
}

MPI_Datatype makeRowType(const DArray2 &array) {
    MPI_Datatype row_type;
    MPI_Type_contiguous(array.size(1), MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);
    return row_type;
}

MPI_Datatype makeDataVectorType(int num_of_blocks, int block_size, int stride) {
    MPI_Datatype send_type;
    MPI_Type_vector(num_of_blocks, block_size, stride, MPI_DOUBLE, &send_type);
    MPI_Type_commit(&send_type);
    return send_type;
}

MPI_Datatype makeDataVectorType(const DArray2 &array) {
    return makeDataVectorType(array.size(0), array.size(1), array.size(1) + 2 * array.shadowSize(1));
}

void syncShadowsKPrev(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    if (rank > 0) {
        MPI_Request req[2];
        // Send data column.
        MPI_Isend(&arr(0, 0), 1, col_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &req[0]);
        // Receive in shadow.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), size_t(0)), 1, col_type, rank - 1, TAG_NEXT, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    }
}

void syncShadowsKNext(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    if (rank < size - 1) {
        MPI_Request req[2];
        // Send data column.
        MPI_Isend(&arr(size_t(0), arr.size(1) - 1), 1, col_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &req[0]);
        // Receive in shadow.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), arr.fullSize(1) - 1), 1, col_type, rank + 1, TAG_PREV, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    }
}

void syncShadowsIPrev(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    if (rank > 0) {
        MPI_Request req[2];
        // Send data column.
        MPI_Isend(&arr(0, 0), 1, row_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &req[0]);
        // Receive in shadow.
        MPI_Irecv(&arr.raw(size_t(0), arr.shadowSize(1)), 1, row_type, rank - 1, TAG_NEXT, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    }
}

void syncShadowsINext(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    if (rank < size - 1) {
        MPI_Request req[2];
        // Send data column.
        MPI_Isend(&arr(arr.size(0) - 1, size_t(0)), 1, row_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &req[0]);
        // Receive in shadow.
        MPI_Irecv(&arr.raw(arr.fullSize(0) - 1, arr.shadowSize(1)), 1, row_type, rank + 1, TAG_PREV, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
    }
}

void syncShadowsK(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    syncShadowsKPrev(arr, col_type, rank, size);
    syncShadowsKNext(arr, col_type, rank, size);
}

void syncShadowsI(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    syncShadowsIPrev(arr, row_type, rank, size);
    syncShadowsINext(arr, row_type, rank, size);
}

DArray2 gatherArrayK(const DArray2 &local_data, const BlockDecomposition &k_decomp, int rank, int size, int root) {
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
        data == DArray2(local_data.size(0), k_decomp.fullSize());
        for (int r = 0; r < size; r++) {
            auto recv_type = makeDataVectorType(data.size(0), k_decomp.getBlockSize(r), data.size(1) + 2 * data.shadowSize(1));
            MPI_Irecv(&data(0, k_decomp.getBlockShift(r)), 1, recv_type, r, TAG_GATHER, MPI_COMM_WORLD, &reqs[r + 1]);
            MPI_Type_free(&recv_type);
        }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    return data;
}

DArray2 scatterArrayK(const DArray2 &data, int size_x, const BlockDecomposition &k_decomp, int shadow_x, int shadow_y, int rank, int size, int root) {
    std::vector<MPI_Request> reqs;
    if (rank == root) {
        reqs.resize(size + 1);
    } else {
        reqs.resize(1);
    }

    DArray2 local_data(size_x, k_decomp.getBlockSize(rank), shadow_x, shadow_y);

    auto recv_type = makeDataVectorType(local_data);
    MPI_Irecv(&local_data(0, 0), 1, recv_type, 0, TAG_SCATTER, MPI_COMM_WORLD, &reqs[0]);
    MPI_Type_free(&recv_type);

    if (rank == root) {
        for (int r = 0; r < size; r++) {
            auto send_type = makeDataVectorType(data.size(0), k_decomp.getBlockSize(r), data.size(1) + 2 * data.shadowSize(1));
            MPI_Isend(&data(0, k_decomp.getBlockShift(r)), 1, send_type, r, TAG_SCATTER, MPI_COMM_WORLD, &reqs[r + 1]);
            MPI_Type_free(&send_type);
        }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    return local_data;
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
