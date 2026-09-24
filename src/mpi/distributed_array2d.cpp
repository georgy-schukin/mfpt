#include "distributed_array2d.h"
#include "comm_util.h"

#include <mpi.h>

DistributedArray2D::DistributedArray2D(ShadowedArray2D<double> &&data, const DistributionType &dtype, const BlockDecomposition &decomp, int rank) :
    _data(std::move(data)),
    _distr_type(dtype),
    _decomp(decomp),
    _rank(rank) {
    _num_of_nodes = decomp.numOfBlocks();
    _range = decomp.getRange(rank);
}

DistributedArray2D::DistributedArray2D(size_t sx, size_t sy, int shadow_x, int shadow_y, const DistributionType &dtype, const BlockDecomposition &decomp, int rank) :
    _data(sx, sy, shadow_x, shadow_y),
    _distr_type(dtype),
    _decomp(decomp),
    _rank(rank) {
    _num_of_nodes = decomp.numOfBlocks();
    _range = decomp.getRange(rank);
}

void DistributedArray2D::syncShadows() {
    if (distributionType() == BY_ROWS) {
        auto row_type = makeRowType(_data);
        syncShadowsI(_data, row_type, _rank, _num_of_nodes);
        MPI_Type_free(&row_type);
    } else if (distributionType() == BY_COLS) {
        auto col_type = makeColType(_data);
        syncShadowsK(_data, col_type, _rank, _num_of_nodes);
        MPI_Type_free(&col_type);
    }
}

void DistributedArray2D::combineFrom(const DistributedArray2D &src) {
    if (distributionType() == BY_COLS && src.distributionType() == BY_ROWS) {
        combineKFromI(src.local_data(), local_data(), decomp(), src.decomp(), rank(), numOfNodes());
    } else if (distributionType() == BY_ROWS && src.distributionType() == BY_COLS) {
        combineIFromK(src.local_data(), local_data(), decomp(), src.decomp(), rank(), numOfNodes());
    }
}
