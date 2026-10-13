#include "distributed_array2d.h"
#include "comm_util.h"

#include <mpi.h>

DistributedArray2D::DistributedArray2D(const BlockDecomposition2D &decomp, int rank, int shadow_size) {
    if (decomp.isByRows()) {
        init(decomp.blockSize(0, rank), decomp.blockSize(1, 0), shadow_size, 0, decomp, rank);
    } else {
        init(decomp.blockSize(0, 0), decomp.blockSize(1, rank), 0, shadow_size, decomp, rank);
    }
}

void DistributedArray2D::init(size_t sx, size_t sy, int shadow_x, int shadow_y, const BlockDecomposition2D &decomp, int rank) {
    _data = std::move(DArray2(sx, sy, shadow_x, shadow_y));
    _decomp = decomp;
    _rank = rank;
    if (isDistributedByRows()) {
        _num_of_nodes = decomp.decomp(0).numOfBlocks();
        _ranges = {decomp.decomp(0).range(rank), decomp.decomp(1).range(0)};
        _shadow_type = makeRowType(_data);
    } else {
        _num_of_nodes = decomp.decomp(1).numOfBlocks();
        _ranges = {decomp.decomp(0).range(0), decomp.decomp(1).range(rank)};
        _shadow_type = makeColType(_data);
    }
}

void DistributedArray2D::syncShadows() {
    if (isDistributedByRows()) {
        syncShadowsRows(_data, _shadow_type, _rank, _num_of_nodes);
    } else {
        syncShadowsCols(_data, _shadow_type, _rank, _num_of_nodes);
    }
}

void DistributedArray2D::combineFrom(const DistributedArray2D &src) {
    if (isDistributedByCols() && src.isDistributedByRows()) {
        combineColsFromRows(src.localArray(), localArray(), decomp(), src.decomp(), rank(), numOfNodes());
    } else if (isDistributedByRows() && src.isDistributedByCols()) {
        combineRowsFromCols(src.localArray(), localArray(), decomp(), src.decomp(), rank(), numOfNodes());
    }
}

DArray2 DistributedArray2D::gather(int root_rank) const {
    if (isDistributedByRows()) {
        return gatherArrayRows(localArray(), decomp(), _rank, _num_of_nodes, root_rank);
    } else {
        return gatherArrayCols(localArray(), decomp(), _rank, _num_of_nodes, root_rank);
    }
}
