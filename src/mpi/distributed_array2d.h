#pragma once

#include "defs.h"
#include "block_decomp.h"

#include <mpi.h>

class DistributedArray2D {
public:
    DistributedArray2D(const BlockDecomposition2D &decomp, int rank, int shadow_size = 0);

    double* data() {
        return _data.data();
    }

    const double* data() const {
        return _data.data();
    }

    size_t size(size_t dim) const {
        return _data.size(dim);
    }

    size_t shadowSize(size_t dim) const {
        return _data.shadowSize(dim);
    }

    size_t fullSize(size_t dim) const {
        return _data.fullSize(dim);
    }

    size_t size() const {
        return _data.size();
    }

    template <typename Index>
    size_t at(Index x, Index y) const {
        return _data.at<Index>(x, y);
    }

    double& operator[](size_t index) {
        return _data[index];
    }

    const double& operator[](size_t index) const {
        return _data[index];
    }

    template <typename Index>
    double& operator()(Index x, Index y) {
        return _data(x, y);
    }

    template <typename Index>
    const double& operator()(Index x, Index y) const {
        return _data(x, y);
    }

    typename std::vector<double>::iterator begin() {
        return _data.begin();
    }

    typename std::vector<double>::iterator end() {
        return _data.end();
    }

    const BlockDecomposition2D& decomp2d() const {
        return _decomp;
    }

    const BlockDecomposition& decomp(size_t dim) const {
        return _decomp.decomp(dim);
    }

    const BlockDecomposition& decomp() const {
        return isDistributedByRows() ? decomp(0) : decomp(1);
    }

    const BlockDecomposition::Range& range(size_t dim) const {
        return _ranges[dim];
    }

    const BlockDecomposition::Range& range() const {
        return isDistributedByRows() ? range(0) : range(1);
    }

    int rank() const {
        return _rank;
    }

    int numOfNodes() const {
        return _num_of_nodes;
    }

    const DArray2& localArray() const {
        return _data;
    }

    DArray2& localArray() {
        return _data;
    }

    bool isDistributedByRows() const {
        return _decomp.isByRows();
    }

    bool isDistributedByCols() const {
        return _decomp.isByCols();
    }

    void syncShadows();
    void combineFrom(const DistributedArray2D &src);

    DArray2 gather(int root_rank = 0) const;

private:
    void init(size_t sx, size_t sy, int shadow_x, int shadow_y, const BlockDecomposition2D &decomp, int rank);

private:
    DArray2 _data;    
    BlockDecomposition2D _decomp;
    std::array<BlockDecomposition::Range, 2> _ranges;
    int _rank;
    int _num_of_nodes;
    MPI_Datatype _shadow_type;
};
