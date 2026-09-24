#pragma once

#include "shadowed_array2d.h"
#include "block_decomp.h"

class DistributedArray2D {
public:
    enum DistributionType {
        BY_ROWS = 0,
        BY_COLS
    };

public:
    DistributedArray2D(const BlockDecomposition &dx, size_t sy, int rank, int shadow_size = 0) :
        DistributedArray2D(dx.getBlockSize(rank), sy, shadow_size, 0, BY_ROWS, dx, rank) {
    }

    DistributedArray2D(size_t sx, const BlockDecomposition &dy, int rank, int shadow_size = 0) :
        DistributedArray2D(sx, dy.getBlockSize(rank), 0, shadow_size, BY_COLS, dy, rank) {
    }

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

    DistributionType distributionType() const {
        return _distr_type;
    }

    const BlockDecomposition& decomp() const {
        return _decomp;
    }

    const BlockDecomposition::Range& range() const {
        return _range;
    }

    int rank() const {
        return _rank;
    }

    int numOfNodes() const {
        return _num_of_nodes;
    }

    const ShadowedArray2D<double>& local_data() const {
        return _data;
    }

    void syncShadows();
    void combineFrom(DistributedArray2D &src);

private:
    DistributedArray2D(size_t sx, size_t sy, int shadow_x, int shadow_y, const DistributionType &dtype, const BlockDecomposition &decomp, int rank);

private:
    ShadowedArray2D<double> _data;
    const DistributionType _distr_type;
    const BlockDecomposition &_decomp;
    BlockDecomposition::Range _range;
    int _rank;
    int _num_of_nodes;
};
