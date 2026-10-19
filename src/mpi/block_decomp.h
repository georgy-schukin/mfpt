#pragma once

#include <vector>
#include <array>

class BlockDecomposition {
public:
    class Range {
    public:
        Range() {}
        Range(int s, int e):
            start(s), end(e) {}

        int size() const;
        int localStart(int g_index) const;
        int localEnd(int g_index) const;
        int toGlobal(int l_index) const;
        int toLocal(int g_index) const;
        bool hasIndex(int g_index) const;

    public:
        int start;
        int end;
    };

public:
    BlockDecomposition() {}
    BlockDecomposition(int size, int num_of_parts);
    BlockDecomposition(const std::vector<int> &sizes);

    int blockSize(int block_index) const;
    int blockShift(int block_index) const;
    int numOfBlocks() const;
    int fullSize() const;
    int localStart(int index, int block_index) const;
    int localEnd(int index, int block_index) const;
    int toGlobal(int index, int block_index) const;

    Range range(int block_index) const;

    BlockDecomposition grown(int amount) const;
    BlockDecomposition multuplied(int amount) const;

private:
    static std::vector<int> computeSizes(int size, int num_of_parts);
    static std::vector<int> computeShifts(const std::vector<int> &sizes);

private:
    std::vector<int> sizes;
    std::vector<int> shifts;
};

class BlockDecomposition2D {
public:
    enum DecompositionType {
        BY_ROWS = 0,
        BY_COLS
    };

public:
    BlockDecomposition2D() {}
    BlockDecomposition2D(int sx, const BlockDecomposition &dy) :
        _decomps {BlockDecomposition(sx, 1), dy}, _dec_type(BY_COLS) {
    }
    BlockDecomposition2D(const BlockDecomposition &dx, int sy) :
        _decomps {dx, BlockDecomposition(sy, 1)}, _dec_type(BY_ROWS) {
    }

    const BlockDecomposition& decomp(int dim) const {
        return _decomps[dim];
    }

    int numOfBlocks(int dim) const {
        return decomp(dim).numOfBlocks();
    }

    int blockSize(int dim, int index) const {
        return decomp(dim).blockSize(index);
    }

    int blockShift(int dim, int index) const {
        return decomp(dim).blockShift(index);
    }

    BlockDecomposition::Range range(int dim, int index) const {
        return decomp(dim).range(index);
    }

    DecompositionType decompositionType() const {
        return _dec_type;
    }

    bool isByRows() const {
        return decompositionType() == BY_ROWS;
    }

    bool isByCols() const {
        return decompositionType() == BY_COLS;
    }

    static BlockDecomposition2D byRows(int sx, int sy, int size) {
        return BlockDecomposition2D(BlockDecomposition(sx, size), sy);
    }

    static BlockDecomposition2D byCols(int sx, int sy, int size) {
        return BlockDecomposition2D(sx, BlockDecomposition(sy, size));
    }

private:
    std::array<BlockDecomposition, 2> _decomps;
    DecompositionType _dec_type;
};

