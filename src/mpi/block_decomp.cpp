#include "block_decomp.h"

int BlockDecomposition::Range::size() const {
    return end - start;
}

int BlockDecomposition::Range::localStart(int g_index) const {
    if (g_index > start) {
        return (g_index < end ? g_index - start : end - start);
    } else {
        return 0;
    }
}

int BlockDecomposition::Range::localEnd(int g_index) const {
    if (g_index < end) {
        return (g_index > start ? g_index - start : 0);
    } else {
        return end - start;
    }
}

int BlockDecomposition::Range::toGlobal(int l_index) const {
    return l_index + start;
}

BlockDecomposition::BlockDecomposition(int size, int num_of_parts) {
    const int block_size = size / num_of_parts;
    for (int i = 0; i < num_of_parts; i++) {
        const int bsize = (i < size % num_of_parts) ? block_size + 1 : block_size;
        sizes.push_back(bsize);
    }
    shifts.push_back(0);
    for (const auto &bs: sizes) {
        shifts.push_back(shifts.back() + bs);
    }
}

int BlockDecomposition::getBlockSize(int block_index) const {
    return sizes[block_index];
}

int BlockDecomposition::getBlockShift(int block_index) const {
    return shifts[block_index];
}

int BlockDecomposition::numOfBlocks() const {
    return sizes.size();
}

int BlockDecomposition::localStart(int index, int block_index) const {
    return (block_index == 0 ? index : 0);
}

int BlockDecomposition::localEnd(int index, int block_index) const {
    return (block_index == numOfBlocks() - 1 ? index - getBlockShift(block_index): getBlockSize(block_index));
}

int BlockDecomposition::toGlobal(int index, int block_index) const {
    return index + getBlockShift(block_index);
}

BlockDecomposition::Range BlockDecomposition::getRange(int block_index) const {
    return BlockDecomposition::Range(getBlockShift(block_index), getBlockShift(block_index) + getBlockSize(block_index));
}
