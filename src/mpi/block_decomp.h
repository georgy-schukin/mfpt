#pragma once

#include <vector>

class BlockDecomposition {
public:
    BlockDecomposition(int size, int num_of_parts);

    int getBlockSize(int block_index) const;
    int getBlockShift(int block_index) const;
    int numOfBlocks() const;
    int localStart(int index, int block_index) const;
    int localEnd(int index, int block_index) const;
    int toGlobal(int index, int block_index) const;

private:
    std::vector<int> sizes;
    std::vector<int> shifts;
};
