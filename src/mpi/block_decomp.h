#pragma once

#include <vector>

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
    BlockDecomposition(int size, int num_of_parts);
    BlockDecomposition(const std::vector<int> &sizes);

    int getBlockSize(int block_index) const;
    int getBlockShift(int block_index) const;
    int numOfBlocks() const;
    int fullSize() const;
    int localStart(int index, int block_index) const;
    int localEnd(int index, int block_index) const;
    int toGlobal(int index, int block_index) const;

    Range getRange(int block_index) const;

    BlockDecomposition grown(int amount) const;
    BlockDecomposition multuplied(int amount) const;

private:
    static std::vector<int> computeSizes(int size, int num_of_parts);
    static std::vector<int> computeShifts(const std::vector<int> &sizes);

private:
    std::vector<int> sizes;
    std::vector<int> shifts;
};
