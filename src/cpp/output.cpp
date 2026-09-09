#include "output.h"

#include <iomanip>

using namespace std;

void print(const DArray2 &data, int i_start, int i_end, int j_start, int j_end, std::ofstream &out) {
    out << setw(7) << "";
    for (int i = i_start; i < i_end; i++) {
        out << setw(3) << i + 1;
        if (i < i_end - 1) {
            out << setw(6) << "";
        }
    }
    out << std::endl;
    for (int j = j_end - 1; j >= j_start; j--) {
        out << setw(3) << j + 1 << setw(1) << "";
        for (int i = i_start; i < i_end; i++) {
            out << setw(10) << std::fixed << setprecision(3) << data(i, j);
        }
        out << std::endl;
    }
}

void output(const std::string &header, const DArray2 &data, int i_start, int i_end, int j_start, int j_end, std::ofstream &out) {
    out << "\n";
    out << " " << header << "\n";
    print(data, i_start, i_end, j_start, j_end, out);
}

void output(const std::string &header, const DArray2 &data, const std::array<int, 4> &range, std::ofstream &out) {
    output(header, data, range[0], range[1], range[2], range[3], out);
}

void output(const std::string &header, const DArray2 &data, std::ofstream &out) {
    output(header, data, 0, data.size(0), 0, data.size(1), out);
}
