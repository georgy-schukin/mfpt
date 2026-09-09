#pragma once

#include <chrono>

class Timer {
public:
    Timer() {
        reset();
    }

    void reset() {
        start = std::chrono::steady_clock::now();
    }

    double time() const {
        const auto end = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(end - start).count();
    }
private:
     std::chrono::time_point<std::chrono::steady_clock> start;
};
