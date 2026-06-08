#pragma once

#include <vector>
#include <array>
#include <cstddef>

template <typename T>
class Array2D {
public:
    Array2D() {}
    Array2D(size_t sx, size_t sy) :
        _size {sx, sy},
        _data(sx * sy) {
    }
    Array2D(const std::array<size_t, 2> &sz) :
        _size(sz),
        _data(sz[0] * sz[1]) {
    }

    void populate(const T* raw_data, size_t data_sz) {
        for (size_t i = 0; i < data_sz; i++) {
            _data[i] = raw_data[i];
        }
    }

    T* data() {
        return _data.data();
    }

    const T* data() const {
        return _data.data();
    }

    size_t size(size_t dim) const {
        return _size[dim];
    }

    size_t size() const {
        return _data.size();
    }

    template <typename Index>
    size_t at(Index x, Index y) const {
        return x * _size[1] + y;
    }

    T& operator[](size_t index) {
        return _data[index];
    }

    const T& operator[](size_t index) const {
        return _data[index];
    }

    template <typename Index>
    T& operator()(Index x, Index y) {
        return _data[at(x, y)];
    }

    template <typename Index>
    const T& operator()(Index x, Index y) const {
        return _data[at(x, y)];
    }

    typename std::vector<T>::iterator begin() {
        return _data.begin();
    }

    typename std::vector<T>::iterator end() {
        return _data.end();
    }

private:
    std::array<size_t, 2> _size;
    typename std::vector<T> _data;
};
