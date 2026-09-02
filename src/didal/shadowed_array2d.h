#pragma once

#include "didal/base/array_nd.h"

#include <vector>
#include <array>
#include <cstddef>
#include <algorithm>

template <typename T>
class ShadowedArray2D {
public:
    ShadowedArray2D() {}
    ShadowedArray2D(size_t sx, size_t sy, size_t shadow_sx = 0, size_t shadow_sy = 0) :
        _size {sx, sy},
        _shadow_size {shadow_sx, shadow_sy},
        _data({sx + 2 * shadow_sx, sy + 2 * shadow_sy}, T {}) {
    }
    ShadowedArray2D(size_t sx, size_t sy, const T &value, size_t shadow_sx, size_t shadow_sy) :
        _size {sx, sy},
        _shadow_size {shadow_sx, shadow_sy},
        _data({sx + 2 * shadow_sx, sy + 2 * shadow_sy}, value) {
    }
    ShadowedArray2D(const std::array<size_t, 2> &sz, const std::array<size_t, 2> &shadow_sz = {0, 0}) :
        _size(sz),
        _shadow_size(shadow_sz),
        _data({sz[0] + 2 * shadow_sz[0], sz[1] + 2 * shadow_sz[1]}, T {}) {
    }

    ShadowedArray2D& operator==(const ShadowedArray2D &arr) = default;
    ShadowedArray2D& operator==(ShadowedArray2D &&arr) = default;

    T* data() {
        return _data.data().data();
    }

    const T* data() const {
        return _data.data().data();
    }

    size_t size(size_t dim) const {
        return _size[dim];
    }

    size_t shadowSize(size_t dim) const {
        return _shadow_size[dim];
    }

    size_t fullSize(size_t dim) const {
        return _data.size(dim);
    }

    size_t size() const {
        return _size[0] * _size[1];
    }

    template <typename Index>
    size_t at(Index x, Index y) const {
        return _data.toLinearIndex(x + _shadow_size[0], y + _shadow_size[1]);
    }

    bool isRowMajorOrder() const {
        return true;
    }

    template <typename Index>
    size_t atRaw(Index x, Index y) const {
        return _data.toLinearIndex(x, y);
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

    template <typename Index>
    T& raw(Index x, Index y) {
        return _data[atRaw(x, y)];
    }

    template <typename Index>
    const T& raw(Index x, Index y) const {
        return _data[atRaw(x, y)];
    }

    typename std::vector<T>::iterator begin() {
        return _data.data().begin();
    }

    typename std::vector<T>::iterator end() {
        return _data.data().end();
    }

private:
    std::array<size_t, 2> _size;
    std::array<size_t, 2> _shadow_size;
    typename ddl::ArrayND<T, 2> _data;
};
