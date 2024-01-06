#pragma once

#include "fmt/core.h"
#include <cstddef>
#include <span>
#include <vector>

template<typename T = float>
struct Matrix
{
    const size_t width_;
    const size_t height_;
    std::vector<T> values_;
    T initial_value_;
    constexpr Matrix(size_t width, size_t height, T initial_value)
        : width_(width)
        , height_(height)
        , initial_value_{initial_value}
    {
        values_ = std::vector<T>(width * height, initial_value);
    }
    constexpr std::span<T> row(size_t ix)
    {
        return std::span{values_.begin() + width_ * ix, width_};
    }

    constexpr std::span<T> operator[](size_t ix) { return row(ix); }

    constexpr void gaussian_elimination()
    {
        int up_to = std::min(width_, height_);

        for (size_t column_ix = 0; column_ix < up_to; column_ix++) {
            for (size_t row_ix = column_ix + 1; row_ix < up_to; row_ix++) {
                T m = initial_value_;
                m = row(row_ix)[column_ix] / row(column_ix)[column_ix];
                for (size_t i = column_ix; i < width_; i++) {
                    row(row_ix)[i] -= m * row(column_ix)[i];
                }
            }
        }
    }
};

template<typename T>
void print_matrix(Matrix<T> &mat)
{
    fmt::println("Matrix:");
    for (size_t i = 0; i < mat.height_; ++i) {
        fmt::println("{}", mat[i]);
    }
}
