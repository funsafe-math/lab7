#pragma once

#include "fmt/core.h"
#include "fmt/ranges.h"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <span>
#include <vector>

template<typename T = float>
struct Matrix
{
    const size_t width_;
    const size_t height_;
    std::vector<T> values_;
    constexpr Matrix(size_t width, size_t height, T initial_value)
        : width_(width)
        , height_(height)
    {
        values_ = std::vector<T>(width * height, initial_value);
    }
    constexpr std::span<T> row(size_t ix)
    {
        return std::span{values_.begin() + width_ * ix, width_};
    }

    constexpr const std::span<const T> row(size_t ix) const
    {
        return std::span{values_.begin() + width_ * ix, width_};
    }

    constexpr std::span<T> operator[](size_t ix) { return row(ix); }
    constexpr std::span<const T> operator[](size_t ix) const { return row(ix); }

    constexpr void gaussian_elimination()
    {
        int up_to = std::min(width_, height_);

        // Create triangular matrix
        for (size_t column_ix = 0; column_ix < up_to; column_ix++) {
            for (size_t row_ix = column_ix + 1; row_ix < up_to; row_ix++) {
                T m = row(row_ix)[column_ix] / row(column_ix)[column_ix];
                for (size_t i = column_ix; i < width_; i++) {
                    row(row_ix)[i] -= m * row(column_ix)[i];
                }
            }
        }

        // Create diagonal normalized matrix

        for (size_t row_ix = height_ - 1; row_ix < height_; row_ix--) {
            for (int column_ix = row_ix + 1; column_ix < height_; ++column_ix) {
                T x = row(row_ix)[column_ix];
                for (size_t i = column_ix; i < width_; i++) {
                    row(row_ix)[i] -= x * row(column_ix)[i];
                }
            }

            T m = row(row_ix)[row_ix];
            for (size_t i = row_ix; i < width_; i++) {
                row(row_ix)[i] /= m;
            }
        }
    }
};

template<typename T>
void print_matrix(const Matrix<T> &mat)
{
    for (size_t i = 0; i < mat.height_; ++i) {
        fmt::println("{:: 04.2f}", mat[i]);
    }
}

static Matrix<float> from_file(const char *filename)
{
    std::ifstream file{filename};
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }
    size_t size;
    file >> size;
    std::cout << "Matrix size is " << size << '\n';
    auto width_ = size + 1;
    auto height_ = size;
    Matrix<float> ret{width_, height_, float{}};
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            file >> ret[i][j];
        }
    }

    for (int i = 0; i < size; ++i) {
        file >> ret[i][size];
    }

    file.close();
    return ret;
}
