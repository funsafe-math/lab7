#include "fmt/core.h"
#include "fmt/ranges.h"
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <span>
#include <vector>

struct Matrix
{
    const size_t width_;
    const size_t height_;
    std::vector<float> values_;
    constexpr Matrix(size_t width, size_t height)
        : width_(width)
        , height_(height)
    {
        values_ = std::vector<float>(width * height, 0.0F);
    }
    constexpr std::span<float> row(size_t ix)
    {
        return std::span{values_.begin() + width_ * ix, width_};
    }

    constexpr std::span<float> operator[](size_t ix) { return row(ix); }

    constexpr void gaussian_elimination()
    {
        int up_to = std::min(width_, height_);

        for (size_t column_ix = 0; column_ix < up_to; column_ix++) {
            for (size_t row_ix = column_ix + 1; row_ix < up_to; row_ix++) {
                float m = row(row_ix)[column_ix] / row(column_ix)[column_ix];
                for (size_t i = column_ix; i < width_; i++) {
                    row(row_ix)[i] -= m * row(column_ix)[i];
                }
            }
        }
    }
};

void print_matrix(Matrix &mat)
{
    fmt::println("Matrix:");
    for (size_t i = 0; i < mat.height_; ++i) {
        fmt::println("{}", mat[i]);
    }
}

int main()
{
    std::ifstream file{"/home/tadeuszk/Dokumenty/lab7/example1.txt"};
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }
    size_t size;
    file >> size;
    std::cout << "Matrix size is " << size << '\n';
    Matrix mat{size + 1, size};
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            file >> mat[i][j];
        }
    }

    for (int i = 0; i < size; ++i) {
        file >> mat[i][size];
    }

    file.close();

    // fmt::println("Matrix: {}", mat.values_);
    print_matrix(mat);
    mat.gaussian_elimination();
    print_matrix(mat);
    // fmt::println("Matrix: {}", mat.values_);
    // fmt::println("Vec: {}", vec.values_);
}
