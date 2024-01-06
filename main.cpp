#include "fmt/core.h"
#include "fmt/ranges.h"
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <span>
#include <vector>

enum class Operations {
    ADDEQ,
    DIV,

};

struct Element;

struct Expression
{
    Operations op;
    Element *lhs;
    Element *rhs;
};

struct Production
{
    Element *lhs;
    Expression rhs;
};

struct Element
{
    constexpr Element() = default;

    constexpr Element &operator+=(Element &other)
    {
        Production p{};
        p.lhs = this;
        p.rhs.lhs = this;
        p.rhs.rhs = &other;

        p.rhs.op = Operations::ADDEQ;
        ops.push_back(p);

        return *this;
    }

    constexpr Element operator/(Element &other)
    {
        Production p{};
        p.lhs = this;
        p.rhs.lhs = this;
        p.rhs.rhs = &other;

        p.rhs.op = Operations::DIV;
        ops.push_back(p);

        return *this;
    }
    constexpr Element &operator*(Element &other) { return *this; }
    constexpr Element &operator-=(Element &other) { return *this; }
    static std::vector<Production> ops;
};

std::vector<Production> Element::ops{};

template<typename T = float>
struct Matrix
{
    const size_t width_;
    const size_t height_;
    std::vector<T> values_;
    constexpr Matrix(size_t width, size_t height)
        : width_(width)
        , height_(height)
    {
        values_ = std::vector<T>(width * height, T{});
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
                T m = row(row_ix)[column_ix] / row(column_ix)[column_ix];
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

constexpr Matrix<float> generate_matrix()
{
    Matrix<float> ret{4, 3};
    ret.values_ = {2, 1, 3, 6, 4, 3, 8, 15, 6, 5, 16, 2};
    return ret;
}

constexpr auto solve_at_compile_time()
{
    Matrix mat = generate_matrix();
    mat.gaussian_elimination();
    std::array<float, 4 * 3> arr{};
    std::copy(mat.values_.begin(), mat.values_.end(), arr.begin());
    return arr;
}

constexpr auto solve_at_compile_time2()
{
    Matrix<Element> mat{4, 3};
    mat.gaussian_elimination();
    // constexpr size_t size = Element::ops.size();
    return Element::ops;
}

int main()
{
#if 0
    std::ifstream file{"/home/tadeuszk/Dokumenty/lab7/example1.txt"};
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }
    size_t size;
    file >> size;
    std::cout << "Matrix size is " << size << '\n';
    Matrix<float> mat{size + 1, size};
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            file >> mat[i][j];
        }
    }

    for (int i = 0; i < size; ++i) {
        file >> mat[i][size];
    }

    file.close();

    fmt::println("Matrix: {}", mat.values_);
    print_matrix(mat);
    mat.gaussian_elimination();
    print_matrix(mat);
#endif
    // fmt::println("Matrix: {}", mat.values_);
    fmt::println("Mat: {}", solve_at_compile_time());
    // fmt::println("Mat: {}", solve_at_compile_time2());
    for (const auto &p : solve_at_compile_time2()) {
        std::cout << p.lhs << '\n';
    }
}
