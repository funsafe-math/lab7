#include "trace.hpp"
#include <span>

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

void test_with_matrix()
{
    std::vector<trace::Production> log;
    Matrix<trace::Variable> mat{6, 5, trace::Variable{&log}};
    mat.gaussian_elimination();

    for (auto &v : log) {
        fmt::println("{}", v);
    }
}

void test1()
{
    std::vector<trace::Production> log;
    std::vector<trace::Variable> vec{10, trace::Variable{&log}};
    auto tmp1 = -(vec[2] + 1);
    auto tmp2 = 5 * tmp1 * 2;
    auto tmp3 = tmp2 - 4;
    vec[0] = tmp3;
    vec[0] = 1.12;

    for (const auto &p : log) {
        fmt::println("{}", p);
    }
}

int main(int argc, char *argv[])
{
    test1();
    // test_with_matrix();
    return 0;
}
