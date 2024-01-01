// #include <blaze/Math.h>
#include "fmt/core.h"
#include "fmt/format.h"
#include <array>
#include <iostream>
#include <type_traits>
#include <vector>

/*
 * Goals:
 * - calculate max amount of threads that can be used for the computation
 * - calculate 2D pseudo-matrix of calculations to perform
 * 
 * 2 versions:
 * - dynamic
 * - compile-time-known
 * 
 * Possible optimizations:
 * - When running on a platform supporting SIMD (like modern x86), in a foata layer, group tasks to use SIMD functionality
 * 
 * Look at https://godbolt.org/z/zdEE3fr7o
*/

// Tracing test
enum class Operations { UNKNOWN, ADD, MULTIPLY, SUBSTRACT, DIVIDE };

template<>
struct fmt::formatter<Operations> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(Operations operation, FormatContext &ctx)
    {
        string_view name = "unknown";
        switch (operation) {
        case Operations::UNKNOWN:
            name = "unknown";
            break;
        case Operations::ADD:
            name = "+";
            break;
        case Operations::MULTIPLY:
            name = "*";
            break;
        case Operations::SUBSTRACT:
            name = "-";
            break;
        case Operations::DIVIDE:
            name = "/";
            break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

struct Element;

struct Expression
{
    const Element *lhs;
    const Element *rhs;
    Operations op;
};

template<>
struct fmt::formatter<Expression> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const Expression &expr, FormatContext &ctx)
    {
        return format_to(ctx.out(), "{} {} {}", fmt::ptr(expr.lhs), expr.op, fmt::ptr(expr.rhs));
    }
};

struct BinOp
{
    const Element *lhs;
    Expression rhs;
};

template<>
struct fmt::formatter<BinOp> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const BinOp &op, FormatContext &ctx)
    {
        return format_to(ctx.out(), "{} = {}", fmt::ptr(op.lhs), op.rhs);
    }
};

struct Element
{
    std::vector<BinOp> *log;

    constexpr auto get_id() const { return this; }
    Expression operator+(const Element &other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::ADD;
        expr.rhs = other.get_id();
        return expr;
    }

    Expression operator*(const Element &other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::MULTIPLY;
        expr.rhs = other.get_id();
        return expr;
    }

    Expression operator/(const Element &other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::DIVIDE;
        expr.rhs = other.get_id();
        return expr;
    }

    Expression operator-(const Element &other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::SUBSTRACT;
        expr.rhs = other.get_id();
        return expr;
    }

    Element(std::vector<BinOp> *log_)
        : log{log_}
    {}

    Element &operator=(const Expression &expr)
    {
        BinOp op;
        op.lhs = get_id();
        op.rhs = expr;
        log->push_back(op);
        return *this;
    }
};

int main()
{
    std::vector<BinOp> log{};
    std::vector<Element> vec{100, Element(&log)};

    vec[0] = vec[1] + vec[0];

    for (auto &v : log) {
        fmt::println("{}", v);
    }
}
