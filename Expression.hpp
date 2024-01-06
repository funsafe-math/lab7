#pragma once
#include "fmt/core.h"
#include "fmt/format.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <type_traits>
#include <variant>
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

enum class Operations { UNKNOWN, VALUE, ADD, MULTIPLY, SUBSTRACT, DIVIDE, DESTROY };

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
        case Operations::DESTROY:
            name = "destroy";
            break;
        case Operations::VALUE:
            name = "value";
            break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

struct Element;

struct Expression
{
    const Element *lhs;
    std::variant<const Element *, float> rhs;
    Operations op;
};

template<>
struct fmt::formatter<Expression> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const Expression &expr, FormatContext &ctx)
    {
        if (std::holds_alternative<float>(expr.rhs)) {
            float f = std::get<float>(expr.rhs);
            return format_to(ctx.out(), "{} {} (float){}", fmt::ptr(expr.lhs), expr.op, f);
        }
        const Element *elem_ptr = std::get<const Element *>(expr.rhs);
        return format_to(ctx.out(), "{} {} {}", fmt::ptr(expr.lhs), expr.op, fmt::ptr(elem_ptr));
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
        if (op.rhs.op == Operations::DESTROY) {
            return format_to(ctx.out(), "{} DESTROY", fmt::ptr(op.lhs));
        }
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

    Element &operator+=(const Expression &expr)
    {
        // TODO: improve this
        Element tmp = *this;
        tmp = expr;
        *this = tmp + *this;
        return *this;
    }

    Element &operator-=(const Expression &expr)
    {
        // TODO: improve this
        Element tmp = *this;
        tmp = expr;
        *this = tmp - *this;
        return *this;
    }

    Expression operator+(const float other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::ADD;
        expr.rhs = other;
        return expr;
    }

    Expression operator*(const float other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::MULTIPLY;
        expr.rhs = other;
        return expr;
    }

    Expression operator/(const float other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::DIVIDE;
        expr.rhs = other;
        return expr;
    }

    Expression operator-(const float other)
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::SUBSTRACT;
        expr.rhs = other;
        return expr;
    }

    Element(const Element &other)
        : log{other.log}
    {
        BinOp op;
        op.lhs = this;
        op.rhs.lhs = &other;
        op.rhs.op = Operations::VALUE;
        log->push_back(op);
    };
    Element(Element &&other) = delete;
    // Element(const ElementInitializer &other)
    //     : log{other.log}
    // {}
    ~Element()
    {
        BinOp op;
        op.lhs = this;
        op.rhs.op = Operations::DESTROY;
        log->push_back(op);
    }
};