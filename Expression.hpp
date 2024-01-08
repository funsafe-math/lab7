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

namespace trace {
enum class Operations { UNKNOWN, VALUE, ADD, MULTIPLY, SUBSTRACT, DIVIDE, DESTROY };

struct Element;

struct Expression
{
    const Element *lhs;
    std::variant<const Element *, float> rhs;
    Operations op;
};


struct BinOp
{
    const Element *lhs;
    Expression rhs;
};


struct Element
{
    std::vector<BinOp> *log;

    constexpr auto get_id() const { return this; }

    operator Expression()
    {
        Expression expr;
        expr.lhs = get_id();
        expr.op = Operations::VALUE;
        return expr;
    }

    Element &operator=(const Expression &expr)
    {
        BinOp op;
        op.lhs = get_id();
        op.rhs = expr;
        log->push_back(op);
        return *this;
    }

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


    Element &operator+=(const Expression &expr)
    {
        // TODO: improve this
        Element tmp = create_without_logging();
        tmp = expr;
        *this = *this + tmp;
        return *this;
    }

    Element &operator-=(const Expression &expr)
    {
        // TODO: improve this
        Element tmp = create_without_logging();
        tmp = expr;
        *this = *this - tmp;
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
    // Element(Element &&other) = delete;
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
    // Element create_instance() { return create_without_logging(); }

private:
    Element create_without_logging() { return Element{this->log}; }
};

} // namespace trace

template<>
struct fmt::formatter<trace::Operations> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(trace::Operations operation, FormatContext &ctx)
    {
        string_view name = "unknown";
        switch (operation) {
        case trace::Operations::UNKNOWN:
            name = "unknown";
            break;
        case trace::Operations::ADD:
            name = "+";
            break;
        case trace::Operations::MULTIPLY:
            name = "*";
            break;
        case trace::Operations::SUBSTRACT:
            name = "-";
            break;
        case trace::Operations::DIVIDE:
            name = "/";
            break;
        case trace::Operations::DESTROY:
            name = "destroy";
            break;
        case trace::Operations::VALUE:
            name = "value";
            break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

template<>
struct fmt::formatter<trace::Expression> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const trace::Expression &expr, FormatContext &ctx)
    {
        if (std::holds_alternative<float>(expr.rhs)) {
            float f = std::get<float>(expr.rhs);
            return format_to(ctx.out(), "{} {} (float){}", fmt::ptr(expr.lhs), expr.op, f);
        }
        const trace::Element *elem_ptr = std::get<const trace::Element *>(expr.rhs);
        return format_to(ctx.out(), "{} {} {}", fmt::ptr(expr.lhs), expr.op, fmt::ptr(elem_ptr));
    }
};

template<>
struct fmt::formatter<trace::BinOp> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const trace::BinOp &op, FormatContext &ctx)
    {
        if (op.rhs.op == trace::Operations::DESTROY) {
            return format_to(ctx.out(), "{} DESTROY", fmt::ptr(op.lhs));
        }
        return format_to(ctx.out(), "{} = {}", fmt::ptr(op.lhs), op.rhs);
    }
};
