#pragma once
#ifndef TRACE_HPP
#define TRACE_HPP

#include <cstddef>
#include <variant>
#include <vector>

#include <fmt/core.h>
#include <fmt/format.h>

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
 * 
 * Test idea: execute sequentially shuffling productions in each foata group
*/

namespace trace {

enum class Operations { ADD, SUBSTRACT, MULTIPLY, DIVIDE, VALUE, UNARY_SUBSTRACT };

struct Variable;

struct Expression // Binary or unary operation
{
    Operations op;
    std::variant<const Variable *, float> lhs; // Used to store literal values
    std::variant<const Variable *, float> rhs;

    constexpr Expression(std::variant<const Variable *, float> lhs_,
                         std::variant<const Variable *, float> rhs_,
                         Operations op_) noexcept
        : op{op_}
        , lhs{lhs_}
        , rhs{rhs_}
    {}
    explicit constexpr Expression() noexcept = default;
    inline constexpr static Expression makeBinOp(const Variable *lhs,
                                                 Operations op,
                                                 const Variable *rhs) noexcept
    {
        return {lhs, rhs, op};
    }

    inline constexpr static Expression makeLiteral(const float val) noexcept
    {
        return {val, val, Operations::VALUE};
    }

    constexpr bool isUnaryOp() const noexcept { return op == Operations::VALUE; }
};

struct Production
{
    const Variable *lhs;
    Expression rhs;
};

// Structs to be used in expression evaluations

struct LiteralExpression
{
    Operations op;
    float lhs;
    float rhs;

    constexpr operator Expression() const noexcept { return Expression{lhs, rhs, op}; }

    constexpr LiteralExpression(float literal) noexcept
        : op{Operations::VALUE}
        , lhs{literal}
        , rhs{literal}
    {}
};

struct VariableExpression : public Expression
{
    std::vector<Production> *log;

    constexpr VariableExpression() = delete;
    constexpr VariableExpression(std::variant<const Variable *, float> lhs_,
                                 std::variant<const Variable *, float> rhs_,
                                 Operations op_,
                                 std::vector<Production> *log_) noexcept
        : Expression(lhs_, rhs_, op_)
        , log{log_}
    {}

    inline constexpr static VariableExpression makeUnaryOp(const Variable *val,
                                                           Operations op,
                                                           std::vector<Production> *log) noexcept
    {
        return {val, val, op, log};
    }
};

struct Variable
{
    std::vector<Production> *log;

    constexpr inline Variable(std::vector<Production> *log_) noexcept
        : log{log_}
    {}

    constexpr inline Variable(const VariableExpression &expr) noexcept
        : log{expr.log}
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = expr;
        log->push_back(prod);
    }

    // Arithmetic operators with self

    inline constexpr VariableExpression operator+(const Variable &other) const noexcept
    {
        return {get_id(), other.get_id(), Operations::ADD, log};
    }

    inline constexpr VariableExpression operator-(const Variable &other) const noexcept
    {
        return {get_id(), other.get_id(), Operations::SUBSTRACT, log};
    }

    inline constexpr VariableExpression operator*(const Variable &other) const noexcept
    {
        return {get_id(), other.get_id(), Operations::MULTIPLY, log};
    }

    inline constexpr VariableExpression operator/(const Variable &other) const noexcept
    {
        return {get_id(), other.get_id(), Operations::DIVIDE, log};
    }

    // Unary operators with self
    inline constexpr VariableExpression operator-() const noexcept
    {
        return VariableExpression::makeUnaryOp(get_id(), Operations::UNARY_SUBSTRACT, log);
    }

    inline constexpr const Variable &operator+() const noexcept { return *this; }

    // Arithmetic operators with float
    // TODO: template this

    inline constexpr VariableExpression operator+(const float other) const noexcept
    {
        return {get_id(), other, Operations::ADD, log};
    }

    inline constexpr VariableExpression operator-(const float other) const noexcept
    {
        return {get_id(), other, Operations::SUBSTRACT, log};
    }

    inline constexpr VariableExpression operator*(const float other) const noexcept
    {
        return {get_id(), other, Operations::MULTIPLY, log};
    }

    inline constexpr VariableExpression operator/(const float other) const noexcept
    {
        return {get_id(), other, Operations::DIVIDE, log};
    }

    // Assignment operators with self

    inline constexpr Variable &operator+=(const Variable &other) noexcept
    {
        VariableExpression tmp = *this + other;
        *this = tmp;
        return *this;
    }

    inline constexpr Variable &operator-=(const Variable &other) noexcept
    {
        VariableExpression tmp = *this - other;
        *this = tmp;
        return *this;
    }
    inline constexpr Variable &operator*=(const Variable &other) noexcept
    {
        VariableExpression tmp = *this * other;
        *this = tmp;
        return *this;
    }
    inline constexpr Variable &operator/=(const Variable &other) noexcept
    {
        VariableExpression tmp = *this / other;
        *this = tmp;
        return *this;
    }

    // Assignment operators with Expression

    inline constexpr Variable &operator=(const VariableExpression &expr) noexcept
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = expr;
        expr.log->push_back(prod);
        return *this;
    }

    inline constexpr Variable &operator=(const float val) noexcept
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = Expression::makeLiteral(val);
        log->push_back(prod);
        return *this;
    }

    inline constexpr Variable &operator+=(const VariableExpression &expr) noexcept
    {
        Variable tmp = expr;
        *this = *this + tmp;

        return *this;
    }

    inline constexpr Variable &operator-=(const VariableExpression &expr) noexcept
    {
        Variable tmp = expr;
        *this = *this - tmp;
        return *this;
    }

    inline constexpr Variable &operator*=(const VariableExpression &expr) noexcept
    {
        Variable tmp = expr;
        *this = *this * tmp;
        return *this;
    }

    inline constexpr Variable &operator/=(const VariableExpression &expr) noexcept
    {
        Variable tmp = expr;
        *this = *this / tmp;
        return *this;
    }

    // Assignment operators with float
    inline constexpr Variable &operator+=(const float val) noexcept
    {
        VariableExpression tmp = *this + val;
        *this = tmp;
        return *this;
    }

    inline constexpr Variable &operator-=(const float val) noexcept
    {
        VariableExpression tmp = *this - val;
        *this = tmp;
        return *this;
    }
    inline constexpr Variable &operator*=(const float val) noexcept
    {
        VariableExpression tmp = *this * val;
        *this = tmp;
        return *this;
    }
    inline constexpr Variable &operator/=(const float val) noexcept
    {
        VariableExpression tmp = *this / val;
        *this = tmp;
        return *this;
    }

    inline constexpr Variable(const Variable &other) noexcept
        : log{other.log}
    {
        Production p;
        p.lhs = this;
        p.rhs.lhs = &other;
        p.rhs.op = Operations::VALUE;
        log->push_back(p);
    };
    Variable(Variable &&other) = delete;
    // inline constexpr ~Variable() {}

private:
    inline constexpr Variable create_without_logging() const noexcept
    {
        return Variable{this->log};
    }
    constexpr const Variable *get_id() const noexcept { return this; }
};

inline constexpr VariableExpression operator+(float lhs, const Variable &rhs) noexcept
{
    return {lhs, &rhs, Operations::ADD, rhs.log};
}

inline constexpr VariableExpression operator-(float lhs, const Variable &rhs) noexcept
{
    return {lhs, &rhs, Operations::SUBSTRACT, rhs.log};
}
inline constexpr VariableExpression operator*(float lhs, const Variable &rhs) noexcept
{
    return {lhs, &rhs, Operations::MULTIPLY, rhs.log};
}
inline constexpr VariableExpression operator/(float lhs, const Variable &rhs) noexcept
{
    return {lhs, &rhs, Operations::DIVIDE, rhs.log};
}

// Experimental: Expression-Variable operators allowing for ex. v0 = v1 + v2 * 2;
inline constexpr VariableExpression operator+(const Variable &lhs,
                                              const VariableExpression &rhs) noexcept
{
    Variable tmp = rhs;
    return lhs + tmp;
}

inline constexpr VariableExpression operator-(const Variable &lhs,
                                              const VariableExpression &rhs) noexcept
{
    Variable tmp = rhs;
    return lhs + tmp;
}
inline constexpr VariableExpression operator*(const Variable &lhs,
                                              const VariableExpression &rhs) noexcept
{
    Variable tmp = rhs;
    return lhs + tmp;
}
inline constexpr VariableExpression operator/(const Variable &lhs,
                                              const VariableExpression &rhs) noexcept
{
    Variable tmp = rhs;
    return lhs + tmp;
}

inline constexpr VariableExpression operator+(const VariableExpression &lhs,
                                              const Variable &rhs) noexcept
{
    Variable tmp = rhs;
    return tmp + rhs;
}

inline constexpr VariableExpression operator+(const VariableExpression &lhs, float rhs) noexcept
{
    Variable tmp = lhs;
    return tmp + rhs;
}

inline constexpr VariableExpression operator-(const VariableExpression &lhs, float rhs) noexcept
{
    Variable tmp = lhs;
    return tmp - rhs;
}

inline constexpr VariableExpression operator*(const VariableExpression &lhs, float rhs) noexcept
{
    Variable tmp = lhs;
    return tmp * rhs;
}

inline constexpr VariableExpression operator/(const VariableExpression &lhs, float rhs) noexcept
{
    Variable tmp = lhs;
    return tmp / rhs;
}

inline constexpr VariableExpression operator-(const VariableExpression &operand) noexcept
{
    Variable tmp = operand;
    return -tmp;
}

} // namespace trace

// Formatting
template<>
struct fmt::formatter<trace::Operations> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(trace::Operations op, FormatContext &ctx)
    {
        string_view name = "unknown";
        switch (op) {
        case trace::Operations::ADD:
            name = "+";
            break;
        case trace::Operations::SUBSTRACT:
            name = "-";
            break;
        case trace::Operations::MULTIPLY:
            name = "*";
            break;
        case trace::Operations::DIVIDE:
            name = "/";
            break;
        case trace::Operations::VALUE:
            name = "VALUE";
            break;
        case trace::Operations::UNARY_SUBSTRACT:
            name = "Unary substract";
            break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

template<>
struct fmt::formatter<std::variant<const trace::Variable *, float>> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const std::variant<const trace::Variable *, float> var, FormatContext &ctx)
    {
        if (std::holds_alternative<const trace::Variable *>(var)) {
            return fmt::format_to(ctx.out(), "{}", fmt::ptr(std::get<const trace::Variable *>(var)));
        }
        if (std::holds_alternative<float>(var)) {
            return fmt::format_to(ctx.out(), "{}", std::get<float>(var));
        }
        throw std::runtime_error("Unknown variant");
    }
};

template<>
struct fmt::formatter<trace::Expression> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const trace::Expression &expr, FormatContext &ctx)
    {
        switch (expr.op) {
        // case trace::Operations::UNKNOWN:
        //     throw std::runtime_error("UNKNOWN operation");
        case trace::Operations::ADD:
        case trace::Operations::SUBSTRACT:
        case trace::Operations::MULTIPLY:
        case trace::Operations::DIVIDE:
            return fmt::format_to(ctx.out(), "({} {} {})", expr.lhs, expr.op, expr.rhs);
        case trace::Operations::VALUE:
            return fmt::format_to(ctx.out(), "{}", expr.lhs);
            break;
        case trace::Operations::UNARY_SUBSTRACT:
            return fmt::format_to(ctx.out(), "(- {})", expr.lhs);
            break;
        }
        throw std::runtime_error("UNKNOWN operation");
    }
};

template<>
struct fmt::formatter<trace::Production> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const trace::Production &prod, FormatContext &ctx)
    {
        return fmt::format_to(ctx.out(), "{} := {}", fmt::ptr(prod.lhs), prod.rhs);
    }
};

#endif // TRACE_HPP
