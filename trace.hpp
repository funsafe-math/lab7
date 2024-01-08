#pragma once

#include <cstddef>
#include <fmt/core.h>
#include <fmt/format.h>

// #include <fmt/std.h>
#include <optional>
#include <variant>
#include <vector>

namespace trace {

enum class Operations { UNKNOWN, ADD, SUBSTRACT, MULTIPLY, DIVIDE, VALUE, UNARY_SUBSTRACT };

struct Variable;

struct Expression // Binary or unary operation
{
    Operations op;
    std::variant<std::monostate, const Variable *, float> lhs; // Used to store literal values
    std::variant<std::monostate, const Variable *, float> rhs;

    Expression(std::variant<std::monostate, const Variable *, float> lhs_,
               std::variant<std::monostate, const Variable *, float> rhs_,
               Operations op_)
        : op{op_}
        , lhs{lhs_}
        , rhs{rhs_}
    {}
    explicit Expression() = default;
    inline static Expression makeBinOp(const Variable *lhs, Operations op, const Variable *rhs)
    {
        return {lhs, rhs, op};
    }

    // inline static Expression makeUnaryOp(const Variable *val, Operations op)
    // {
    //     Expression ret;
    //     ret.lhs = val;
    //     ret.op = op;
    //     return ret;
    // }
    inline static Expression makeLiteral(const float val) { return {val, val, Operations::VALUE}; }

    bool isUnaryOp() const { return op == Operations::VALUE; }

    // Please do not call unless you know what you're doing
    // std::optional<const Variable *> any_variable() const
    // {
    //     if (std::holds_alternative<const Variable *>(lhs)) {
    //         return std::get<const Variable *>(lhs);
    //     }
    //     if (std::holds_alternative<const Variable *>(rhs)) {
    //         return std::get<const Variable *>(rhs);
    //     }
    //     return {};
    // }

    // Evaluate if both operands are literal types
    // std::optional<float> evaluate() const
    // {
    //     if (any_variable().has_value())
    //         return {};
    //     float lhs = std::get<float>(this->lhs);
    //     float rhs = std::get<float>(this->rhs);
    //     switch (op) {
    //     case Operations::UNKNOWN:
    //         return {};
    //     case Operations::ADD:
    //         return (lhs + rhs);
    //     case Operations::SUBSTRACT:
    //         return (lhs - rhs);
    //     case Operations::MULTIPLY:
    //         return (lhs * rhs);
    //     case Operations::DIVIDE:
    //         return (lhs / rhs);
    //     case Operations::VALUE:
    //         return (lhs);
    //     }
    // }
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

    operator Expression() const { return Expression{lhs, rhs, op}; }

    LiteralExpression(float literal)
        : op{Operations::VALUE}
        , lhs{literal}
        , rhs{literal}
    {}
};

struct VariableExpression : public Expression
{
    std::vector<Production> *log;

    VariableExpression() = delete;
    VariableExpression(std::variant<std::monostate, const Variable *, float> lhs_,
                       std::variant<std::monostate, const Variable *, float> rhs_,
                       Operations op_,
                       std::vector<Production> *log_)
        : Expression(lhs_, rhs_, op_)
        , log{log_}
    {}

    inline static VariableExpression makeUnaryOp(const Variable *val,
                                                 Operations op,
                                                 std::vector<Production> *log)
    {
        return {val, val, op, log};
    }
};

struct Variable
{
    std::vector<Production> *log;

    Variable(std::vector<Production> *log_)
        : log{log_}
    {}

    Variable(const VariableExpression &expr)
        : log{expr.log}
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = expr;
        log->push_back(prod);
    }

    // Arithmetic operators with self

    VariableExpression operator+(const Variable &other) const
    {
        return {get_id(), other.get_id(), Operations::ADD, log};
    }

    VariableExpression operator-(const Variable &other) const
    {
        return {get_id(), other.get_id(), Operations::SUBSTRACT, log};
    }

    VariableExpression operator*(const Variable &other) const
    {
        return {get_id(), other.get_id(), Operations::MULTIPLY, log};
    }

    VariableExpression operator/(const Variable &other) const
    {
        return {get_id(), other.get_id(), Operations::DIVIDE, log};
    }

    // Unary operators with self
    VariableExpression operator-() const
    {
        return VariableExpression::makeUnaryOp(get_id(), Operations::UNARY_SUBSTRACT, log);
    }

    const Variable &operator+() const { return *this; }

    // Arithmetic operators with float
    // TODO: template this

    VariableExpression operator+(const float other) const
    {
        return {get_id(), other, Operations::ADD, log};
    }

    VariableExpression operator-(const float other) const
    {
        return {get_id(), other, Operations::SUBSTRACT, log};
    }

    VariableExpression operator*(const float other) const
    {
        return {get_id(), other, Operations::MULTIPLY, log};
    }

    VariableExpression operator/(const float other) const
    {
        return {get_id(), other, Operations::DIVIDE, log};
    }

    // Assignment operators

    Variable &operator=(const VariableExpression &expr)
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = expr;
        expr.log->push_back(prod);
        return *this;
    }

    Variable &operator=(const float val)
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = Expression::makeLiteral(val);
        log->push_back(prod);
        return *this;
    }

    Variable &operator+=(const VariableExpression &expr)
    {
        Variable tmp = expr;
        *this = *this + tmp;

        return *this;
    }

    Variable &operator-=(const VariableExpression &expr)
    {
        Variable tmp = expr;
        *this = *this - tmp;
        return *this;
    }

    Variable &operator*=(const VariableExpression &expr)
    {
        Variable tmp = expr;
        *this = *this * tmp;
        return *this;
    }

    Variable &operator/=(const VariableExpression &expr)
    {
        Variable tmp = expr;
        *this = *this / tmp;
        return *this;
    }

    Variable(const Variable &other)
        : log{other.log}
    {
        Production p;
        p.lhs = this;
        p.rhs.lhs = &other;
        p.rhs.op = Operations::VALUE;
        log->push_back(p);
    };
    Variable(Variable &&other) = delete;
    // Variable(const VariableInitializer &other)
    //     : log{other.log}
    // {}
    ~Variable() {}

private:
    Variable create_without_logging() const { return Variable{this->log}; }
    constexpr const Variable *get_id() const { return this; }
};

VariableExpression operator+(float lhs, const Variable &rhs)
{
    return {lhs, &rhs, Operations::ADD, rhs.log};
}

VariableExpression operator-(float lhs, const Variable &rhs)
{
    return {lhs, &rhs, Operations::SUBSTRACT, rhs.log};
}
VariableExpression operator*(float lhs, const Variable &rhs)
{
    return {lhs, &rhs, Operations::MULTIPLY, rhs.log};
}
VariableExpression operator/(float lhs, const Variable &rhs)
{
    return {lhs, &rhs, Operations::DIVIDE, rhs.log};
}

// Experimental: Expression-Variable operators allowing for ex. v0 = v1 + v2 * 2;
VariableExpression operator+(const Variable &lhs, const VariableExpression &rhs)
{
    Variable tmp = rhs;
    return lhs + tmp;
}

VariableExpression operator-(const Variable &lhs, const VariableExpression &rhs)
{
    Variable tmp = rhs;
    return lhs + tmp;
}
VariableExpression operator*(const Variable &lhs, const VariableExpression &rhs)
{
    Variable tmp = rhs;
    return lhs + tmp;
}
VariableExpression operator/(const Variable &lhs, const VariableExpression &rhs)
{
    Variable tmp = rhs;
    return lhs + tmp;
}

VariableExpression operator+(const VariableExpression &lhs, const Variable &rhs)
{
    Variable tmp = rhs;
    return tmp + rhs;
}

VariableExpression operator+(const VariableExpression &lhs, float rhs)
{
    Variable tmp = lhs;
    return tmp + rhs;
}

VariableExpression operator-(const VariableExpression &lhs, float rhs)
{
    Variable tmp = lhs;
    return tmp - rhs;
}

VariableExpression operator*(const VariableExpression &lhs, float rhs)
{
    Variable tmp = lhs;
    return tmp * rhs;
}

VariableExpression operator/(const VariableExpression &lhs, float rhs)
{
    Variable tmp = lhs;
    return tmp / rhs;
}

VariableExpression operator-(const VariableExpression &operand)
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
    auto format(trace::Operations op, FormatContext &ctx)
    {
        string_view name = "unknown";
        switch (op) {
        case trace::Operations::UNKNOWN:
            name = "UNKNOWN";
            break;
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
struct fmt::formatter<std::variant<std::monostate, const trace::Variable *, float>>
    : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const std::variant<std::monostate, const trace::Variable *, float> &var,
                FormatContext &ctx)
    {
        if (std::holds_alternative<const trace::Variable *>(var)) {
            return format_to(ctx.out(), "{}", fmt::ptr(std::get<const trace::Variable *>(var)));
        }
        if (std::holds_alternative<float>(var)) {
            return format_to(ctx.out(), "{}", std::get<float>(var));
        }
        throw std::runtime_error("Unknown variant");
    }
};

template<>
struct fmt::formatter<trace::Expression> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const trace::Expression &expr, FormatContext &ctx)
    {
        switch (expr.op) {
        case trace::Operations::UNKNOWN:
            throw std::runtime_error("UNKNOWN operation");
        case trace::Operations::ADD:
        case trace::Operations::SUBSTRACT:
        case trace::Operations::MULTIPLY:
        case trace::Operations::DIVIDE:
            return format_to(ctx.out(), "({} {} {})", expr.lhs, expr.op, expr.rhs);
        case trace::Operations::VALUE:
            return format_to(ctx.out(), "{}", expr.lhs);
            break;
        case trace::Operations::UNARY_SUBSTRACT:
            return format_to(ctx.out(), "(- {})", expr.lhs);
            break;
        }
        throw std::runtime_error("UNKNOWN operation");
    }
};

template<>
struct fmt::formatter<trace::Production> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const trace::Production &prod, FormatContext &ctx)
    {
        return format_to(ctx.out(), "{} := {}", fmt::ptr(prod.lhs), prod.rhs);
    }
};
