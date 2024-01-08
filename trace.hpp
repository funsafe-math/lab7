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
    inline static Expression makeBinOp(const Variable *lhs, Operations op, const Variable *rhs)
    {
        Expression ret;
        ret.lhs = lhs;
        ret.op = op;
        ret.rhs = rhs;
        return ret;
    }

    // inline static Expression makeUnaryOp(const Variable *val, Operations op)
    // {
    //     Expression ret;
    //     ret.lhs = val;
    //     ret.op = op;
    //     return ret;
    // }
    inline static Expression makeLiteral(const float val)
    {
        Expression ret;
        ret.lhs = val;
        ret.rhs = val;
        ret.op = Operations::VALUE;
        return ret;
    }

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

    operator Expression() const { return Expression{op, lhs, rhs}; }

    LiteralExpression(float literal)
        : op{Operations::VALUE}
        , lhs{literal}
        , rhs{literal}
    {}
};

struct VariableExpression : public Expression
{
    std::vector<Production> *log;

    inline static VariableExpression makeUnaryOp(const Variable *val,
                                                 Operations op,
                                                 std::vector<Production> *log)
    {
        VariableExpression ret;
        ret.lhs = val;
        ret.op = op;
        ret.log = log;
        return ret;
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
        *this = static_cast<Expression>(expr);
    }

    // Arithmetic operators with self

    VariableExpression operator+(const Variable &other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::ADD;
        expr.rhs = other.get_id();
        expr.log = log;
        return expr;
    }

    VariableExpression operator*(const Variable &other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::MULTIPLY;
        expr.rhs = other.get_id();
        expr.log = log;
        return expr;
    }

    VariableExpression operator/(const Variable &other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::DIVIDE;
        expr.rhs = other.get_id();
        return expr;
    }

    VariableExpression operator-(const Variable &other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::SUBSTRACT;
        expr.rhs = other.get_id();
        return expr;
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
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::ADD;
        expr.rhs = other;
        expr.log = log;
        return expr;
    }

    VariableExpression operator*(const float other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::MULTIPLY;
        expr.rhs = other;
        expr.log = log;
        return expr;
    }

    VariableExpression operator/(const float other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::DIVIDE;
        expr.rhs = other;
        expr.log = log;
        return expr;
    }

    VariableExpression operator-(const float other) const
    {
        VariableExpression expr{};
        expr.lhs = get_id();
        expr.op = Operations::SUBSTRACT;
        expr.rhs = other;
        expr.log = log;
        return expr;
    }

    // Assignment operators

    Variable &operator=(const Expression &expr)
    {
        Production prod;
        prod.lhs = get_id();
        prod.rhs = expr;
        log->push_back(prod);
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

    Variable &operator+=(const Expression &expr)
    {
        Variable tmp = create_without_logging();
        tmp = expr;
        *this = *this + tmp;
        return *this;
    }

    Variable &operator-=(const Expression &expr)
    {
        Variable tmp = create_without_logging();
        tmp = expr;
        *this = *this - tmp;
        return *this;
    }

    Variable &operator*=(const Expression &expr)
    {
        Variable tmp = create_without_logging();
        tmp = expr;
        *this = *this * tmp;
        return *this;
    }

    Variable &operator/=(const Expression &expr)
    {
        Variable tmp = create_without_logging();
        tmp = expr;
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

    friend Expression operator+(const Variable &lhs, const Expression &rhs);
    friend Expression operator-(const Variable &lhs, const Expression &rhs);
    friend Expression operator*(const Variable &lhs, const Expression &rhs);
    friend Expression operator/(const Variable &lhs, const Expression &rhs);

    friend Expression operator+(const Expression &lhs, const Expression &rhs);
    friend Expression operator-(const Expression &lhs, const Expression &rhs);
    friend Expression operator*(const Expression &lhs, const Expression &rhs);
    friend Expression operator/(const Expression &lhs, const Expression &rhs);

    friend Expression operator+(const Expression &lhs, const Variable &rhs);

private:
    Variable create_without_logging() const { return Variable{this->log}; }
    constexpr const Variable *get_id() const { return this; }
};

VariableExpression operator+(float lhs, const Variable &rhs)
{
    VariableExpression expr{};
    expr.lhs = lhs;
    expr.op = Operations::ADD;
    expr.rhs = &rhs;
    expr.log = rhs.log;
    return expr;
}

VariableExpression operator-(float lhs, const Variable &rhs)
{
    VariableExpression expr{};
    expr.lhs = lhs;
    expr.op = Operations::SUBSTRACT;
    expr.rhs = &rhs;
    expr.log = rhs.log;
    return expr;
}
VariableExpression operator*(float lhs, const Variable &rhs)
{
    VariableExpression expr{};
    expr.lhs = lhs;
    expr.op = Operations::MULTIPLY;
    expr.rhs = &rhs;
    expr.log = rhs.log;
    return expr;
}
VariableExpression operator/(float lhs, const Variable &rhs)
{
    VariableExpression expr{};
    expr.lhs = lhs;
    expr.op = Operations::DIVIDE;
    expr.rhs = &rhs;
    expr.log = rhs.log;
    return expr;
}

// Experimental: Expression-Variable operators allowing for ex. v0 = v1 + v2 * 2;
Expression operator+(const Variable &lhs, const Expression &rhs)
{
    Variable tmp = lhs.create_without_logging();
    tmp = rhs;
    return lhs + tmp;
}

Expression operator-(const Variable &lhs, const Expression &rhs)
{
    Variable tmp = lhs.create_without_logging();
    tmp = rhs;
    return lhs + tmp;
}
Expression operator*(const Variable &lhs, const Expression &rhs)
{
    Variable tmp = lhs.create_without_logging();
    tmp = rhs;
    return lhs + tmp;
}
Expression operator/(const Variable &lhs, const Expression &rhs)
{
    Variable tmp = lhs.create_without_logging();
    tmp = rhs;
    return lhs + tmp;
}

Expression operator+(const Expression &lhs, const Variable &rhs)
{
    Variable tmp = rhs.create_without_logging();
    tmp = lhs;
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
