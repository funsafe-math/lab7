#pragma once

#include <cstddef>
#include <span>

#include "trace.hpp"

namespace analysis {

struct TemporaryVariable
{
    size_t number;
    TemporaryVariable(size_t number_)
        : number{number_}
    {}
    float interpret(std::span<float> temps) const { return temps[number]; }
};

struct FinalVariable
{
    std::ptrdiff_t index;
    FinalVariable(std::ptrdiff_t index_)
        : index{index_}
    {}

    float interpret(std::span<float> finals) const { return finals[index]; }
};

struct LiteralVariable
{
    float value;
    LiteralVariable(float value_)
        : value{value_}
    {}
    float interpret() const { return value; }
};

struct Variable : public std::variant<TemporaryVariable, FinalVariable, LiteralVariable>
{
    using std::variant<TemporaryVariable, FinalVariable, LiteralVariable>::variant;
    float interpret(std::span<float> finals, std::span<float> temps) const
    {
        return std::visit(
            [&](auto x) -> float {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, TemporaryVariable>)
                    return x.interpret(temps);
                if constexpr (std::is_same_v<T, FinalVariable>)
                    return x.interpret(finals);
                if constexpr (std::is_same_v<T, LiteralVariable>)
                    return x.interpret();
            },
            *this);
    }
};

// You cant have literalVariable on production lhs
using WritableVariable = std::variant<TemporaryVariable, FinalVariable>;

struct BinOp
{
    Variable lhs;
    Variable rhs;
    trace::Operations op;
    float interpret(std::span<float> finals, std::span<float> temps) const
    {
        auto interpret = [&](Variable var) { return var.interpret(finals, temps); };
        switch (op) {
        case trace::Operations::ADD:
            return interpret(lhs) + interpret(rhs);
        case trace::Operations::SUBSTRACT:
            return interpret(lhs) - interpret(rhs);
        case trace::Operations::MULTIPLY:
            return interpret(lhs) * interpret(rhs);
        case trace::Operations::DIVIDE:
            return interpret(lhs) / interpret(rhs);
        case trace::Operations::VALUE:
            return interpret(lhs);
        case trace::Operations::UNARY_SUBSTRACT: // Should be an unaryOp, but whatever
            return -interpret(lhs);
        }
    }
};

struct UnaryOp
{
    Variable var;
    trace::Operations op; // must be a unary operation
    float interpret(std::span<float> finals, std::span<float> temps) const
    {
        auto interpret = [&](Variable var) { return var.interpret(finals, temps); };
        switch (op) {
        case trace::Operations::ADD:
            return std::numeric_limits<float>::quiet_NaN();
        case trace::Operations::SUBSTRACT:
            return std::numeric_limits<float>::quiet_NaN();
        case trace::Operations::MULTIPLY:
            return std::numeric_limits<float>::quiet_NaN();
        case trace::Operations::DIVIDE:
            return std::numeric_limits<float>::quiet_NaN();
        case trace::Operations::VALUE:
            return interpret(var);
        case trace::Operations::UNARY_SUBSTRACT:
            return -interpret(var);
        }
    }
};

struct Expression : public std::variant<Variable, BinOp, UnaryOp>
{
    using std::variant<Variable, BinOp, UnaryOp>::variant;

    float interpret(std::span<float> finals, std::span<float> temps) const
    {
        return std::visit(
            [&](auto x) -> float {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, Variable>)
                    return x.interpret(finals, temps);
                if constexpr (std::is_same_v<T, BinOp>)
                    return x.interpret(finals, temps);
                if constexpr (std::is_same_v<T, UnaryOp>)
                    return x.interpret(finals, temps);
            },
            *this);
    }
};

struct Production
{
    WritableVariable lhs;
    Expression rhs;
    Production(const WritableVariable &lhs_, const Expression &rhs_)
        : lhs{lhs_}
        , rhs{rhs_}
    {}
};

} // namespace analysis

template<>
struct fmt::formatter<analysis::Variable> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(analysis::Variable variable, FormatContext &ctx)
    {
        // TODO: finish
        return std::visit(
            [&](const auto x) -> auto {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, std::monostate>) {
                    throw std::runtime_error("This should never happen");
                    return format_to(ctx.out(), "Error");
                }
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>)
                    return format_to(ctx.out(), "temporary_{}", x.number);
                if constexpr (std::is_same_v<T, analysis::FinalVariable>)
                    return format_to(ctx.out(), "final_{}", x.index);
                if constexpr (std::is_same_v<T, analysis::LiteralVariable>)
                    return format_to(ctx.out(), "(literal: {})", x.value);
            },
            variable);
    }
};

template<>
struct fmt::formatter<analysis::WritableVariable> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(analysis::WritableVariable variable, FormatContext &ctx)
    {
        // TODO: finish
        return std::visit(
            [&](const auto x) -> auto {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>)
                    return format_to(ctx.out(), "temporary_{}", x.number);
                if constexpr (std::is_same_v<T, analysis::FinalVariable>)
                    return format_to(ctx.out(), "final_{}", x.index);
            },
            variable);
        // return formatter<string_view>::format(name, ctx);
    }
};

template<>
struct fmt::formatter<analysis::BinOp> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const analysis::BinOp &binary_operation, FormatContext &ctx)
    {
        return format_to(ctx.out(),
                         "({} {} {})",
                         binary_operation.lhs,
                         binary_operation.op,
                         binary_operation.rhs);
    }
};

template<>
struct fmt::formatter<analysis::Expression> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const analysis::Expression &expr, FormatContext &ctx)
    {
        if (std::holds_alternative<analysis::Variable>(expr)) {
            return format_to(ctx.out(), "{}", std::get<analysis::Variable>(expr));
        }
        if (std::holds_alternative<analysis::BinOp>(expr)) {
            return format_to(ctx.out(), "{}", std::get<analysis::BinOp>(expr));
        }
        throw std::runtime_error("This should never happen");
    }
};

template<>
struct fmt::formatter<analysis::Production> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const analysis::Production &prod, FormatContext &ctx)
    {
        return format_to(ctx.out(), "{} := {}", prod.lhs, prod.rhs);
    }
};
