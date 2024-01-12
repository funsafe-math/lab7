#pragma once
#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <cstddef>
#include <span>

#include "trace.hpp"

namespace analysis {

struct TemporaryVariable
{
    size_t number;
    constexpr TemporaryVariable(size_t number_)
        : number{number_}
    {}
    constexpr float interpret(std::span<float> temps) const { return temps[number]; }
    constexpr auto operator<=>(const TemporaryVariable &) const = default;
};

struct FinalVariable
{
    std::ptrdiff_t index;
    constexpr FinalVariable(std::ptrdiff_t index_)
        : index{index_}
    {}

    constexpr float interpret(std::span<float> finals) const { return finals[index]; }
    constexpr auto operator<=>(const FinalVariable &) const = default;
};

struct LiteralVariable
{
    float value;
    constexpr LiteralVariable(float value_)
        : value{value_}
    {}
    constexpr float interpret() const { return value; }
    constexpr auto operator<=>(const LiteralVariable &) const = default;
};

// You cant have literalVariable on production lhs
using WritableVariable = std::variant<TemporaryVariable, FinalVariable>;

struct Variable : public std::variant<TemporaryVariable, FinalVariable, LiteralVariable>
{
    using std::variant<TemporaryVariable, FinalVariable, LiteralVariable>::variant;
    constexpr auto operator<=>(const Variable &) const = default;
    constexpr float interpret(std::span<float> finals, std::span<float> temps) const
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

    constexpr bool operator==(const TemporaryVariable w) const
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, TemporaryVariable>)
                    return x == w;
                if constexpr (std::is_same_v<T, FinalVariable>)
                    return false;
                if constexpr (std::is_same_v<T, LiteralVariable>)
                    return false;
            },
            *this);
    }
    constexpr bool operator==(const FinalVariable w) const
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, TemporaryVariable>)
                    return false;
                if constexpr (std::is_same_v<T, FinalVariable>)
                    return x == w;
                if constexpr (std::is_same_v<T, LiteralVariable>)
                    return false;
            },
            *this);
    }

    constexpr bool operator==(const WritableVariable w) const
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, TemporaryVariable>)
                    return *this == x;
                if constexpr (std::is_same_v<T, FinalVariable>)
                    return *this == x;
            },
            w);
    }
};

struct BinOp
{
    Variable lhs;
    Variable rhs;
    trace::Operations op;
    constexpr float interpret(std::span<float> finals, std::span<float> temps) const
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

    constexpr bool contains(const WritableVariable w) const { return lhs == w || rhs == w; }

    constexpr auto operator<=>(const BinOp &) const = default;
};

struct UnaryOp
{
    Variable var;
    trace::Operations op; // must be a unary operation
    constexpr float interpret(std::span<float> finals, std::span<float> temps) const noexcept
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

    constexpr bool contains(const WritableVariable w) { return var == w; }

    constexpr auto operator<=>(const UnaryOp &) const = default;
};

struct Expression : public std::variant<Variable, BinOp, UnaryOp>
{
    using std::variant<Variable, BinOp, UnaryOp>::variant;

    float interpret(std::span<float> finals, std::span<float> temps) const noexcept
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
    constexpr bool contains(const WritableVariable w) const noexcept
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, Variable>)
                    return x == w;
                if constexpr (std::is_same_v<T, BinOp>)
                    return x.contains(w);
                if constexpr (std::is_same_v<T, UnaryOp>)
                    return x.contains(w);
            },
            *this);
    }

    auto operator<=>(const Expression &) const = default;
};

struct Production
{
    WritableVariable lhs;
    Expression rhs;
    constexpr Production(const WritableVariable lhs_, const Expression &rhs_) noexcept
        : lhs{lhs_}
        , rhs{rhs_}
    {}

    constexpr bool is_dependent(const Production &other) const noexcept
    {
        // return lhs == other.lhs || this->rhs.contains(other.lhs) || other.rhs.contains(this->lhs);
        return this->rhs.contains(other.lhs) || other.rhs.contains(this->lhs);
    }
    constexpr auto operator<=>(const Production &) const noexcept = default;
};

} // namespace analysis

template<>
struct fmt::formatter<analysis::Variable> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(analysis::Variable variable, FormatContext &ctx)
    {
        // TODO: finish
        return std::visit(
            [&](const auto x) -> auto {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, std::monostate>) {
                    throw std::runtime_error("This should never happen");
                    return fmt::format_to(ctx.out(), "Error");
                }
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>)
                    return fmt::format_to(ctx.out(), "temporary_{}", x.number);
                if constexpr (std::is_same_v<T, analysis::FinalVariable>)
                    return fmt::format_to(ctx.out(), "final_{}", x.index);
                if constexpr (std::is_same_v<T, analysis::LiteralVariable>)
                    return fmt::format_to(ctx.out(), "(literal: {})", x.value);
            },
            variable);
    }
};

template<>
struct fmt::formatter<analysis::WritableVariable> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(analysis::WritableVariable variable, FormatContext &ctx)
    {
        // TODO: finish
        return std::visit(
            [&](const auto x) -> auto {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>)
                    return fmt::format_to(ctx.out(), "temporary_{}", x.number);
                if constexpr (std::is_same_v<T, analysis::FinalVariable>)
                    return fmt::format_to(ctx.out(), "final_{}", x.index);
            },
            variable);
        // return formatter<string_view>::format(name, ctx);
    }
};

template<>
struct fmt::formatter<analysis::BinOp> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const analysis::BinOp &binary_operation, FormatContext &ctx)
    {
        return fmt::format_to(ctx.out(),
                              "({} {} {})",
                              binary_operation.lhs,
                              binary_operation.op,
                              binary_operation.rhs);
    }
};

template<>
struct fmt::formatter<analysis::UnaryOp> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const analysis::UnaryOp &unary_operation, FormatContext &ctx)
    {
        return fmt::format_to(ctx.out(), "({}{})", unary_operation.op, unary_operation.var);
    }
};

template<>
struct fmt::formatter<analysis::Expression> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const analysis::Expression &expr, FormatContext &ctx)
    {
        if (std::holds_alternative<analysis::Variable>(expr)) {
            return fmt::format_to(ctx.out(), "{}", std::get<analysis::Variable>(expr));
        }
        if (std::holds_alternative<analysis::BinOp>(expr)) {
            return fmt::format_to(ctx.out(), "{}", std::get<analysis::BinOp>(expr));
        }

        if (std::holds_alternative<analysis::UnaryOp>(expr)) {
            return fmt::format_to(ctx.out(), "{}", std::get<analysis::UnaryOp>(expr));
        }
        throw std::runtime_error("This should never happen");
    }
};

template<>
struct fmt::formatter<analysis::Production> : formatter<string_view>
{
    template<typename FormatContext>
    constexpr auto format(const analysis::Production &prod, FormatContext &ctx)
    {
        return fmt::format_to(ctx.out(), "{} := {}", prod.lhs, prod.rhs);
    }
};
#endif // ANALYSIS_HPP
