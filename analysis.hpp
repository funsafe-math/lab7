#pragma once

#include <cstddef>

#include "trace.hpp"

namespace analysis {

struct TemporaryVariable
{
    size_t number;
    TemporaryVariable(size_t number_)
        : number{number_}
    {}
};

struct FinalVariable
{
    std::ptrdiff_t index;
    FinalVariable(std::ptrdiff_t index_)
        : index{index_}
    {}
};

struct LiteralVariable
{
    float value;
    LiteralVariable(float value_)
        : value{value_}
    {}
};

using Variable = std::variant<std::monostate, TemporaryVariable, FinalVariable, LiteralVariable>;

struct BinOp
{
    Variable lhs;
    Variable rhs;
    trace::Operations op;
};

struct UnaryOp
{
    Variable var;
    trace::Operations op; // must be a unary operation
};

using Expression = std::variant<Variable, BinOp, UnaryOp>;

struct Production
{
    Variable lhs{};
    Expression rhs{};
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
                    return format_to(ctx.out(), "(literal: {}", x.value);
            },
            variable);
        // return formatter<string_view>::format(name, ctx);
    }
};

// template<>
// struct fmt::formatter<analysis::Variable> : formatter<string_view>
// {
//     template<typename FormatContext>
//     auto format(const analysis::Variable &variable, FormatContext &ctx)
//     {
//         return format_to(ctx.out(), "{}_{}", variable.type, variable.index);
//     }
// };

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
