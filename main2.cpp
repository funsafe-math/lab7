// #include <blaze/Math.h>
#include "Expression.hpp"
#include "Matrix.hpp"

#include "fmt/core.h"
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

namespace analysis {

enum class VariableType { UNKNOWN, TEMPORARY, FINAL, LITERAL };

struct Variable
{
    size_t index{};
    VariableType type{};
};

struct BinOp
{
    Variable lhs;
    Variable rhs;
    Operations op;
};

using Expression = std::variant<Variable, BinOp>;

struct Production
{
    Variable lhs{};
    Expression rhs{};
};

} // namespace analysis

template<>
struct fmt::formatter<analysis::VariableType> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(analysis::VariableType type, FormatContext &ctx)
    {
        string_view name = "unknown";
        switch (type) {
        case analysis::VariableType::UNKNOWN:
            name = "UNKNOWN";
            break;
        case analysis::VariableType::TEMPORARY:
            name = "temporary";
            break;
        case analysis::VariableType::FINAL:
            name = "final";
            break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

template<>
struct fmt::formatter<analysis::Variable> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(const analysis::Variable &variable, FormatContext &ctx)
    {
        return format_to(ctx.out(), "{}_{}", variable.type, variable.index);
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
namespace analysis {

struct Problem
{
    std::vector<::BinOp> log{};

    Element get_element() { return Element(&log); }
    void parse(std::span<Element> data, const std::function<std::string(size_t)> &namer)
    {
        std::vector<Production> productions{};
        std::map<::Element, size_t> tmp_element_to_ix{};
        auto generate_variable = [&](const ::Element *lhs) -> Variable {
            if (lhs - data.data() < data.size()) {
                Variable ret{};
                ret.type = VariableType::FINAL;
                ret.index = lhs - data.data();
                return ret;
            } else {
                Variable ret{};
                ret.type = VariableType::TEMPORARY;
                ret.index = 0; // TODO: fix
                return ret;
            }
        };

        auto generate_rhs = [&](const ::BinOp &op) -> decltype(Production::rhs) {
            if (op.rhs.op == Operations::VALUE) {
                Variable var = generate_variable(op.rhs.lhs);
                return var;
            } else {
                BinOp bin_op;
                bin_op.lhs = generate_variable(op.rhs.lhs);
                bin_op.op = op.rhs.op;
                bin_op.rhs = generate_variable(std::get<const Element *>(op.rhs.rhs));
                return bin_op;
            }
        };

        for (const auto &op : log) {
            Production production;
            production.lhs = generate_variable(op.lhs);
            production.rhs = generate_rhs(op);
            productions.push_back(production);
        }

        for (const auto &prod : productions) {
            fmt::println("{}", prod);
        }
    }
};

} // namespace analysis

int main()
{
    analysis::Problem problem{};
    // std::vector<Element> vec{100, &log};

    // vec[0] = vec[1] + vec[0];
    // vec[1] = vec[1] + 0.1;

    Matrix<Element> mat{6, 5, problem.get_element()};
    mat.gaussian_elimination();

    for (auto &v : problem.log) {
        fmt::println("{}", v);
    }
    fmt::println("##########################");
    problem.parse(mat.values_, [&](size_t ix) {
        return fmt::format("{}_{}", ix / mat.width_, ix % mat.width_);
    });
}
