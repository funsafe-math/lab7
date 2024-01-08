// #include <blaze/Math.h>
#include "Expression.hpp"
#include "Matrix.hpp"
#include "trace.hpp"

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
    trace::Operations op;
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
        case analysis::VariableType::LITERAL:
            name = "literal";
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
    std::vector<trace::Production> log{};

    trace::Variable get_element() { return trace::Variable(&log); }

    void parse(std::span<trace::Variable> data)
    {
        //     std::vector<Production> productions{};
        //     std::map<const trace::Element *, size_t> tmp_element_to_ix{};
        //     size_t biggest_tmp_ix = 0;
        //     auto is_final = [&](const trace::Element *lhs) -> bool {
        //         return (lhs >= data.data() && ((lhs - data.data()) < data.size()));
        //     };
        //     auto generate_variable = [&](const trace::Element *lhs) -> Variable {
        //         if (is_final(lhs)) {
        //             Variable ret{};
        //             ret.type = VariableType::FINAL;
        //             ret.index = lhs - data.data();
        //             return ret;
        //         } else {
        //             Variable ret{};
        //             ret.type = VariableType::TEMPORARY;
        //             ret.index = 0; // TODO: fix
        //             if (tmp_element_to_ix.contains(lhs)) {
        //                 ret.index = tmp_element_to_ix.at(lhs);
        //             } else {
        //                 ret.index = biggest_tmp_ix;
        //                 tmp_element_to_ix[lhs] = biggest_tmp_ix++;
        //             }

        //             return ret;
        //         }
        //     };

        //     auto generate_rhs = [&](const trace::BinOp &op) -> decltype(Production::rhs) {
        //         if (op.rhs.op == trace::Operations::VALUE) {
        //             Variable var = generate_variable(op.rhs.lhs);
        //             return var;
        //         } else {
        //             BinOp bin_op;
        //             bin_op.lhs = generate_variable(op.rhs.lhs);
        //             bin_op.op = op.rhs.op;
        //             bin_op.rhs = generate_variable(std::get<const trace::Element *>(op.rhs.rhs));
        //             return bin_op;
        //         }
        //     };

        //     productions.reserve(log.size());
        //     for (const auto &op : log) {
        //         if (op.rhs.op == trace::Operations::DESTROY) {
        //             if (!is_final(op.lhs) && tmp_element_to_ix.contains(op.lhs)) {
        //                 tmp_element_to_ix.at(op.lhs) = biggest_tmp_ix++;
        //             }
        //             continue;
        //         }
        //         Production production;
        //         production.lhs = generate_variable(op.lhs);
        //         production.rhs = generate_rhs(op);
        //         productions.push_back(production);
        //     }

        //     for (const auto &prod : productions) {
        //         fmt::println("{}", prod);
        //     }
    }
};

} // namespace analysis

void test_with_matrix()
{
    analysis::Problem problem{};
    // std::vector<Element> vec{100, &log};

    // vec[0] = vec[1] + vec[0];
    // vec[1] = vec[1] + 0.1;

    Matrix<trace::Variable> mat{6, 5, problem.get_element()};
    // problem.log.clear();
    mat.gaussian_elimination();

    for (auto &v : problem.log) {
        fmt::println("{}", v);
    }
    fmt::println("##########################");
    problem.parse(mat.values_);
}

void simple_test()
{
    analysis::Problem problem{};
    std::vector<trace::Variable> vec{10, problem.get_element()};
    vec[0] += vec[1] + vec[2];

    problem.parse(vec);
}

int main()
{
    test_with_matrix();
    // simple_test();
}
