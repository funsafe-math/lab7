#include <map>
#include <span>

#include "Matrix.hpp"
#include "analysis.hpp"

template<class... Ts>
struct overloaded : Ts...
{
    using Ts::operator()...;
};
// explicit deduction guide (not needed as of C++20)
template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

namespace analysis {

struct Problem
{
    std::vector<trace::Production> log{};

    trace::Variable get_element() { return trace::Variable(&log); }

    void parse(std::span<trace::Variable> data)
    {
        std::vector<Production> productions{};
        std::map<const trace::Variable *, size_t> tmp_variable_to_ix{};
        size_t next_tmp_ix = 0;

        auto is_final = [&](const trace::Variable *lhs) -> bool {
            return (lhs >= data.data() && ((lhs - data.data()) < data.size()));
        };

        auto generate_variable_from_ptr = [&](const trace::Variable *lhs) -> Variable {
            if (is_final(lhs)) {
                return FinalVariable{lhs - data.data()};
            } else {
                size_t index = 0;
                // ret.index = 0; // TODO: fix
                if (tmp_variable_to_ix.contains(lhs)) {
                    auto index = tmp_variable_to_ix.at(lhs);
                    return TemporaryVariable{index};
                } else {
                    auto index = next_tmp_ix++;
                    tmp_variable_to_ix[lhs] = index;
                    return TemporaryVariable{index};
                }
            }
        };

        auto generate_variable_from_variant =
            [&](const std::variant<std::monostate, const trace::Variable *, float> &var) -> Variable {
            return std::visit(
                [&](const auto &val) -> Variable {
                    using T = std::decay_t<decltype(val)>;
                    if constexpr (std::is_same_v<T, std::monostate>)
                        return Variable{};
                    if constexpr (std::is_same_v<T, const trace::Variable *>)
                        return Variable{generate_variable_from_ptr(val)};
                    if constexpr (std::is_same_v<T, float>)
                        return Variable{LiteralVariable{val}};
                },
                var);
        };

        auto generate_rhs = [&](const trace::Expression &expr) -> Expression {
            switch (expr.op) {
            case trace::Operations::ADD:
            case trace::Operations::SUBSTRACT:
            case trace::Operations::MULTIPLY:
            case trace::Operations::DIVIDE:
                return BinOp{generate_variable_from_variant(expr.lhs),
                             generate_variable_from_variant(expr.rhs),
                             expr.op};
            case trace::Operations::VALUE:
                return generate_variable_from_variant(expr.lhs);
            case trace::Operations::UNARY_SUBSTRACT:
                return UnaryOp{generate_variable_from_variant(expr.lhs), expr.op};
            }
        };

        productions.reserve(log.size());
        for (const trace::Production &prod : log) {
            // First process rhs, then lhs - tmp variable indexes are changed on write
            Production production;
            production.rhs = generate_rhs(prod.rhs);
            // if (prod.rhs.)
            // TODO: deal with tmp variable incrementation
            production.lhs = generate_variable_from_ptr(prod.lhs);

            if (std::holds_alternative<TemporaryVariable>(production.lhs)) {
                auto &tmp = std::get<TemporaryVariable>(production.lhs);
                tmp.number = next_tmp_ix++;
                tmp_variable_to_ix[prod.lhs] = tmp.number;
            }
            productions.push_back(production);
        }

        //

        for (const Production &prod : productions) {
            fmt::println("{}", prod);
        }
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
