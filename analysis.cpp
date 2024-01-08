#include <assert.h>
#include <map>
#include <span>
#include <vector>

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
    std::vector<Production> productions{};

    trace::Variable get_element() { return trace::Variable(&log); }

    void parse(std::span<trace::Variable> data)
    {
        std::map<const trace::Variable *, size_t> tmp_variable_to_ix{};
        size_t next_tmp_ix = 0;

        auto is_final = [&](const trace::Variable *lhs) -> bool {
            return (lhs >= data.data() && ((lhs - data.data()) < data.size()));
        };

        auto generate_writable_variable_from_ptr =
            [&](const trace::Variable *lhs) -> WritableVariable {
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
            [&](const std::variant<const trace::Variable *, float> var) -> Variable {
            return std::visit(
                [&](const auto &val) -> Variable {
                    using T = std::decay_t<decltype(val)>;
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
            throw std::runtime_error("Invalid variant");
        };

        productions.reserve(log.size());
        for (const trace::Production &prod : log) {
            // First process rhs, then lhs - tmp variable indexes are changed on write
            Expression rhs = generate_rhs(prod.rhs);
            // if (prod.rhs.)
            // TODO: deal with tmp variable incrementation
            WritableVariable lhs = generate_writable_variable_from_ptr(prod.lhs);

            if (std::holds_alternative<TemporaryVariable>(lhs)) {
                auto &tmp = std::get<TemporaryVariable>(lhs);
                tmp.number = next_tmp_ix++;
                tmp_variable_to_ix[prod.lhs] = tmp.number;
            }
            productions.emplace_back(lhs, rhs);
        }

        //

        for (const Production &prod : productions) {
            fmt::println("{}", prod);
        }

        produce_relations();
        print_relations();
        produce_FNF();
        printFNF();
    }

    using RelationT = std::vector<std::pair<size_t, size_t>>;

    RelationT D{};
    RelationT I{};

    void produce_relations()
    {
        for (size_t a_ix = 0; a_ix < productions.size(); ++a_ix) {
            const Production &a = productions.at(a_ix);
            for (size_t b_ix = a_ix; b_ix < productions.size(); ++b_ix) {
                const Production &b = productions.at(b_ix);
                if (a.is_dependent(b)) {
                    D.emplace_back(a_ix, b_ix);
                } else {
                    I.emplace_back(a_ix, b_ix);
                }
            }
        }
    }

    void print_relations() const
    {
        auto print_arr = [&](const RelationT &set) {
            for (const auto [a_ix, b_ix] : set) {
                if (productions.at(a_ix) == productions.at(b_ix)) {
                    continue;
                }
                fmt::println("{{ {},\t{} }}", productions.at(a_ix), productions.at(b_ix));
            }
        };
        fmt::println("Dependent (identity and swapped skipped):");
        print_arr(D);
        fmt::println("Independent:");
        print_arr(I);
    }

    std::vector<std::vector<size_t>> FNF{};

    void produce_FNF()
    {
        constexpr size_t empty_token = std::numeric_limits<size_t>::max();
        std::vector<std::vector<size_t>> stacks{productions.size(), std::vector<size_t>()};
        // populate stack

        // TODO: use vector<bool> (will use only a bit instead of 8 bytes per entry)
        for (size_t a_ix = productions.size() - 1; a_ix < productions.size(); a_ix--) {
            stacks.at(a_ix).push_back(a_ix);
            const Production &a = productions.at(a_ix);
            for (size_t b_ix = 0; b_ix < productions.size(); ++b_ix) {
                const Production &b = productions.at(b_ix);
                if (a_ix == b_ix)
                    continue;
                if (a.is_dependent(b)) {
                    stacks.at(b_ix).push_back(empty_token);
                }
            }
        }

        // empty the stack
        std::vector<size_t> group{};
        while (true) {
            group.clear();
            std::vector<size_t> empty_stacks{};
            for (size_t i = 0; i < stacks.size(); ++i) {
                auto &stack = stacks.at(i);
                if (stack.empty()) {
                    empty_stacks.push_back(i);
                    continue;
                }

                if (stack.back() == empty_token) {
                    continue;
                }

                auto to_be_pushed_back = stack.back();
                group.push_back(to_be_pushed_back);
                stack.pop_back();
            }

            for (const auto &name : empty_stacks) {
                // stacks.erase(name);
            }

            for (const size_t a_name : group) {
                const Production &a = productions.at(a_name);
                for (size_t b_name = 0; b_name < stacks.size(); ++b_name) {
                    auto &stack = stacks.at(b_name);
                    if (a_name == b_name) {
                        continue;
                    }
                    if (stack.empty())
                        continue;
                    const Production &b = productions.at(b_name);
                    if (a.is_dependent(b)
                        // && stack.back() == empty_token
                    ) {
                        assert(stack.back() == empty_token);
                        stack.pop_back();
                    }
                }
            }
            if (group.empty())
                return;
            // fmt::println("Group: {}", group);
            FNF.emplace_back(std::move(group));
            group = std::vector<size_t>();
            if (empty_stacks.size() == stacks.size()) {
                assert(false); // We should never get here
                return;
            }
        }
    }

    void printFNF() const
    {
        for (const std::vector<size_t> &group : FNF) {
            fmt::println("### Group seperator ###");
            for (size_t ix : group) {
                fmt::println("{}", productions.at(ix));
            }

            fmt::println("### Group seperator ###");
        }
    }

    static void interpret_production(const Production &prod,
                                     std::span<float> &data,
                                     std::vector<float> &temporary_variables)
    {
        auto parse_lhs = [&](WritableVariable var) -> float & {
            return std::visit(
                [&](const auto x) -> float & {
                    using T = std::decay_t<decltype(x)>;
                    if constexpr (std::is_same_v<T, TemporaryVariable>)
                        return temporary_variables[x.number];
                    if constexpr (std::is_same_v<T, FinalVariable>)
                        return data[x.index];
                },
                var);
        };

        float &lhs = parse_lhs(prod.lhs);
        lhs = prod.rhs.interpret(data, temporary_variables);
    }

    // Only for testing purposes
    void interpret(std::span<float> data) const
    {
        // Interpret the pseudo-assembly
        size_t n_temporary_variables = 0;
        for (const Production &prod : productions) {
            if (std::holds_alternative<TemporaryVariable>(prod.lhs)) {
                auto &tmp = std::get<TemporaryVariable>(prod.lhs);
                n_temporary_variables = std::max(n_temporary_variables, tmp.number);
            }
        }

        std::vector<float> temporary_variables(n_temporary_variables + 1, float{});

        for (const Production &prod : productions) {
            interpret_production(prod, data, temporary_variables);
        }
    }
};

} // namespace analysis

constexpr Matrix<float> generate_matrix()
{
    Matrix<float> ret{4, 3, 0};
    ret.values_ = {2, 1, 3, 6, 4, 3, 8, 15, 6, 5, 16, 2};
    return ret;
}

void test_implementation()
{
    auto correct_solution = generate_matrix();
    correct_solution.gaussian_elimination();

    analysis::Problem problem{};
    Matrix<trace::Variable> training{4, 3, problem.get_element()};
    problem.log.clear();
    training.gaussian_elimination();
    problem.parse(training.values_);

    auto generated_solution = generate_matrix();
    problem.interpret(generated_solution.values_);

    for (size_t i = 0; i < correct_solution.values_.size(); i++) {
        fmt::println("Correct: {:5}, \tGenerated: {:5}",
                     correct_solution.values_[i],
                     generated_solution.values_[i]);
        if (correct_solution.values_[i] != generated_solution.values_[i]) {
            throw std::runtime_error("Incorrect");
        }
    }
    fmt::println("Correct:");
    print_matrix(correct_solution);
    fmt::println("Generated:");
    print_matrix(generated_solution);
}

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
    problem.log.clear();
    for (auto &v : vec) {
        v = v + 1.0;
    }
    vec[0] *= 2.0;
    vec[1] += 1.0;

    problem.parse(vec);
}

int main()
{
    // test_implementation();
    // test_with_matrix();
    simple_test();
}
