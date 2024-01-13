#pragma once

#include "analysis.hpp"
#include "graph.hpp"
#include "trace.hpp"
#include <algorithm>
#include <cassert>
#include <execution>
#include <stdexcept>
#include <thread>
#include <vector>

namespace analysis {

template<typename K, typename V>
struct myMap
{
    std::vector<std::pair<K, V>> store{};

    constexpr bool contains(const K key)
    {
        // TODO: replace with std::binary_search
        std::pair<K, V> p{key, V{}};
        return std::binary_search(store.begin(),
                                  store.end(),
                                  p,
                                  [](std::pair<K, V> a, std::pair<K, V> b) {
                                      return a.first < b.first;
                                  });
    }

    constexpr V &at(const K key)
    {
        auto it = std::find_if(store.begin(), store.end(), [key](std::pair<K, V> pair) {
            return pair.first == key;
        });
        if (it == store.end()) {
            throw std::invalid_argument("Invalid arg provided");
        }
        return it->second;
    }

    constexpr V &operator[](const K key)
    {
        auto it = std::find_if(store.begin(), store.end(), [key](std::pair<K, V> pair) {
            return pair.first == key;
        });
        if (it == store.end()) {
            insert(key, V{});
        }
        return at(key);
    }

    constexpr void insert(const K key, const V value)
    {
        // store.push_back({key, value});
        // std::sort(store.begin(), store.end());
        // key
        std::pair p{key, value};

        auto cmp = [](std::pair<K, V> a, std::pair<K, V> b) { return a.first < b.first; };
        store.insert(std::upper_bound(store.begin(), store.end(), p, cmp), p);
    }
};

struct Problem
{
    std::vector<trace::Production> log{};
    std::vector<Production> productions{};

    constexpr trace::Variable get_element() { return trace::Variable(&log); }

    /**
     * @brief Generate a vector of productions from log
     * @param data - range of variables to be marked final, rather than temporary
     * @param log - raw trace::Productions list to be processed
     * @return vector of productions
     */
    static constexpr std::vector<Production> parse_pure(const std::span<trace::Variable> data,
                                                        const std::span<trace::Production> log)
    {
        // std::map<const trace::Variable *, size_t> tmp_variable_to_ix{};
        myMap<size_t, size_t> tmp_variable_to_ix{};

        // HACK: allow compiler to compare pointers from different allocations at compile time
        std::vector<const trace::Variable *> pointers{};

        auto ptr_to_ix = [&](const trace::Variable *ptr) -> size_t {
            auto it = std::find(pointers.begin(), pointers.end(), ptr);
            if (it == pointers.end()) {
                pointers.push_back(ptr);
                return pointers.size() - 1;
            }
            return std::distance(pointers.begin(), it);
        };

        size_t next_tmp_ix = 0;

        auto is_final = [&](const trace::Variable *lhs) constexpr -> bool {
            // auto diff = lhs - data.data();
            // return (diff >= 0 && diff < data.size());
            // return (&(*data.begin()) <= lhs && lhs < &(*data.end()));
            for (size_t i = 0; i < data.size(); i++) {
                const auto s = std::span(lhs, 1);
                auto subspan = data.subspan(i, 1);
                if (subspan.data() == s.data()) {
                    return true;
                }
            }
            return false;
        };

        auto generate_writable_variable_from_ptr =
            [&](const trace::Variable *lhs) constexpr -> WritableVariable {
            if (is_final(lhs)) {
                return FinalVariable{lhs - data.data()};
            } else {
                size_t index = 0;
                if (tmp_variable_to_ix.contains(ptr_to_ix(lhs))) {
                    auto index = tmp_variable_to_ix.at(ptr_to_ix(lhs));
                    return TemporaryVariable{index};
                } else {
                    auto index = next_tmp_ix++;
                    tmp_variable_to_ix[ptr_to_ix(lhs)] = index;
                    return TemporaryVariable{index};
                }
            }
        };

        auto generate_variable_from_ptr = [&](const trace::Variable *lhs) constexpr -> Variable {
            if (is_final(lhs)) {
                return FinalVariable{lhs - data.data()};
            } else {
                size_t index = 0;
                if (tmp_variable_to_ix.contains(ptr_to_ix(lhs))) {
                    auto index = tmp_variable_to_ix.at(ptr_to_ix(lhs));
                    return TemporaryVariable{index};
                } else {
                    auto index = next_tmp_ix++;
                    tmp_variable_to_ix[ptr_to_ix(lhs)] = index;
                    return TemporaryVariable{index};
                }
            }
        };

        auto generate_variable_from_variant =
            [&](const std::variant<const trace::Variable *, float> var) constexpr -> Variable {
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

        std::vector<Production> productions{};
        productions.reserve(log.size());
        for (const trace::Production &prod : log) {
            // First process rhs, then lhs - tmp variable indexes are changed on write
            Expression rhs = generate_rhs(prod.rhs);
            WritableVariable lhs = generate_writable_variable_from_ptr(prod.lhs);

            // Deal with tmp variable incrementation
            if (std::holds_alternative<TemporaryVariable>(lhs)) {
                auto &tmp = std::get<TemporaryVariable>(lhs);
                tmp.number = next_tmp_ix++;
                tmp_variable_to_ix[ptr_to_ix(prod.lhs)] = tmp.number;
            }
            productions.push_back({lhs, rhs});
        }

        return productions;
    }

    constexpr void parse(std::span<trace::Variable> data)
    {
        productions = parse_pure(data, log);
        // produce_relations();
        // print_relations();
        // produce_FNF_multithreaded();
        produce_FNF();
        // printFNF();

        // Graph graph = produce_dependence_graph();
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
                    // D.emplace_back(b_ix, a_ix);
                } else {
                    I.emplace_back(a_ix, b_ix);
                    // I.emplace_back(b_ix, a_ix);
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

    void produce_FNF_multithreaded() noexcept
    {
        constexpr bool empty_token = false;
        constexpr bool non_empty_token = true;
        std::vector<std::vector<bool>> stacks{productions.size(), std::vector<bool>()};
        // populate stack

        uint n_threads = std::thread::hardware_concurrency();
        std::vector<decltype(stacks)> ministacks(n_threads, stacks);

        std::vector<std::thread> threads{};

        for (size_t i = 0; i < n_threads; ++i) {
            auto func = [&](size_t thread_ix) {
                size_t chunk_size = (productions.size() + n_threads) / n_threads;
                size_t starting_ix = chunk_size * thread_ix;
                size_t ending_ix = chunk_size * (thread_ix + 1);
                fmt::println("{} Processing {} items, from {} to {}",
                             thread_ix,
                             ending_ix - starting_ix,
                             starting_ix,
                             ending_ix);

                ending_ix = std::min(ending_ix, productions.size());
                std::vector<std::vector<bool>> &stacks = ministacks.at(thread_ix);
                for (auto &stack : stacks)
                    stack.reserve(productions.size());

                for (size_t _a_ix = starting_ix; _a_ix < ending_ix; _a_ix++) {
                    size_t a_ix = productions.size() - 1 - _a_ix;
                    stacks.at(a_ix).push_back(non_empty_token);
                    const Production a = productions.at(a_ix);
                    // TODO: probably it would be okay to just go from a_ix,max or to 0,a_ix
                    for (size_t b_ix = 0; b_ix < productions.size(); ++b_ix) {
                        // const Production b = productions.at(b_ix);
                        const Production b = productions[b_ix];
                        if (a_ix == b_ix)
                            continue;
                        if (a.is_dependent(b)) {
                            stacks.at(b_ix).push_back(empty_token);
                        }
                    }
                }
            };
            threads.emplace_back(std::thread{func, i});
        }

        for (auto &th : threads) {
            th.join();
        }

        for (size_t i = 0; i < productions.size(); ++i) {
            std::vector<bool> &master = stacks.at(i);
            master.reserve(productions.size());
            for (const auto &st : ministacks) {
                const std::vector<bool> &child = st.at(i);
                master.insert(master.end(), child.begin(), child.end());
            }
        }

        // empty the stack
        std::vector<size_t> group{};
        while (true) {
            group.clear();
            // std::vector<size_t> empty_stacks{};
            for (size_t i = 0; i < stacks.size(); ++i) {
                auto &stack = stacks.at(i);
                if (stack.empty()) {
                    // empty_stacks.push_back(i);
                    continue;
                }

                if (stack.back() == empty_token) {
                    continue;
                }

                group.push_back(i);
                stack.pop_back();
            }

            // for (const auto &name : empty_stacks) {
            //     // stacks.erase(name);
            // }

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
            // if (empty_stacks.size() == stacks.size()) {
            //     assert(false); // We should never get here
            //     return;
            // }
        }
    }

    constexpr void produce_FNF()
    {
        constexpr bool empty_token = false;
        constexpr bool non_empty_token = true;
        std::vector<std::vector<bool>> stacks{productions.size(), std::vector<bool>()};
        // populate stack

        for (size_t a_ix = productions.size() - 1; a_ix < productions.size(); a_ix--) {
            stacks.at(a_ix).push_back(non_empty_token);
            const Production a = productions.at(a_ix);
            for (size_t b_ix = 0; b_ix < productions.size(); ++b_ix) {
                // const Production b = productions.at(b_ix);
                const Production b = productions[b_ix];
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
            // std::vector<size_t> empty_stacks{};
            for (size_t i = 0; i < stacks.size(); ++i) {
                auto &stack = stacks.at(i);
                if (stack.empty()) {
                    // empty_stacks.push_back(i);
                    continue;
                }

                if (stack.back() == empty_token) {
                    continue;
                }

                group.push_back(i);
                stack.pop_back();
            }

            // for (const auto &name : empty_stacks) {
            //     // stacks.erase(name);
            // }

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
            // if (empty_stacks.size() == stacks.size()) {
            //     assert(false); // We should never get here
            //     return;
            // }
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

    // Only for testing purposes
    void interpret_multithreaded(std::span<float> data) const
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

        // for (const Production &prod : productions) {
        //     interpret_production(prod, data, temporary_variables);
        // }
        for (const std::vector<size_t> &group : FNF) {
            std::for_each(std::execution::par_unseq, group.begin(), group.end(), [&](size_t ix) {
                const Production &prod = productions.at(ix);
                interpret_production(prod, data, temporary_variables);
            });
            //     for (const size_t ix : group) {
            //     const Production &prod = productions.at(ix);
            //     interpret_production(prod, data, temporary_variables);
            // }
        }
    }

    constexpr Graph produce_dependence_graph()
    {
        Graph ret{};
        ret.nodes = productions;

        for (int a_ix = 0; a_ix < ret.nodes.size(); ++a_ix) {
            const Production &a = ret.nodes.at(a_ix);
            for (int b_ix = a_ix + 1; b_ix < ret.nodes.size(); ++b_ix) {
                const Production &b = ret.nodes.at(b_ix);
                if (a.is_dependent(b)) {
                    ret.edges.push_back({a_ix, b_ix});
                }
            }
        }

        return ret;
    }

    // void generate_c_code()
    // {
    //     // TODO: count how many temporary variables should exist at a time
    //     size_t n_temporary_variables = 0;
    //     for (const Production &prod : productions) {
    //         if (std::holds_alternative<TemporaryVariable>(prod.lhs)) {
    //             auto &tmp = std::get<TemporaryVariable>(prod.lhs);
    //             n_temporary_variables = std::max(n_temporary_variables, tmp.number);
    //         }
    //     }

    //     struct VariableLifetime
    //     {
    //         size_t start = std::numeric_limits<size_t>::max();
    //         size_t end = std::numeric_limits<size_t>::min();
    //         auto operator<=>(const VariableLifetime &) const = default;
    //     };

    //     std::vector<VariableLifetime> temporary_variables_life(n_temporary_variables + 1,
    //                                                            VariableLifetime{});

    //     size_t index{0};
    //     for (const Production &prod : productions) {
    //         auto add_to_vec = [&](TemporaryVariable tmp) {
    //             VariableLifetime &life = temporary_variables_life.at(tmp.number);
    //             life.start = std::min(index, life.start);
    //             life.end = std::max(index, life.start);
    //             index++;
    //         };
    //         if (std::holds_alternative<TemporaryVariable>(prod.lhs)) {
    //             auto &tmp = std::get<TemporaryVariable>(prod.lhs);
    //             add_to_vec(tmp);
    //         }

    //         auto match_variable = [&](Variable var) {
    //             std::visit(
    //                 [&](auto x) {
    //                     using T = std::decay_t<decltype(x)>;
    //                     if constexpr (std::is_same_v<T, TemporaryVariable>)
    //                         add_to_vec(x);
    //                 },
    //                 var);
    //         };

    //         std::visit(
    //             [&](auto x) {
    //                 using T = std::decay_t<decltype(x)>;
    //                 if constexpr (std::is_same_v<T, Variable>)
    //                     match_variable(x);
    //                 if constexpr (std::is_same_v<T, BinOp>) {
    //                     match_variable(x.lhs);
    //                     match_variable(x.rhs);
    //                 }
    //                 if constexpr (std::is_same_v<T, UnaryOp>)
    //                     match_variable(x.var);
    //             },
    //             prod.rhs);
    //     }

    //     std::vector<int> temporary_variables_per_step(index, 0);
    //     for (const VariableLifetime &life : temporary_variables_life) {
    //         if (life == VariableLifetime{})
    //             continue;
    //         temporary_variables_per_step.at(life.start)++;
    //         temporary_variables_per_step.at(life.end)--;
    //     }

    //     int64_t max_variables_at_a_time{0};
    //     int64_t variables_at_a_time{0};
    //     for (auto v : temporary_variables_per_step) {
    //         variables_at_a_time += v;
    //         max_variables_at_a_time = std::max(variables_at_a_time, max_variables_at_a_time);
    //     }
    //     // max_variables_at_a_time += 1;

    //     assert(max_variables_at_a_time >= 0);
    //     fmt::println("Max temporary variables used: {}", max_variables_at_a_time);

    //     // for (const Production &prod : productions) {
    //     //     interpret_production(prod, data, temporary_variables);
    //     // }
    // }
};

} // namespace analysis
