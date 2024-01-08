#include <assert.h>
#include <functional>
#include <iomanip>
#include <map>
#include <span>
#include <sstream>
#include <stack>
#include <unordered_set>
#include <vector>

#include "Matrix.hpp"
#include "analysis.hpp"

namespace analysis {

struct Graph
{
    std::vector<Production> nodes{};
    std::vector<std::pair<size_t, size_t>> edges{};

    std::string as_dot()
    {
        std::string ret;
        auto inserter = std::back_inserter(ret);
        fmt::format_to(inserter, "digraph g{{\n");
        for (const auto [from, to] : edges) {
            fmt::format_to(inserter, "  {} -> {};\n", from, to);
        }
        for (size_t i = 0; i < nodes.size(); ++i) {
            fmt::format_to(inserter, "  {}[label=\"{}\"];\n", i, nodes.at(i));
        }
        fmt::format_to(inserter, "}}");
        return ret;
    }

    std::vector<size_t> topological_sort()
    {
        std::vector<size_t> sortedNodes;
        std::vector<bool> visited(nodes.size(), false);

        std::function<void(int)> dfs = [&](int nodeIndex) {
            visited[nodeIndex] = true;

            for (const auto &edge : edges) {
                if (edge.first == nodeIndex && !visited[edge.second]) {
                    dfs(edge.second);
                }
            }

            sortedNodes.push_back(nodeIndex);
        };

        for (int i = 0; i < nodes.size(); ++i) {
            if (!visited[i]) {
                dfs(i);
            }
        }

        std::reverse(sortedNodes.begin(), sortedNodes.end());

        return sortedNodes;
    }

    bool hasPath(const size_t startNode, const size_t endNode)
    {
        std::unordered_set<size_t> visited;
        std::stack<size_t> stack;
        stack.push(startNode);

        while (!stack.empty()) {
            size_t currentNode = stack.top();
            stack.pop();

            if (currentNode == endNode) {
                return true;
            }

            visited.insert(currentNode);

            for (const auto &edge : edges) {
                if (edge.first == currentNode && !visited.contains(edge.second)) {
                    stack.push(edge.second);
                }
            }
        }

        return false;
    }

    Graph transitive_reduction()
    {
        Graph ret{};
        ret.nodes = this->nodes;
        std::vector<size_t> t = topological_sort();

        for (int i = 0; i < t.size(); ++i) {
            for (int j = t.size() - 1; j >= 0; --j) {
                const std::pair<size_t, size_t> edge{t[j], t[i]};
                if (std::ranges::find(edges, edge) != edges.end()) {
                    if (!ret.hasPath(edge.first, edge.second)) {
                        ret.edges.push_back(edge);
                    }
                }
            }
        }

        return ret;
    }
};

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

        // for (const Production &prod : productions) {
        //     fmt::println("{}", prod);
        // }

        // produce_relations();
        // print_relations();
        produce_FNF();
        // printFNF();

        Graph graph = produce_dependence_graph();
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
        for (const std::vector<size_t> group : FNF) {
            for (const size_t ix : group) {
                const Production &prod = productions.at(ix);
                interpret_production(prod, data, temporary_variables);
            }
        }
    }

    Graph produce_dependence_graph()
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
};

} // namespace analysis

std::string url_encode(const std::string &value)
{
    std::ostringstream escaped;
    escaped.fill('0');
    escaped << std::hex;

    for (char c : value) {
        // Keep alphanumeric and other accepted characters intact
        if (isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~') {
            escaped << c;
            continue;
        }

        // Any other characters are percent-encoded
        escaped << std::uppercase;
        escaped << '%' << std::setw(2) << int((unsigned char) c);
        escaped << std::nouppercase;
    }

    return escaped.str();
}

constexpr Matrix<float> generate_matrix()
{
    Matrix<float> ret{4, 3, 0};
    ret.values_ = {2, 1, 3, 6, 4, 3, 8, 15, 6, 5, 16, 2};
    return ret;
}

constexpr Matrix<float> generate_bigger_matrix()
{
    Matrix<float> ret{15, 16, 0};
    ret.values_
        = {0.3598381494, 0.8409342858, 0.2078457000, 0.3076592230, 0.5746280031, 0.3193307736,
           0.3614672566, 0.7279397259, 0.9533114347, 0.0125489458, 0.3668514230, 0.0822051884,
           0.6728263458, 0.2414697087, 0.8907866449, 0.0990509280, 0.8448001488, 0.7335959037,
           0.4410333892, 0.0004641690, 0.2571903929, 0.2827449953, 0.1344851542, 0.8379886466,
           0.9408241372, 0.7638297249, 0.7754726642, 0.3612390876, 0.4475851775, 0.6895253650,
           0.0611367836, 0.5240200388, 0.3451626587, 0.9028476472, 0.9895839527, 0.4061565322,
           0.5826884804, 0.1881432569, 0.8365774747, 0.2832440536, 0.0622547836, 0.5527200972,
           0.3878068150, 0.8118781254, 0.9632170592, 0.6823218116, 0.4663162978, 0.8145247654,
           0.1093526633, 0.6058778448, 0.1784634465, 0.9460446959, 0.9433581349, 0.9184359435,
           0.7168164496, 0.6029857273, 0.4788331646, 0.0059115238, 0.0046892389, 0.7971097761,
           0.9368129044, 0.5912644209, 0.9445518145, 0.3562214748, 0.6420745397, 0.5789400685,
           0.5413422112, 0.7382183675, 0.1538579404, 0.6128888873, 0.7117459899, 0.9483826822,
           0.5363043470, 0.9832042774, 0.1617978656, 0.5352875036, 0.1358311358, 0.5173959241,
           0.3176994038, 0.4072910492, 0.2258247967, 0.8745417488, 0.9285354458, 0.2594445028,
           0.4034480317, 0.2240388131, 0.6246113122, 0.7297402803, 0.6790216676, 0.6353992729,
           0.4114219336, 0.0139796539, 0.7853549691, 0.6315148523, 0.8068024749, 0.5601392847,
           0.3525633491, 0.2650392302, 0.0876001591, 0.6359691166, 0.2433253886, 0.7660147408,
           0.9022450545, 0.0176539539, 0.5241959248, 0.2785960288, 0.3268870842, 0.7827245496,
           0.2605054719, 0.4557001634, 0.2953058685, 0.2070036651, 0.5655408288, 0.7126389005,
           0.2165884607, 0.0632256706, 0.1091152343, 0.7542658297, 0.5096093700, 0.1980688909,
           0.5963010280, 0.5848972130, 0.7296366203, 0.1246918985, 0.5640114466, 0.5862422241,
           0.9399262573, 0.5094661313, 0.2112823140, 0.4058047461, 0.9899567962, 0.8544495118,
           0.8785187571, 0.3366953742, 0.4019945394, 0.3812462110, 0.2216356655, 0.5895190368,
           0.6389691158, 0.1072504053, 0.0889349987, 0.3579798574, 0.6291786687, 0.3718725101,
           0.6017271235, 0.8620720106, 0.5571317762, 0.6181986809, 0.9641170863, 0.6162011510,
           0.0948736058, 0.0646352080, 0.0830413107, 0.0143082636, 0.4444321599, 0.0114416364,
           0.1665537668, 0.9500438661, 0.7895594488, 0.9186577436, 0.0865250001, 0.5972672054,
           0.0137319643, 0.0647808592, 0.4413314562, 0.6646218739, 0.8314134551, 0.5496519350,
           0.9347766200, 0.4076652435, 0.8008724220, 0.6403601923, 0.8993992383, 0.4330066054,
           0.0393853043, 0.9755776891, 0.5354885381, 0.8611343274, 0.8454461847, 0.9415865965,
           0.7115397145, 0.7280899996, 0.8476419914, 0.6123640034, 0.9775477666, 0.6008731913,
           0.7461419270, 0.0064400154, 0.7794512955, 0.2574379281, 0.1902886965, 0.0548470032,
           0.7312017753, 0.9393672943, 0.8445678922, 0.5973289158, 0.2003489442, 0.1237908044,
           0.2352716256, 0.9280816961, 0.3936953738, 0.9808339366, 0.5414606953, 0.8976471080,
           0.5914889976, 0.6838673279, 0.0963694169, 0.3441834018, 0.1717851197, 0.5220487224,
           0.9991857941, 0.6751889463, 0.5910245339, 0.3281276648, 0.9021790219, 0.7498552792,
           0.5425374148, 0.0792096039, 0.7344354286, 0.8225782968, 0.1015977843, 0.6541523947,
           0.8510025307, 0.1997580345, 0.7838164701, 0.6384948161, 0.8550308482, 0.4348680819,
           0.0930752870, 0.7562935026, 0.1305411126, 0.3408784191, 0.3158685483, 0.9779708263,
           0.9447345766, 0.8344664526, 0.0308326939, 0.6042533948, 0.8313023673, 0.6430521759};

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

void multithreaded_test_implementation()
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

void big_multithreaded_test_implementation()
{
    auto correct_solution = generate_bigger_matrix();
    fmt::println("Starting matrix:");
    print_matrix(correct_solution);
    correct_solution.gaussian_elimination();

    analysis::Problem problem{};
    Matrix<trace::Variable> training{correct_solution.width_,
                                     correct_solution.height_,
                                     problem.get_element()};
    problem.log.clear();
    training.gaussian_elimination();
    problem.parse(training.values_);

    auto generated_solution = generate_bigger_matrix();
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

    Matrix<trace::Variable> mat{16, 15, problem.get_element()};
    problem.log.clear();
    mat.gaussian_elimination();

    for (auto &v : problem.log) {
        fmt::println("{}", v);
    }
    fmt::println("##########################");
    problem.parse(mat.values_);

    // fmt::println("Graph: ");
    // auto graph = problem.produce_dependence_graph();
    // std::string dot = graph.transitive_reduction().as_dot();
    // fmt::println("{}", graph.transitive_reduction().as_dot());
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
    // multithreaded_test_implementation();
    big_multithreaded_test_implementation();
    // simple_test();
}
