#include <algorithm>
#include <assert.h>
#include <bitset>
#include <chrono>
#include <fmt/chrono.h>
#include <functional>
#include <iomanip>
#include <map>
#include <random>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <stack>
#include <unordered_set>
#include <vector>

#include "Matrix.hpp"
#include "analysis.hpp"
#include "graph.hpp"
#include "problem.hpp"

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
    Matrix<float> ret{16, 15, 0};
    ret.values_ = {
        0.3598381494, 0.8409342858, 0.2078457000, 0.3076592230, 0.5746280031, 0.3193307736,
        0.3614672566, 0.7279397259, 0.9533114347, 0.0125489458, 0.3668514230, 0.0822051884,
        0.6728263458, 0.2414697087, 0.8907866449, 0.6384948161, 0.0990509280, 0.8448001488,
        0.7335959037, 0.4410333892, 0.0004641690, 0.2571903929, 0.2827449953, 0.1344851542,
        0.8379886466, 0.9408241372, 0.7638297249, 0.7754726642, 0.3612390876, 0.4475851775,
        0.6895253650, 0.8550308482, 0.0611367836, 0.5240200388, 0.3451626587, 0.9028476472,
        0.9895839527, 0.4061565322, 0.5826884804, 0.1881432569, 0.8365774747, 0.2832440536,
        0.0622547836, 0.5527200972, 0.3878068150, 0.8118781254, 0.9632170592, 0.4348680819,
        0.6823218116, 0.4663162978, 0.8145247654, 0.1093526633, 0.6058778448, 0.1784634465,
        0.9460446959, 0.9433581349, 0.9184359435, 0.7168164496, 0.6029857273, 0.4788331646,
        0.0059115238, 0.0046892389, 0.7971097761, 0.0930752870, 0.9368129044, 0.5912644209,
        0.9445518145, 0.3562214748, 0.6420745397, 0.5789400685, 0.5413422112, 0.7382183675,
        0.1538579404, 0.6128888873, 0.7117459899, 0.9483826822, 0.5363043470, 0.9832042774,
        0.1617978656, 0.7562935026, 0.5352875036, 0.1358311358, 0.5173959241, 0.3176994038,
        0.4072910492, 0.2258247967, 0.8745417488, 0.9285354458, 0.2594445028, 0.4034480317,
        0.2240388131, 0.6246113122, 0.7297402803, 0.6790216676, 0.6353992729, 0.1305411126,
        0.4114219336, 0.0139796539, 0.7853549691, 0.6315148523, 0.8068024749, 0.5601392847,
        0.3525633491, 0.2650392302, 0.0876001591, 0.6359691166, 0.2433253886, 0.7660147408,
        0.9022450545, 0.0176539539, 0.5241959248, 0.3408784191, 0.2785960288, 0.3268870842,
        0.7827245496, 0.2605054719, 0.4557001634, 0.2953058685, 0.2070036651, 0.5655408288,
        0.7126389005, 0.2165884607, 0.0632256706, 0.1091152343, 0.7542658297, 0.5096093700,
        0.1980688909, 0.3158685483, 0.5963010280, 0.5848972130, 0.7296366203, 0.1246918985,
        0.5640114466, 0.5862422241, 0.9399262573, 0.5094661313, 0.2112823140, 0.4058047461,
        0.9899567962, 0.8544495118, 0.8785187571, 0.3366953742, 0.4019945394, 0.9779708263,
        0.3812462110, 0.2216356655, 0.5895190368, 0.6389691158, 0.1072504053, 0.0889349987,
        0.3579798574, 0.6291786687, 0.3718725101, 0.6017271235, 0.8620720106, 0.5571317762,
        0.6181986809, 0.9641170863, 0.6162011510, 0.9447345766, 0.0948736058, 0.0646352080,
        0.0830413107, 0.0143082636, 0.4444321599, 0.0114416364, 0.1665537668, 0.9500438661,
        0.7895594488, 0.9186577436, 0.0865250001, 0.5972672054, 0.0137319643, 0.0647808592,
        0.4413314562, 0.8344664526, 0.6646218739, 0.8314134551, 0.5496519350, 0.9347766200,
        0.4076652435, 0.8008724220, 0.6403601923, 0.8993992383, 0.4330066054, 0.0393853043,
        0.9755776891, 0.5354885381, 0.8611343274, 0.8454461847, 0.9415865965, 0.0308326939,
        0.7115397145, 0.7280899996, 0.8476419914, 0.6123640034, 0.9775477666, 0.6008731913,
        0.7461419270, 0.0064400154, 0.7794512955, 0.2574379281, 0.1902886965, 0.0548470032,
        0.7312017753, 0.9393672943, 0.8445678922, 0.6042533948, 0.5973289158, 0.2003489442,
        0.1237908044, 0.2352716256, 0.9280816961, 0.3936953738, 0.9808339366, 0.5414606953,
        0.8976471080, 0.5914889976, 0.6838673279, 0.0963694169, 0.3441834018, 0.1717851197,
        0.5220487224, 0.8313023673, 0.9991857941, 0.6751889463, 0.5910245339, 0.3281276648,
        0.9021790219, 0.7498552792, 0.5425374148, 0.0792096039, 0.7344354286, 0.8225782968,
        0.1015977843, 0.6541523947, 0.8510025307, 0.1997580345, 0.7838164701, 0.6430521759,
    };

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
    for (const auto trace : problem.productions) {
        fmt::println("{}", trace);
    }

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
    problem.interpret_multithreaded(generated_solution.values_);

    for (size_t i = 0; i < correct_solution.values_.size(); i++) {
        fmt::println("Correct: {:5}, \tGenerated: {:5}",
                     correct_solution.values_[i],
                     generated_solution.values_[i]);
        if (correct_solution.values_[i] != generated_solution.values_[i]) {
            std::cout << std::endl;
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
    auto tick = std::chrono::high_resolution_clock::now();
    problem.interpret_multithreaded(generated_solution.values_);
    auto tock = std::chrono::high_resolution_clock::now();
    fmt::println("Multithreaded interpretation took {}",
                 std::chrono::duration_cast<std::chrono::microseconds>(tock - tick));

    {
        auto generated_solution = generate_bigger_matrix();
        auto tick = std::chrono::high_resolution_clock::now();
        problem.interpret(generated_solution.values_);
        auto tock = std::chrono::high_resolution_clock::now();
        fmt::println("Single-threaded interpretation took {}",
                     std::chrono::duration_cast<std::chrono::microseconds>(tock - tick));
    }

    for (size_t i = 0; i < correct_solution.values_.size(); i++) {
        fmt::println("Correct: {:5}, \tGenerated: {:5}",
                     correct_solution.values_[i],
                     generated_solution.values_[i]);
        if (correct_solution.values_[i] != generated_solution.values_[i]) {
            std::cout << std::endl;
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

    constexpr size_t width = 41;
    constexpr size_t height = 40;

    Matrix<trace::Variable> mat{width, height, problem.get_element()};
    problem.log.clear();
    mat.gaussian_elimination();

    // for (auto &v : problem.log) {
    //     fmt::println("{}", v);
    // }
    fmt::println("##########################");
    problem.parse(mat.values_);

    {
        Matrix<float> mat{width, height, 0};
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(1.0, 2.0);
        for (auto &v : mat.values_) {
            v = dis(gen);
        }

        auto tick = std::chrono::high_resolution_clock::now();
        problem.interpret_multithreaded(mat.values_);
        auto tock = std::chrono::high_resolution_clock::now();
        fmt::println("Multithreaded interpretation took {}",
                     std::chrono::duration_cast<std::chrono::microseconds>(tock - tick));
    }

    {
        Matrix<float> mat{width, height, 0};
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(1.0, 2.0);
        for (auto &v : mat.values_) {
            v = dis(gen);
        }

        auto tick = std::chrono::high_resolution_clock::now();
        problem.interpret(mat.values_);
        auto tock = std::chrono::high_resolution_clock::now();
        fmt::println("Singlethreaded interpretation took {}",
                     std::chrono::duration_cast<std::chrono::microseconds>(tock - tick));
        print_matrix(mat);
    }

    // fmt::println("Graph: ");
    // auto graph = problem.produce_dependence_graph();
    // std::string dot = graph.transitive_reduction().as_dot();
    // fmt::println("{}", graph.transitive_reduction().as_dot());
    fmt::println("Productions processed: {}", problem.log.size());

    auto elem = std::max_element(problem.FNF.begin(), problem.FNF.end(), [](auto a, auto b) {
        return a.size() < b.size();
    });
    fmt::println("For a {}x{} matrix, you can use at most {} threads",
                 mat.width_,
                 mat.height_,
                 elem->size());
}

// Constepxr test
constexpr size_t count_instructions()
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
    return problem.productions.size();
}

constexpr size_t matrix_count_instructions()
{
    analysis::Problem problem{};
    // std::vector<trace::Variable> vec{10, problem.get_element()};
    Matrix<trace::Variable> mat{8, 7, problem.get_element()};
    problem.log.clear();
    mat.gaussian_elimination();

    problem.parse(mat.values_);
    return problem.productions.size();
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


constexpr bool test_map()
{
    analysis::myMap<int, int> map{};
    map[2] = 5;
    map[5] = 1;
    map[6] = 1;
    map[7] = 1;
    map[8] = 1;
    map[9] = 1;
    if (!map.contains(7))
        return false;
    map[3] = 2;
    if (!map.contains(3))
        return false;
    map[3] = 3;

    return map[3] == 3;
}

static_assert(test_map() == true);

void generate_c_code()
{
    analysis::Problem problem{};
    // std::vector<Element> vec{100, &log};

    // vec[0] = vec[1] + vec[0];
    // vec[1] = vec[1] + 0.1;

    Matrix<trace::Variable> mat{6, 5, problem.get_element()};
    problem.log.clear();
    trace::Variable v = -mat.values_[0];
    mat.gaussian_elimination();
    // mat.values_[mat.values_.size() - 1] = -v;

    for (auto &v : problem.log) {
        fmt::println("{}", v);
    }
    fmt::println("##########################");
    problem.parse(mat.values_);
    // problem.generate_c_code();
}

// Store production variables efficiently
struct SuperProduction
{
    std::set<analysis::TemporaryVariable> lhs_tmp;
    std::set<analysis::FinalVariable> lhs_final;
    std::set<analysis::TemporaryVariable> rhs_tmp;
    std::set<analysis::FinalVariable> rhs_final;

    bool is_depent_lhs(const analysis::WritableVariable v)
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>) {
                    return rhs_tmp.contains(x);
                }
                if constexpr (std::is_same_v<T, analysis::FinalVariable>) {
                    return rhs_final.contains(x);
                }
            },
            v);
    }

    bool is_depent_rhs(const analysis::Variable v)
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>) {
                    return lhs_tmp.contains(x);
                }
                if constexpr (std::is_same_v<T, analysis::FinalVariable>) {
                    return lhs_final.contains(x);
                }
                if constexpr (std::is_same_v<T, analysis::LiteralVariable>) {
                    return false; // ignore
                }
            },
            v);
    }

    bool is_dependent_rhs_expr(const analysis::Expression &e)
    {
        return std::visit(
            [&](auto x) -> bool {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::Variable>) {
                    return is_depent_rhs(x);
                }
                if constexpr (std::is_same_v<T, analysis::BinOp>) {
                    return is_dependent_rhs(x);
                }
                if constexpr (std::is_same_v<T, analysis::UnaryOp>) {
                    return is_dependent_rhs(x);
                }
            },
            e);
    }

    bool is_dependent_rhs(const analysis::BinOp &v)
    {
        return is_depent_rhs(v.lhs) || is_depent_rhs(v.rhs);
    }
    bool is_dependent_rhs(const analysis::UnaryOp &op) { return is_depent_rhs(op.var); }

    bool is_dependent(const analysis::Production &p)
    {
        return is_depent_lhs(p.lhs) || is_dependent_rhs_expr(p.rhs);
    }

    constexpr void add_lhs(const analysis::WritableVariable v)
    {
        std::visit(
            [&](auto x) -> void {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>) {
                    lhs_tmp.insert(x);
                    return;
                }
                if constexpr (std::is_same_v<T, analysis::FinalVariable>) {
                    lhs_final.insert(x);
                    return;
                }
            },
            v);
    }

    constexpr void add_rhs_variable(const analysis::Variable v)
    {
        std::visit(
            [&](auto x) -> void {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::TemporaryVariable>) {
                    rhs_tmp.insert(x);
                    return;
                }
                if constexpr (std::is_same_v<T, analysis::FinalVariable>) {
                    rhs_final.insert(x);
                }
                if constexpr (std::is_same_v<T, analysis::LiteralVariable>) {
                    return; // ignore
                }
            },
            v);
    }

    void add_rhs_bin_op(const analysis::BinOp &v)
    {
        add_rhs_variable(v.lhs);
        add_rhs_variable(v.rhs);
    }
    void add_rhs_unary_op(const analysis::UnaryOp &op) { add_rhs_variable(op.var); }

    void add_rhs(const analysis::Expression &e)
    {
        std::visit(
            [&](auto x) -> void {
                using T = std::decay_t<decltype(x)>;
                if constexpr (std::is_same_v<T, analysis::Variable>) {
                    return add_rhs_variable(x);
                }
                if constexpr (std::is_same_v<T, analysis::BinOp>) {
                    return add_rhs_bin_op(x);
                }
                if constexpr (std::is_same_v<T, analysis::UnaryOp>) {
                    return add_rhs_unary_op(x);
                }
            },
            e);
    }

    void add(const analysis::Production &p)
    {
        add_lhs(p.lhs);
        add_rhs(p.rhs);
    }
};

std::vector<std::vector<size_t>> better_fnf(const std::span<analysis::Production> productions)
{
    std::vector<std::vector<size_t>> groups{};
    std::vector<size_t> group{};
    group.reserve(productions.size());
    SuperProduction group_sp{};
    SuperProduction not_added{};

    std::vector<size_t> remaining_productions{};
    for (size_t i = 0; i < productions.size(); ++i) {
        remaining_productions.push_back(i);
    }
    // If preceding production was not added, this means that someone in the group is dependent on it
    auto can_be_added = [&](const size_t p_ix) {
        // Can be optimized by creating a super-production from current group
        const analysis::Production &p = productions[p_ix];

        if (group_sp.is_dependent(p)) {
            return false;
        }
        if (not_added.is_dependent(p)) {
            return false;
        }
        return true;
    };

    auto step = [&]() {
        assert(remaining_productions.size());
        const size_t a_ix = remaining_productions.front();
        const analysis::Production &a = productions[a_ix];
        group.push_back(a_ix);
        group_sp.add(a);
        for (size_t b_ix : remaining_productions | std::ranges::views::drop(1)) {
            const analysis::Production &b = productions[b_ix];
            if (can_be_added(b_ix)) {
                group.push_back(b_ix);
                group_sp.add(b);
            } else {
                not_added.add(b);
            }
        }
        assert(group.size());
        // fmt::println("Group contents:");
        // for (const auto i : group) {
        //     fmt::println("{} : {}", i, productions[i]);
        // }
        // std::cout << std::flush;
    };

    while (remaining_productions.size() > 0) {
        const size_t starting_size = remaining_productions.size();
        step();

        // remove indexes from current group
        for (size_t i : group) {
            auto it = std::find(remaining_productions.begin(), remaining_productions.end(), i);
            remaining_productions.erase(it);
        }
        // add group to groups
        groups.emplace_back(std::move(group));
        group.clear(); // to make sure
        group_sp = SuperProduction();
        not_added = SuperProduction(); // somehow seems faster than .clear()
        assert(remaining_productions.size() < starting_size);
    }

    return groups;
}

void test_batter_fnf()
{
    analysis::Problem problem{};
    Matrix<trace::Variable> training{4, 3, problem.get_element()};
    problem.log.clear();
    training.gaussian_elimination();
    problem.parse(training.values_);
    // for (const auto trace : problem.productions) {
    //     fmt::println("{}", trace);
    // }

    auto generated_solution = generate_matrix();
    problem.interpret(generated_solution.values_);

    std::vector<std::vector<size_t>> fnf = better_fnf(problem.productions);
    // std::vector<std::vector<size_t>> fnf = sorting_fnf(problem.productions);
    // Compare solutions
    for (size_t i = 0; i < problem.FNF.size(); ++i) {
        if (fnf.size() <= i) {
            fmt::println("test fnf does not have group at index {}", i);
            return;
        }
        std::vector<size_t> &good = problem.FNF.at(i);
        std::vector<size_t> &test = problem.FNF.at(i);
        if (good.size() != test.size()) {
            fmt::println("test fnf.at({}) has size {} but correct has size {}",
                         i,
                         test.size(),
                         good.size());
        }

        std::sort(good.begin(), good.end());
        std::sort(test.begin(), test.end());
        for (int j = 0; j < good.size(); ++j) {
            assert(good.at(j) == test.at(j));
        }
    }

    fmt::println("both FNFs are the same");
}

void test_with_matrix_better_fnf()
{
    analysis::Problem problem{};
    // std::vector<Element> vec{100, &log};

    // vec[0] = vec[1] + vec[0];
    // vec[1] = vec[1] + 0.1;

    constexpr size_t width = 41;
    constexpr size_t height = 40;

    Matrix<trace::Variable> mat{width, height, problem.get_element()};
    problem.log.clear();
    mat.gaussian_elimination();

    // for (auto &v : problem.log) {
    //     fmt::println("{}", v);
    // }
    fmt::println("##########################");
    // problem.parse(mat.values_);
    problem.productions = analysis::Problem::parse_pure(mat.values_, problem.log);
    problem.FNF = better_fnf(problem.productions);

    {
        Matrix<float> mat{width, height, 0};
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(1.0, 2.0);
        for (auto &v : mat.values_) {
            v = dis(gen);
        }

        auto tick = std::chrono::high_resolution_clock::now();
        problem.interpret_multithreaded(mat.values_);
        auto tock = std::chrono::high_resolution_clock::now();
        fmt::println("Multithreaded interpretation took {}",
                     std::chrono::duration_cast<std::chrono::microseconds>(tock - tick));
    }

    {
        Matrix<float> mat{width, height, 0};
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(1.0, 2.0);
        for (auto &v : mat.values_) {
            v = dis(gen);
        }

        auto tick = std::chrono::high_resolution_clock::now();
        problem.interpret(mat.values_);
        auto tock = std::chrono::high_resolution_clock::now();
        fmt::println("Singlethreaded interpretation took {}",
                     std::chrono::duration_cast<std::chrono::microseconds>(tock - tick));
        print_matrix(mat);
    }

    // fmt::println("Graph: ");
    // auto graph = problem.produce_dependence_graph();
    // std::string dot = graph.transitive_reduction().as_dot();
    // fmt::println("{}", graph.transitive_reduction().as_dot());
    fmt::println("Productions processed: {}", problem.log.size());

    auto elem = std::max_element(problem.FNF.begin(), problem.FNF.end(), [](auto a, auto b) {
        return a.size() < b.size();
    });
    fmt::println("For a {}x{} matrix, you can use at most {} threads",
                 mat.width_,
                 mat.height_,
                 elem->size());
}

int main()
{
    // test_batter_fnf();
    // test_implementation();
    // test_with_matrix();
    test_with_matrix_better_fnf();
    // generate_c_code();
    // multithreaded_test_implementation();
    // big_multithreaded_test_implementation();

    // std::bitset<matrix_count_instructions()> test{};
    // fmt::println("Used instructions: {}", test.size());

    // return std::
    // std::bitset<count_instructions()> test{};
    // simple_test();
}
