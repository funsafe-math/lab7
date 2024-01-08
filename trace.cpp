#include "trace.hpp"

void test1()
{
    std::vector<trace::Production> log;
    std::vector<trace::Variable> vec{10, trace::Variable{&log}};
    auto tmp1 = -(vec[2] + 1);
    auto tmp2 = 5 * tmp1 * 2;
    auto tmp3 = tmp2 - 4;
    vec[0] = 1.12;

    for (const auto &p : log) {
        fmt::println("{}", p);
    }
}

int main(int argc, char *argv[])
{
    test1();
    return 0;
}
