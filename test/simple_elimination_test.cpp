#include "../src/Matrix.hpp"
#include <fstream>
#include <iostream>

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fmt::print(stderr, "Usage: {} <path to matrix file> <patch to correct output>\n", argv[0]);
        return 1;
    }

    Matrix<float> input = from_file(argv[1]);
    Matrix<float> correct_output = from_file(argv[2]);

    input.gaussian_elimination();

    fmt::println("Ours:");
    print_matrix(input);

    fmt::println("Correct:");
    print_matrix(correct_output);
    // print_matrix(input);

    fmt::println("Ours:    {}", input.values_);
    fmt::println("Correct: {}", correct_output.values_);

    std::cout << std::endl;
    for (size_t i = 0; i < input.values_.size(); ++i) {
        if (std::abs(input.values_.at(i) - correct_output.values_.at(i)) > 0.1) {
            throw std::runtime_error("Incorrect value");
        }
    }
    return 0;
}
