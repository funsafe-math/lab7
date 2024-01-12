#pragma once

#include "analysis.hpp"

namespace analysis {
struct Production;
struct Graph
{
    std::vector<Production> nodes{};
    std::vector<std::pair<size_t, size_t>> edges{};

    std::string as_dot();

    std::vector<size_t> topological_sort();

    bool hasPath(const size_t startNode, const size_t endNode);

    Graph transitive_reduction();
};

} // namespace analysis
