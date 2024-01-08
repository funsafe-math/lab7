#include "graph.hpp"

#include <functional>
#include <stack>
#include <unordered_set>

std::string analysis::Graph::as_dot()
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

std::vector<size_t> analysis::Graph::topological_sort()
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

bool analysis::Graph::hasPath(const size_t startNode, const size_t endNode)
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

analysis::Graph analysis::Graph::transitive_reduction()
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
