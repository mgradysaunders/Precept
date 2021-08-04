#pragma once

#include <vector>

template <typename Node, typename NodeAlloc, typename Int>
inline size_t GrowPoolVector(
    std::vector<Node, NodeAlloc>& pool, Int Node::*next) {
    auto size = pool.size();
    pool.resize(size == 0 ? 32 : 2 * size);
    for (auto node = size; node < pool.size(); node++)
        pool[node].*next = node + 1;
    pool.back().*next = -1;
    return size;
}
