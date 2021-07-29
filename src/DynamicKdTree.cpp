#include <iostream>
#include <pre-graphics/DynamicKdTree>

namespace pre {

template <size_t Dim>
typename DynamicKdTree<Dim>::Int DynamicKdTree<Dim>::private_allocate() {
    if (free_ == Nil) {
        Int nodes_size = nodes_.size();
        if (nodes_.size() == 0)
            nodes_.resize(32);
        else
            nodes_.resize(nodes_.size() * 2);
        for (Int node = nodes_size; node < Int(nodes_.size()); node++) {
            nodes_[node].next = node + 1;
            nodes_[node].height = -1;
        }
        nodes_.back().next = Nil;
        free_ = nodes_size;
    }
    Int node = free_;
    free_ = nodes_[node].next;
    nodes_[node] = Node();
    return node;
}

template <size_t Dim>
void DynamicKdTree<Dim>::private_insert(Int leaf) {
    if (root_ == Nil) {
        root_ = leaf;
        nodes_[root_].parent = Nil;
        return;
    }
    Int node = root_;
    while (1) {
        Node& node_ref = nodes_[node];
        Node& leaf_ref = nodes_[leaf];
        Int& child = leaf_ref.point < node_ref //
                             ? node_ref.child0
                             : node_ref.child1;
        if (child == Nil) {
            child = leaf;
            leaf_ref.parent = node;
            leaf_ref.axis = pre::abs(node_ref.point - leaf_ref.point).argmax();
            break;
        }
        node = child;
    }

    node = nodes_[leaf].parent;
    while (node != Nil) {
        Node& node_ref = nodes_[node];
        Node& child0_ref = nodes_[node_ref.child0];
        Node& child1_ref = nodes_[node_ref.child1];
        if (node_ref.child0 != Nil)
            node_ref.height = std::max(node_ref.height, 1 + child0_ref.height);
        if (node_ref.child1 != Nil)
            node_ref.height = std::max(node_ref.height, 1 + child1_ref.height);
        node = node_ref.parent;
    }

    if (max_imbalance() > 4)
        private_rebalance(root_);
}

template <size_t Dim>
void DynamicKdTree<Dim>::private_remove(Int node) {
    // TODO
}

template <size_t Dim>
void DynamicKdTree<Dim>::private_rebalance(Int node) {
    // TODO
}

template class DynamicKdTree<2>;

template class DynamicKdTree<3>;

} // namespace pre
