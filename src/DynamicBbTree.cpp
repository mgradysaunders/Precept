#include <iostream>
#include <pre-graphics/DynamicBbTree>

namespace pre {

template <size_t Dim>
typename DynamicBbTree<Dim>::Int DynamicBbTree<Dim>::private_allocate() {
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
void DynamicBbTree<Dim>::private_insert(Int leaf) {
    if (root_ == Nil) {
        root_ = leaf;
        nodes_[root_].parent = Nil;
        return;
    }
    Box leaf_box = nodes_[leaf].box;
    Int node = root_;
    while (nodes_[node].is_branch()) {
        Node& node_ref = nodes_[node];
        Node& child0_ref = nodes_[node_ref.child0];
        Node& child1_ref = nodes_[node_ref.child1];
        float area = node_ref.box.surface_area();
        float combined_area = (leaf_box | node_ref.box).surface_area();
        float cost = 2 * combined_area;
        float cost_child0 = 2 * (combined_area - area);
        float cost_child1 = 2 * (combined_area - area);
        cost_child0 += (leaf_box | child0_ref.box).surface_area();
        cost_child1 += (leaf_box | child1_ref.box).surface_area();
        if (child0_ref.is_branch())
            cost_child0 -= child0_ref.box.surface_area();
        if (child1_ref.is_branch())
            cost_child1 -= child1_ref.box.surface_area();
        if (cost < cost_child0 and //
            cost < cost_child1)
            break;
        else
            node = cost_child0 < cost_child1 ? node_ref.child0
                                             : node_ref.child1;
    }
    Int old_parent = nodes_[node].parent;
    Int new_parent = private_allocate();
    Node& new_parent_ref = nodes_[new_parent];
    new_parent_ref.box = leaf_box | nodes_[node].box;
    new_parent_ref.parent = old_parent;
    new_parent_ref.child0 = node;
    new_parent_ref.child1 = leaf;
    new_parent_ref.height = nodes_[node].height + 1;
    nodes_[node].parent = new_parent;
    nodes_[leaf].parent = new_parent;
    if (old_parent != Nil) {
        if (nodes_[old_parent].child0 == node)
            nodes_[old_parent].child0 = new_parent;
        else
            nodes_[old_parent].child1 = new_parent;
    }
    else
        root_ = new_parent; // Sibling was root

    // Traverse back to root, rebalancing and updating heights.
    node = nodes_[leaf].parent;
    while (node != Nil) {
        node = private_balance(node);
        Node& node_ref = nodes_[node];
        Node& child0_ref = nodes_[node_ref.child0];
        Node& child1_ref = nodes_[node_ref.child1];
        node_ref.height = std::max(child0_ref.height, child1_ref.height);
        node_ref.height++;
        node_ref.box = child0_ref.box | child1_ref.box;
        node = node_ref.parent;
    }
}

template <size_t Dim>
void DynamicBbTree<Dim>::private_remove(Int leaf) {
    if (root_ == leaf) {
        root_ = Nil;
        return;
    }
    Int parent = nodes_[leaf].parent;
    Int sibling = nodes_[parent].child0 == leaf ? nodes_[parent].child1
                                                : nodes_[parent].child0;
    if (Int grandparent = nodes_[parent].parent; grandparent != Nil) {
        if (nodes_[grandparent].child0 == parent)
            nodes_[grandparent].child0 = sibling;
        else
            nodes_[grandparent].child1 = sibling;
        nodes_[sibling].parent = grandparent;
        private_deallocate(parent);

        Int node = grandparent;
        while (node != Nil) {
            node = private_balance(node);
            Node& node_ref = nodes_[node];
            Node& child0_ref = nodes_[node_ref.child0];
            Node& child1_ref = nodes_[node_ref.child1];
            node_ref.box = child0_ref.box | child1_ref.box;
            node_ref.height = std::max(child0_ref.height, child1_ref.height);
            node_ref.height++;
            node = node_ref.parent;
        }
    }
    else {
        root_ = sibling;
        nodes_[sibling].parent = Nil;
        private_deallocate(parent);
    }
}

template <size_t Dim>
typename DynamicBbTree<Dim>::Int DynamicBbTree<Dim>::private_balance(
        Int node) {
    // Rotate B above its parent A.
    auto rotate = [&](Int b, Int a) {
        Node& a_ref = nodes_[a];
        Node& b_ref = nodes_[b];
        Node& c_ref = nodes_[a_ref.child0 == b ? a_ref.child1 : a_ref.child0];
        Int d = b_ref.child0;
        Int e = b_ref.child1;
        if (nodes_[d].height < nodes_[e].height)
            std::swap(d, e);
        Node& d_ref = nodes_[d];
        Node& e_ref = nodes_[e];
        b_ref.child0 = a;
        b_ref.child1 = d;
        b_ref.parent = a_ref.parent;
        a_ref.parent = b;
        if (b_ref.parent != Nil) {
            if (nodes_[b_ref.parent].child0 == a)
                nodes_[b_ref.parent].child0 = b;
            else
                nodes_[b_ref.parent].child1 = b;
        }
        else
            root_ = b;
        if (a_ref.child0 == b)
            a_ref.child0 = e;
        else
            a_ref.child1 = e;
        e_ref.parent = a;
        a_ref.box = c_ref.box | e_ref.box;
        b_ref.box = a_ref.box | d_ref.box;
        a_ref.height = 1 + std::max(c_ref.height, e_ref.height);
        b_ref.height = 1 + std::max(a_ref.height, d_ref.height);
    };
    // Do tree rotation if necessary.
    Node& node_ref = nodes_[node];
    if (node_ref.is_leaf() or //
        node_ref.height < 2)
        return node;
    Int child0 = node_ref.child0;
    Int child1 = node_ref.child1;
    Int imbalance = nodes_[child1].height - nodes_[child0].height;
    if (imbalance > 1) {
        rotate(child1, node);
        return child1;
    }
    if (imbalance < -1) {
        rotate(child0, node);
        return child0;
    }
    return node;
}

template class DynamicBbTree<2>;

template class DynamicBbTree<3>;

} // namespace pre
