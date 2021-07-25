#include <pre-graphics/DynamicBbTree>

namespace pre {

#if 0
template <size_t Dim>
void DynamicBbTree<Dim>::rebuild_optimal() {
    if (root_ == Nil)
        return;
    std::vector<Int> todo;
    todo.reserve(nodes_.size());
    for (Int node = 0; node < Int(nodes_.size()); node++) {
        if (nodes_[node].height < 0)
            continue;
        if (nodes_[node].is_leaf()) {
            nodes_[node].parent = Nil;
            todo.push_back(node);
        }
        else {
            private_deallocate(node);
        }
    }
    while (todo.size() > 1) {
        float min_cost = INFINITY;
        Int mini = -1;
        Int minj = -1;
        for (Int i = 0; i < Int(todo.size()); ++i) {
            for (Int j = i + 1; j < Int(todo.size()); ++j) {
                Box boxi = nodes_[todo[i]].box;
                Box boxj = nodes_[todo[j]].box;
                float cost = (boxi | boxj).surface_area();
                if (min_cost > cost)
                    min_cost = cost, mini = i, minj = j;
            }
        }
        Int parent = private_allocate();
        Int child0 = todo[mini];
        Int child1 = todo[minj];
        Node& child0_ref = nodes_[child0];
        Node& child1_ref = nodes_[child1];
        Node& parent_ref = nodes_[parent];
        parent_ref.child0 = child0;
        parent_ref.child1 = child1;
        parent_ref.height = pre::max(child0_ref.height, child1_ref.height);
        parent_ref.height++;
        parent_ref.box = child0_ref.box | child1_ref.box;
        parent_ref.parent = Nil;
        child0_ref.parent = parent;
        child1_ref.parent = parent;
        todo[minj] = todo.back();
        todo[mini] = parent;
        todo.pop_back();
    }
    root_ = todo[0];
}
#endif

template <size_t Dim>
typename DynamicBbTree<Dim>::Int DynamicBbTree<
        Dim>::private_allocate() {
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
void DynamicBbTree<Dim>::private_deallocate(Int node) {
    nodes_[node].next = free_;
    nodes_[node].height = -1;
    free_ = node;
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
    new_parent_ref.parent = old_parent;
    new_parent_ref.box = leaf_box | nodes_[node].box;
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

    node = nodes_[leaf].parent;
    while (node != Nil) {
        node = private_balance(node);
        Node& node_ref = nodes_[node];
        Node& child0_ref = nodes_[node_ref.child0];
        Node& child1_ref = nodes_[node_ref.child1];
        node_ref.height = 1 + pre::max(child0_ref.height, child1_ref.height);
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
            node_ref.height = pre::max(child0_ref.height, child1_ref.height);
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
        Int a) {
    Node& a_ref = nodes_[a];
    if (a_ref.is_leaf() or a_ref.height < 2)
        return a;
    Int b = a_ref.child0;
    Int c = a_ref.child1;
    Node& b_ref = nodes_[b];
    Node& c_ref = nodes_[c];
    Int imbalance = c_ref.height - b_ref.height;
    // Rotate C up to A.
    if (imbalance > 1) {
        Int f = c_ref.child0;
        Int g = c_ref.child1;
        Node& f_ref = nodes_[f];
        Node& g_ref = nodes_[g];
        // Swap A and C.
        c_ref.child0 = a;
        c_ref.parent = a_ref.parent;
        a_ref.parent = c;
        if (c_ref.parent != Nil) {
            if (nodes_[c_ref.parent].child0 == a)
                nodes_[c_ref.parent].child0 = c;
            else
                nodes_[c_ref.parent].child1 = c;
        }
        else
            root_ = c;
        // Rotate.
        if (f_ref.height > g_ref.height) {
            c_ref.child1 = f;
            a_ref.child1 = g;
            g_ref.parent = a;
            a_ref.box = b_ref.box | g_ref.box;
            c_ref.box = a_ref.box | f_ref.box;
            a_ref.height = 1 + pre::max(b_ref.height, g_ref.height);
            c_ref.height = 1 + pre::max(a_ref.height, f_ref.height);
        }
        else {
            c_ref.child1 = g;
            a_ref.child1 = f;
            f_ref.parent = a;
            a_ref.box = b_ref.box | f_ref.box;
            c_ref.box = a_ref.box | g_ref.box;
            a_ref.height = 1 + pre::max(b_ref.height, f_ref.height);
            c_ref.height = 1 + pre::max(a_ref.height, g_ref.height);
        }
        return c;
    }
    // Rotate B up to A.
    if (imbalance < -1) {
        Int d = b_ref.child0;
        Int e = b_ref.child1;
        Node& d_ref = nodes_[d];
        Node& e_ref = nodes_[e];
        // Swap A and B.
        b_ref.child0 = a;
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
        // Rotate.
        if (d_ref.height > e_ref.height) {
            b_ref.child1 = d;
            a_ref.child0 = e;
            e_ref.parent = a;
            a_ref.box = c_ref.box | e_ref.box;
            b_ref.box = a_ref.box | d_ref.box;
            a_ref.height = 1 + pre::max(c_ref.height, e_ref.height);
            b_ref.height = 1 + pre::max(a_ref.height, d_ref.height);
        }
        else {
            b_ref.child1 = e;
            a_ref.child0 = d;
            d_ref.parent = a;
            a_ref.box = c_ref.box | d_ref.box;
            b_ref.box = a_ref.box | e_ref.box;
            a_ref.height = 1 + pre::max(c_ref.height, d_ref.height);
            b_ref.height = 1 + pre::max(a_ref.height, e_ref.height);
        }
        return b;
    }
    return a;
}

template class DynamicBbTree<2>;

template class DynamicBbTree<3>;

} // namespace pre
