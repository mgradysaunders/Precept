#include <pre-graphics/DynamicKdTree>

namespace pre {

template <size_t Dim>
typename DynamicKdTree<Dim>::Box DynamicKdTree<Dim>::region(Int node) const {
    Box box = box_;
    Array<bool, Dim> min_done = {};
    Array<bool, Dim> max_done = {};
    Point point = nodes_[node].point;
    Int parent = Nil;
    for (parent = nodes_[node].parent; parent != Nil;
         parent = nodes_[parent].parent) {
        const Node& ref = nodes_[parent];
        if (point < ref) {
            minimize(box[1][ref.axis], ref.threshold());
            max_done[ref.axis] = true;
        }
        else {
            maximize(box[0][ref.axis], ref.threshold());
            min_done[ref.axis] = true;
        }
        if (min_done.all() and max_done.all())
            break;
    }
    return box;
}

template <size_t Dim>
typename DynamicKdTree<Dim>::Nearest DynamicKdTree<Dim>::nearest(
        Point point) const {
    Nearest near;
    GrowableStack<Int> todo;
    if (root_ != Nil)
        todo.push(root_);
    while (not todo.empty()) {
        Int node = todo.pop();
        const Node& node_ref = nodes_[node];
        if (not node_ref.dead) {
            float dist2 = distance2(node_ref.point, point);
            if (near.dist2 > dist2) {
                near.dist2 = dist2;
                near.node = node;
            }
        }
        float diff = node_ref.point[node_ref.axis] - point[node_ref.axis];
        Int child0 = node_ref.child0;
        Int child1 = node_ref.child1;
        if (diff > 0) {
            diff = -diff;
            std::swap(child0, child1);
        }
        if (child0 != Nil and (diff > 0 or diff * diff < near.dist2))
            todo.push(child0);
        if (child1 != Nil and (diff < 0 or diff * diff < near.dist2))
            todo.push(child1);
    }
    return near;
}
template <size_t Dim>
void DynamicKdTree<Dim>::nearest(
        Point point, IteratorRange<Nearest*> near) const {
    if (near.size() == 0)
        return;
    if (near.size() == 1) {
        near[0] = nearest(point);
        return;
    }
    near.fill(Nearest());
    Nearest* heap = near.begin();
    GrowableStack<Int> todo;
    if (root_ != Nil)
        todo.push(root_);
    while (not todo.empty()) {
        Int node = todo.pop();
        const Node& node_ref = nodes_[node];
        float dist2 = distance2(node_ref.point, point);
        if (heap != near.end() or //
            dist2 < near[0].dist2) {
            if (heap == near.end())
                std::pop_heap(near.begin(), heap--);
            heap->node = node;
            heap->dist2 = dist2;
            std::push_heap(near.begin(), ++heap);
        }
        float diff = node_ref.point[node_ref.axis] - point[node_ref.axis];
        Int child0 = node_ref.child0;
        Int child1 = node_ref.child1;
        if (diff > 0) {
            diff = -diff;
            std::swap(child0, child1);
        }
        if (child0 != Nil and
            (heap != near.end() or diff > 0 or diff * diff < near[0].dist2))
            todo.push(child0);
        if (child1 != Nil and
            (heap != near.end() or diff < 0 or diff * diff < near[0].dist2))
            todo.push(child1);
    }
    std::sort_heap(near.begin(), heap);
}

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
    node_count_++;
    nodes_[node] = Node();
    return node;
}

template <size_t Dim>
typename DynamicKdTree<Dim>::Int DynamicKdTree<Dim>::private_select_axis(
        Int node) {
    Node& node_ref = nodes_[node];
    Point cost = Point(Inf<float>);
    for (Int axis = 0; axis < Int(Dim); axis++) {
        node_ref.axis = axis;
        Box box = region(node);
        Box box0 = box;
        Box box1 = box;
        box0[0][axis] = node_ref.threshold();
        box1[1][axis] = node_ref.threshold();
        float cost0 = box0.surface_area() / box0.volume();
        float cost1 = box1.surface_area() / box1.volume();
        if (std::isfinite(cost0) and //
            std::isfinite(cost1))
            cost[axis] = std::fmax(cost0, cost1);
        else
            cost[axis] = std::isfinite(cost0) ? cost0 : cost1;
    }
    return cost.argmin();
}

template <size_t Dim>
void DynamicKdTree<Dim>::private_insert(Int leaf) {
    box_ |= nodes_[leaf].point;
    if (root_ == Nil) {
        root_ = leaf;
        nodes_[root_].parent = Nil;
        return;
    }
    Int node = root_;
    while (1) {
        Node& leaf_ref = nodes_[leaf];
        Node& node_ref = nodes_[node];
        if (node_ref.axis == -1)
            node_ref.axis = private_select_axis(node);
        Int* child0 = &node_ref.child0;
        Int* child1 = &node_ref.child1;
        Int* child = leaf_ref.point < node_ref ? child0 : child1;
        if (*child == Nil) {
            *child = leaf;
            leaf_ref.parent = node;
            break;
        }
        node = *child;
    }
    // Update heights.
    Int imbalance = 0;
    node = nodes_[leaf].parent;
    while (node != Nil) {
        Node& node_ref = nodes_[node];
        Int child0 = node_ref.child0;
        Int child1 = node_ref.child1;
        Int height0 = child0 != Nil ? nodes_[child0].height : 0;
        Int height1 = child1 != Nil ? nodes_[child1].height : 0;
        maximize(imbalance, std::abs(height1 - height0));
        node_ref.height = std::max(height0, height1);
        node_ref.height++;
        node = node_ref.parent;
    }
    if (imbalance > 4)
        private_rebalance();
}

template <size_t Dim>
bool DynamicKdTree<Dim>::private_remove(Int node) {
    Node& node_ref = nodes_[node];
    Int child0 = node_ref.child0;
    Int child1 = node_ref.child1;
    if (child0 != Nil and child1 != Nil) {
        node_ref.dead = true;
        dead_count_++;
        if (dead_count_ > node_count_ / 2)
            private_rebalance();
        return false;
    }
    Int child = child0 != Nil ? child0 : child1;
    Int parent = node_ref.parent;
    if (parent != Nil) {
        Node& parent_ref = nodes_[parent];
        if (parent_ref.child0 == node)
            parent_ref.child0 = child;
        else
            parent_ref.child1 = child;
    }
    if (child != Nil)
        nodes_[child].parent = parent;
    // Update heights.
    node = parent;
    while (node != Nil) {
        Node& node_ref = nodes_[node];
        child0 = node_ref.child0;
        child1 = node_ref.child1;
        Int height0 = child0 != Nil ? nodes_[child0].height : 0;
        Int height1 = child1 != Nil ? nodes_[child1].height : 0;
        node_ref.height = std::max(height0, height1);
        node_ref.height++;
        node = node_ref.parent;
    }
    return true;
}

template <size_t Dim>
void DynamicKdTree<Dim>::private_rebalance() {
    std::vector<Int> nodes;
    nodes.reserve(nodes_.size());
    box_ = {};
    for (Int node = 0; node < Int(nodes_.size()); node++) {
        if (nodes_[node].dead) { // Dead?
            private_deallocate(node);
            continue;
        }
        Node& node_ref = nodes_[node];
        if (node_ref.height >= 0) {
            node_ref.parent = Nil;
            node_ref.child0 = Nil;
            node_ref.child1 = Nil;
            node_ref.height = 0;
            node_ref.axis = Nil;
            node_ref.dead = false;
            nodes.push_back(node);
            box_ |= node_ref.point;
        }
    }
    dead_count_ = 0;
    root_ = private_rebalance({nodes.data(), nodes.data() + nodes.size()});
    rebalance_count_++;
}

template <size_t Dim>
typename DynamicKdTree<Dim>::Int DynamicKdTree<Dim>::private_rebalance(
        IteratorRange<Int*> nodes) {
    if (nodes.size() == 0)
        return Nil;
    if (nodes.size() == 1)
        return nodes[0];
    Box box;
    for (Int node : nodes)
        box |= nodes_[node].point;
    Int axis = box.extent().argmax();
    Int* middle = nodes.begin() + nodes.size() / 2;
    std::nth_element(
            nodes.begin(), middle, nodes.end(), [&](Int lhs, Int rhs) {
                return nodes_[lhs].point[axis] < nodes_[rhs].point[axis];
            });
    Int child0 = private_rebalance({nodes.begin(), middle});
    Int child1 = private_rebalance({middle + 1, nodes.end()});
    Node& middle_ref = nodes_[*middle];
    middle_ref.axis = axis;
    middle_ref.child0 = child0;
    middle_ref.child1 = child1;
    middle_ref.height = 0;
    Node& child0_ref = nodes_[child0];
    Node& child1_ref = nodes_[child1];
    if (child0 != Nil) {
        child0_ref.parent = *middle;
        middle_ref.height = 1 + child0_ref.height;
    }
    if (child1 != Nil) {
        child1_ref.parent = *middle;
        if (middle_ref.height < 1 + child1_ref.height)
            middle_ref.height = 1 + child1_ref.height;
    }
    return *middle;
}

template class DynamicKdTree<2>;

template class DynamicKdTree<3>;

} // namespace pre
