#include <pre/memory>
#include <pre-graphics/ImmutableKdTree>

namespace pre {

template <size_t Dim>
class ImmutableKdTreeBuilder {
  public:
    using KdTree = ImmutableKdTree<Dim>;
    using Point = typename KdTree::Point;
    using Box = typename KdTree::Box;
    using Item = typename KdTree::Item;
    using Items = typename KdTree::Items;

    struct Node {
        Point point;
        Node* child0 = nullptr;
        Node* child1 = nullptr;
        ssize_t index = 0;
        ssize_t axis = 0;
    };

    Node* root = nullptr;
    HeapArena<> node_arena = {};

  public:
    /// Build.
    void build(Items& items) {
        root = build_range({items.data(), items.data() + items.size()});
    }

    /// Build range recursively.
    Node* build_range(IteratorRange<Item*> items) {
        if (items.size() == 0)
            return nullptr;
        Node* node = new (node_arena) Node();
        if (items.size() == 1) {
            node->point = items[0].point;
            node->index = items[0].index;
            return node;
        }
        Box box;
        for (Item& item : items)
            box |= item.point;
        int axis = box.extent().argmax();
        Item* middle = items.begin() + items.size() / 2;
        std::nth_element(
                items.begin(), middle, items.end(),
                [&](const Item& lhs, const Item& rhs) {
                    return lhs.point[axis] < rhs.point[axis];
                });
        node->point = middle->point;
        node->index = middle->index;
        node->child0 = build_range({items.begin(), middle});
        node->child1 = build_range({middle + 1, items.end()});
        node->axis = axis;
        return node;
    }

    /// Collapse.
    static void collapse(Node* from, auto& nodes) {
        ASSERT(from);
        auto& node = nodes.emplace_back();
        node.point = from->point;
        node.index = from->index;
        node.right = 0;
        node.left = 0;
        node.axis = from->axis;
        if (from->child0) {
            node.left = 1;
            collapse(from->child0, nodes);
        }
        if (from->child1) {
            node.right = &nodes.back() - &node + 1;
            collapse(from->child1, nodes);
        }
    }
};

template <size_t Dim>
void ImmutableKdTree<Dim>::build(Items& items) {
    // Run builder.
    ImmutableKdTreeBuilder<Dim> builder;
    builder.build(items);

    // Collapse.
    nodes.clear();
    nodes.reserve(items.size());
    ImmutableKdTreeBuilder<Dim>::collapse(builder.root, nodes);
}

template class ImmutableKdTree<2>;

template class ImmutableKdTree<3>;

} // namespace pre
