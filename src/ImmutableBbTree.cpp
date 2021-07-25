#include <pre-graphics/ImmutableBbTree>

namespace pre {

/// An immutable bounding volume hierarchy builder.
template <size_t Dim>
class ImmutableBbTreeBuilder {
    static_assert(Dim == 2 or Dim == 3);

  public:
    using Bvh = ImmutableBbTree<Dim>;
    using Box = typename Bvh::Box;
    using Item = typename Bvh::Item;
    using Items = typename Bvh::Items;

    struct Node {
        Box box = {};           ///< Box.
        Node* left = nullptr;   ///< If branch, left child.
        Node* right = nullptr;  ///< If branch, right child.
        ssize_t split_axis = 0; ///< If branch, split axis.
        ssize_t first_item = 0; ///< If leaf, first item index.
        ssize_t item_count = 0; ///< If leaf, item count.
    };

    Node* root = {};
    ssize_t leaf_limit = 4;
    ssize_t node_count = 0;
    HeapArena<> node_arena = {};

  public:
    /// Build.
    void build(Items& items) {
        ssize_t first_item = 0;
        root = build_range(
                first_item, {items.data(), items.data() + items.size()});
        ASSERT(first_item == items.size());
    }

    /// Build range recursively.
    Node* build_range(ssize_t& first_item, IteratorRange<Item*> items) {
        Node* node = new (node_arena) Node();
        node_count++;
        Box box;
        Box box_center;
        for (const Item& item : items) {
            box |= item.box;
            box_center |= item.box_center;
        }
        ssize_t split_axis = box_center.extent().argmax();
        ssize_t item_count = items.size();
        if (item_count <= leaf_limit) { // Leaf?
            *node = {box, nullptr, nullptr, 0, first_item, item_count};
            first_item += item_count;
        }
        else {
            Item* split = find_split_sah(box_center, split_axis, items);
            *node = {
                    box,
                    build_range(first_item, {items.begin(), split}),
                    build_range(first_item, {split, items.end()}),
                    split_axis,
                    0,
                    0};
        }
        return node;
    }

    /// Find split using surface area heuristic.
    static Item* find_split_sah(
            const Box& box_center,
            ssize_t split_axis,
            IteratorRange<Item*> items) {
        constexpr ssize_t Nbins = 8;
        using Bin = std::pair<Box, ssize_t>;
        if (box_center.min()[split_axis] == //
            box_center.max()[split_axis])
            return find_split_equal_counts(split_axis, items);

        // Initialize bins.
        std::array<Bin, Nbins> bins = {};
        for (const Item& item : items) {
            float cen = item.box_center[split_axis];
            float cen_min = box_center.min()[split_axis];
            float cen_max = box_center.max()[split_axis];
            int pos = std::min<int>(
                    Nbins * ((cen - cen_min) / (cen_max - cen_min)),
                    Nbins - 1);
            bins[pos].first |= item.box;
            bins[pos].second++;
        }

        // Initialize sweeps.
        std::array<Bin, Nbins - 1> lsweep;
        std::array<Bin, Nbins - 1> rsweep;
        {
            auto itrlsweep = lsweep.begin(), itrlbins = bins.begin();
            auto itrrsweep = rsweep.rbegin(), itrrbins = bins.rbegin();
            *itrlsweep++ = *itrlbins++;
            *itrrsweep++ = *itrrbins++;
            for (; itrlsweep < lsweep.end();
                 ++itrlsweep, ++itrlbins, ++itrrsweep, ++itrrbins) {
                itrlsweep->first = (itrlsweep - 1)->first | itrlbins->first;
                itrrsweep->first = (itrrsweep - 1)->first | itrrbins->first;
                itrlsweep->second = (itrlsweep - 1)->second + itrlbins->second;
                itrrsweep->second = (itrrsweep - 1)->second + itrrbins->second;
            }
        }

        // Compute costs.
        Array<float, Nbins - 1> costs;
        {
            auto itrcosts = costs.begin();
            auto itrlsweep = lsweep.begin();
            auto itrrsweep = rsweep.begin();
            for (; itrrsweep < rsweep.end();
                 ++itrcosts, ++itrlsweep, ++itrrsweep)
                *itrcosts =
                        itrlsweep->first.surface_area() * itrlsweep->second +
                        itrrsweep->first.surface_area() * itrrsweep->second;
        }

        // Partition.
        int costs_argmin = costs.argmin();
        Item* split = std::partition(
                items.begin(), items.end(), [=](const Item& item) {
                    float cen = item.box_center[split_axis];
                    float cen_min = box_center.min()[split_axis];
                    float cen_max = box_center.max()[split_axis];
                    int pos = std::min<int>(
                            Nbins * ((cen - cen_min) / (cen_max - cen_min)),
                            Nbins - 1);
                    return pos <= costs_argmin;
                });

        // Partition successful?
        if (split != items.begin() and //
            split != items.end())
            return split;
        else
            return find_split_equal_counts(split_axis, items);
    }

    /// Find split using equal counts.
    static Item* find_split_equal_counts(
            ssize_t split_axis, IteratorRange<Item*> items) {
        Item* split = items.begin() + items.size() / 2;
        std::nth_element(
                items.begin(), split, items.end(),
                [=](const Item& item0, //
                    const Item& item1) -> bool {
                    return item0.box_center[split_axis] <
                           item1.box_center[split_axis];
                });
        return split;
    }

    /// Collapse.
    static void collapse(Node* from, auto& nodes) {
        ASSERT(from);
        nodes.emplace_back();
        auto& node = nodes.back();
        node.box = from->box;
        if (from->item_count > 0) {
            ASSERT(not from->left);
            ASSERT(not from->right);
            node.first = from->first_item;
            node.count = from->item_count;
            node.split_axis = 0;
        }
        else {
            collapse(from->left, nodes);
            node.right = &nodes.back() - &node + 1;
            node.count = 0;
            node.split_axis = from->split_axis;
            collapse(from->right, nodes);
        }
    }
};

template <size_t Dim>
inline void ImmutableBbTree<Dim>::build(int leaf_limit, Items& items) {
    if (leaf_limit < 1)
        leaf_limit = 1;

    // Run builder.
    ImmutableBbTreeBuilder<Dim> builder;
    builder.leaf_limit = leaf_limit;
    builder.build(items);

    // Collapse.
    nodes.clear();
    nodes.reserve(builder.node_count);
    ImmutableBbTreeBuilder<Dim>::collapse(builder.root, nodes);
}

template class ImmutableBbTree<2>;

template class ImmutableBbTree<3>;

} // namespace pre
