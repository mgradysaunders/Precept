#include <pre-graphics/AdjacencyLists>

#include "GrowPoolVector.h"

namespace pre {

typename AdjacencyLists::Int AdjacencyLists::private_allocate_node() {
    if (free_node_ == Nil)
        free_node_ = GrowPoolVector(nodes_, &Node::next);
    Int node = free_node_;
    free_node_ = nodes_[node].next;
    nodes_[node] = Node();
    return node;
}

typename AdjacencyLists::Int AdjacencyLists::private_allocate_link() {
    if (free_link_ == Nil)
        free_link_ = GrowPoolVector(links_, &Link::next);
    Int link = free_link_;
    free_link_ = links_[link].next;
    links_[link] = Link();
    return link;
}

void AdjacencyLists::private_connect(Int node0, Int node1, Int Node::*list) {
    ASSERT(node0 < Int(nodes_.size()));
    ASSERT(node1 < Int(nodes_.size()));
    if (node0 == Nil or node1 == Nil)
        return;
    Node& node0_ref = nodes_[node0];
    if (node0_ref.*list == Nil) {
        node0_ref.*list = private_allocate_link();
        links_[node0_ref.*list].node = node1;
        return;
    }
    Int tail = node0_ref.*list;
    while (1) {
        Int next = links_[tail].next;
        if (next == Nil)
            break;
        tail = next;
        if (links_[tail].node == node1)
            return;
    }
    Int link = private_allocate_link();
    links_[tail].next = link;
    links_[link].node = node1;
}

void AdjacencyLists::private_disconnect(
    Int node0, Int node1, Int Node::*list) {
    ASSERT(node0 < Int(nodes_.size()));
    ASSERT(node1 < Int(nodes_.size()));
    if (node0 == Nil or node1 == Nil)
        return;
    Int prev = Nil;
    Int link = nodes_[node0].*list;
    while (link != Nil) {
        if (links_[link].node == node1)
            break;
        prev = link;
        link = links_[link].next;
    }
    if (link != Nil) {
        Int next = links_[link].next;
        if (prev != Nil)
            links_[prev].next = next;
        if (nodes_[node0].*list == link)
            nodes_[node0].*list = next;
        private_deallocate_link(link);
    }
}

} // namespace pre
