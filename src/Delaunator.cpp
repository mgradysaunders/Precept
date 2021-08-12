#include <pre-graphics/Delaunator>

namespace pre {

void Delaunator::build() {
    triangles.clear();
    edge_triangles.clear();
    boundary_edges.clear();
    if (points.size() < 3)
        return;

    // Center.
    Pointd center = {};
    for (const Point& point : points)
        if (pre::isfinite(point).all())
            center += point;
    center /= points.size();

    // Sort by distance.
    using PointQueue = std::vector<std::pair<double, Int>>;
    PointQueue point_queue;
    point_queue.reserve(points.size());
    for (const Point& point : points)
        point_queue.emplace_back(
            distance2(point, center), &point - &points[0]);
    std::sort(
        point_queue.begin(), //
        point_queue.end(),   //
        [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });
    auto point_itr = point_queue.begin();

    // Form first triangle.
    Triangle first_tri;
    double first_tri_area = 0;
    while (point_itr + 2 < point_queue.end()) {
        first_tri[0] = point_itr[0].second;
        first_tri[1] = point_itr[1].second;
        first_tri[2] = point_itr[2].second;
        first_tri_area = signed_area(first_tri);
        if (std::fabs(first_tri_area) > 1e-9 * point_queue.back().first)
            break;
        ++point_itr;
    }

    // First triangle initialization failed?
    if (point_itr + 3 >= point_queue.end())
        throw std::runtime_error("first triangle initialization failed!");

    // First triangle area negative?
    if (first_tri_area < 0)
        first_tri.flip_winding();

    // Add first triangle.
    triangles.reserve(points.size());
    triangles.push_back(first_tri);
    Edge first_edges[3] = {
        {first_tri[0], first_tri[1]},
        {first_tri[1], first_tri[2]},
        {first_tri[2], first_tri[0]}};
    edge_triangles[first_edges[0]].push(0);
    edge_triangles[first_edges[1]].push(0);
    edge_triangles[first_edges[2]].push(0);
    boundary_edges.insert(first_edges[0]);
    boundary_edges.insert(first_edges[1]);
    boundary_edges.insert(first_edges[2]);

    // Add point, update boundary.
    auto add_point = [&](Int v) {
        std::set<Edge> boundary_edges_to_add;
        std::vector<Edge> boundary_edges_to_remove;
        if (not pre::isfinite(points[v]).all())
            return; // Ignore
        for (Edge edge : boundary_edges) {
            Int a = edge[0];
            Int b = edge[1];
            if (signed_area(Triangle(v, a, b)) < 0) {
                // Form counter-clockwise triangle.
                Int f = triangles.size();
                edge_triangles[edge].push(f);
                edge_triangles[Edge(v, a)].push(f);
                edge_triangles[Edge(b, v)].push(f);
                triangles.emplace_back(v, b, a);

                // Record boundary updates.
                boundary_edges_to_add.insert(Edge(a, v));
                boundary_edges_to_add.insert(Edge(v, b));
                boundary_edges_to_remove.push_back(edge);
            }
        }
        for (Edge edge : boundary_edges_to_remove)
            boundary_edges.erase(boundary_edges.find(edge));
        for (Edge edge : boundary_edges_to_add)
            if (not edge_triangles[edge].full())
                boundary_edges.insert(edge);
    };

    // ... and flip edges as needed to maintain Delaunay condition.
    auto and_flip_edges = [&] {
        std::set<Edge> flip;
        for (const auto& [edge, tris] : edge_triangles)
            flip.insert(edge);
        while (not flip.empty()) {
            std::set<Edge> next;
            for (Edge edge : flip) {
                EdgeIterator itr = edge_triangles.find(edge);
                if (delaunay_condition(itr))
                    continue;
                // Flip edge.
                Int a = itr->first[0], f = itr->second[0];
                Int b = itr->first[1], g = itr->second[1];
                Int p = triangles[f].opposite(itr->first);
                Int q = triangles[g].opposite(itr->first);
                triangles[f] = Triangle(p, a, q);
                triangles[g] = Triangle(p, q, b);
                edge_triangles.erase(itr);
                edge_triangles[Edge(p, q)] = TrianglePair(f, g);
                edge_triangles[Edge(p, b)].replace(f, g);
                edge_triangles[Edge(q, a)].replace(g, f);
                if (signed_area(triangles[f]) < 0) {
                    triangles[f].flip_winding();
                    triangles[g].flip_winding();
                }

                // Maybe flip neighboring edges on next iteration.
                next.insert(Edge(p, b));
                next.insert(Edge(q, a));
                next.insert(Edge(p, a));
                next.insert(Edge(q, b));
            }
            flip.swap(next);
        }
    };

    // Add remaining points.
    for (point_itr += 3; //
         point_itr < point_queue.end(); point_itr++) {
        add_point(point_itr->second);
        and_flip_edges();
    }

    for (Edge edge : boundary_edges) {
        TrianglePair& tris = edge_triangles.at(edge);
        ASSERT(
            tris[0] == Nil or //
            tris[1] == Nil);

        if (tris[0] == Nil)
            std::swap(tris[0], tris[1]);

        // Force boundary edge to be first edge in triangle
        Triangle& tri = triangles[tris[0]];
        if (tri[0] != edge[0] or //
            tri[1] != edge[1]) {
            tri.cycle();
            if (tri[0] != edge[0] or //
                tri[1] != edge[1])
                tri.cycle();
        }
        ASSERT(tri[0] == edge[0] and tri[1] == edge[1]);
    }
}

} // namespace pre
