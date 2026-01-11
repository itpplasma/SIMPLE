#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <iterator>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Surface_mesh.h>

namespace {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Segment = Kernel::Segment_3;
using SurfaceMesh = CGAL::Surface_mesh<Point>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh>;
using Traits = CGAL::AABB_traits<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

struct Wall {
    SurfaceMesh mesh;
    Tree tree;

    explicit Wall(SurfaceMesh&& mesh_in)
        : mesh(std::move(mesh_in)),
          tree(faces(mesh).begin(), faces(mesh).end(), mesh) {
        tree.accelerate_distance_queries();
    }
};

static Point to_point(const double p[3]) { return Point(p[0], p[1], p[2]); }

static void scale_mesh(SurfaceMesh& mesh, double s) {
    if (!(s > 0.0)) {
        throw std::runtime_error("scale must be > 0");
    }
    if (s == 1.0) {
        return;
    }
    for (auto v : mesh.vertices()) {
        const Point p = mesh.point(v);
        mesh.point(v) = Point(p.x() * s, p.y() * s, p.z() * s);
    }
}

static Point first_point_along_segment(const Segment& seg, const Point& p0) {
    const Point a = seg.source();
    const Point b = seg.target();
    const double da = CGAL::squared_distance(p0, a);
    const double db = CGAL::squared_distance(p0, b);
    return (da <= db) ? a : b;
}
}  // namespace

extern "C" void* stl_wall_create(const char* filename, double scale_to_m) {
    try {
        if (filename == nullptr || std::strlen(filename) == 0) {
            return nullptr;
        }

        SurfaceMesh mesh;
        if (!CGAL::IO::read_polygon_mesh(filename, mesh)) {
            return nullptr;
        }
        if (num_faces(mesh) == 0 || num_vertices(mesh) == 0) {
            return nullptr;
        }

        scale_mesh(mesh, scale_to_m);

        auto* wall = new Wall(std::move(mesh));
        return static_cast<void*>(wall);
    } catch (...) {
        return nullptr;
    }
}

extern "C" void stl_wall_destroy(void* h) {
    if (h == nullptr) {
        return;
    }
    auto* wall = static_cast<Wall*>(h);
    delete wall;
}

extern "C" int stl_wall_first_hit_segment(void* h, const double p0_m[3],
                                          const double p1_m[3],
                                          double hit_m[3]) {
    if (h == nullptr || p0_m == nullptr || p1_m == nullptr || hit_m == nullptr) {
        return 0;
    }

    const auto* wall = static_cast<const Wall*>(h);
    const Point p0 = to_point(p0_m);
    const Point p1 = to_point(p1_m);
    const Segment seg(p0, p1);

    using Intersection = Tree::Intersection_and_primitive_id<Segment>::Type;
    std::vector<Intersection> intersections;
    try {
        wall->tree.all_intersections(seg, std::back_inserter(intersections));
    } catch (...) {
        return 0;
    }

    if (intersections.empty()) {
        return 0;
    }

    double best_d2 = std::numeric_limits<double>::infinity();
    Point best_hit(0.0, 0.0, 0.0);

    for (const auto& entry : intersections) {
        const auto& variant_geom = entry.first;
        Point cand;
        if (const Point* ipoint = std::get_if<Point>(&variant_geom)) {
            cand = *ipoint;
        } else if (const Segment* iseg = std::get_if<Segment>(&variant_geom)) {
            cand = first_point_along_segment(*iseg, p0);
        } else {
            continue;
        }

        const double d2 = CGAL::squared_distance(p0, cand);
        if (d2 < best_d2) {
            best_d2 = d2;
            best_hit = cand;
        }
    }

    if (!std::isfinite(best_d2)) {
        return 0;
    }

    hit_m[0] = best_hit.x();
    hit_m[1] = best_hit.y();
    hit_m[2] = best_hit.z();
    return 1;
}
