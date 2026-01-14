#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <exception>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

#include <boost/variant/get.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/number_utils.h>
// Support both CGAL 5.x and 4.x APIs
#if CGAL_VERSION_MAJOR >= 5
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#else
#include <CGAL/IO/STL_reader.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#endif
#include <CGAL/squared_distance_3.h>
#include <CGAL/Surface_mesh.h>

namespace {
using Kernel = CGAL::Simple_cartesian<double>;
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

template <typename T, typename V>
auto variant_get_if_impl(const V* v, int) -> decltype(std::get_if<T>(v)) {
    return std::get_if<T>(v);
}

template <typename T, typename V>
const T* variant_get_if_impl(const V* v, long) {
    return boost::get<T>(v);
}

template <typename T, typename V>
const T* variant_get_if(const V* v) {
    return variant_get_if_impl<T>(v, 0);
}

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

static bool face_normal_unit(const SurfaceMesh& mesh, SurfaceMesh::Face_index f,
                             double normal_out[3]) {
    if (normal_out == nullptr) {
        return false;
    }
    if (f == SurfaceMesh::null_face()) {
        return false;
    }

    const auto h0 = halfedge(f, mesh);
    if (h0 == SurfaceMesh::null_halfedge()) {
        return false;
    }

    const auto h1 = next(h0, mesh);
    const auto h2 = next(h1, mesh);

    const Point p0 = mesh.point(source(h0, mesh));
    const Point p1 = mesh.point(target(h0, mesh));
    const Point p2 = mesh.point(target(h1, mesh));

    const auto u = p1 - p0;
    const auto v = p2 - p0;
    const auto n = CGAL::cross_product(u, v);
    const double n2 = CGAL::to_double(n.squared_length());
    if (!(n2 > 0.0)) {
        return false;
    }
    const double inv = 1.0 / std::sqrt(n2);
    normal_out[0] = CGAL::to_double(n.x()) * inv;
    normal_out[1] = CGAL::to_double(n.y()) * inv;
    normal_out[2] = CGAL::to_double(n.z()) * inv;
    return true;
}
}  // namespace

extern "C" void* stl_wall_create(const char* filename, double scale_to_m) {
    try {
        if (filename == nullptr || std::strlen(filename) == 0) {
            return nullptr;
        }

        SurfaceMesh mesh;
#if CGAL_VERSION_MAJOR >= 5
        if (!CGAL::IO::read_polygon_mesh(filename, mesh)) {
            return nullptr;
        }
#else
        // CGAL 4.x API for reading STL
        std::vector<std::array<double, 3>> points;
        std::vector<std::array<int, 3>> triangles;
        std::ifstream input(filename, std::ios::binary);
        if (!input || !CGAL::read_STL(input, points, triangles)) {
            return nullptr;
        }
        std::vector<Point> pts;
        pts.reserve(points.size());
        for (const auto& p : points) {
            pts.emplace_back(p[0], p[1], p[2]);
        }
        std::vector<std::vector<std::size_t>> faces;
        faces.reserve(triangles.size());
        for (const auto& t : triangles) {
            faces.push_back({static_cast<std::size_t>(t[0]),
                           static_cast<std::size_t>(t[1]),
                           static_cast<std::size_t>(t[2])});
        }
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(pts, faces, mesh);
#endif
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
        if (const Point* ipoint = variant_get_if<Point>(&variant_geom)) {
            cand = *ipoint;
        } else if (const Segment* iseg = variant_get_if<Segment>(&variant_geom)) {
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

extern "C" int stl_wall_first_hit_segment_with_normal(void* h, const double p0_m[3],
                                                      const double p1_m[3],
                                                      double hit_m[3],
                                                      double normal[3]) {
    if (h == nullptr || p0_m == nullptr || p1_m == nullptr || hit_m == nullptr ||
        normal == nullptr) {
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
    SurfaceMesh::Face_index best_face = SurfaceMesh::null_face();

    for (const auto& entry : intersections) {
        const auto& variant_geom = entry.first;
        const auto face = entry.second;
        Point cand;
        if (const Point* ipoint = variant_get_if<Point>(&variant_geom)) {
            cand = *ipoint;
        } else if (const Segment* iseg = variant_get_if<Segment>(&variant_geom)) {
            cand = first_point_along_segment(*iseg, p0);
        } else {
            continue;
        }

        const double d2 = CGAL::squared_distance(p0, cand);
        if (d2 < best_d2) {
            best_d2 = d2;
            best_hit = cand;
            best_face = face;
        }
    }

    if (!std::isfinite(best_d2)) {
        return 0;
    }

    double n_unit[3] = {0.0, 0.0, 0.0};
    if (!face_normal_unit(wall->mesh, best_face, n_unit)) {
        return 0;
    }

    hit_m[0] = best_hit.x();
    hit_m[1] = best_hit.y();
    hit_m[2] = best_hit.z();

    normal[0] = n_unit[0];
    normal[1] = n_unit[1];
    normal[2] = n_unit[2];
    return 1;
}
