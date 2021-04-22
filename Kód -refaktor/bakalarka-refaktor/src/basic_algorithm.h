#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H

#include "function.h"
#include "mesh.h"
#include "triangle.h"

#include <set>

class BasicAlgorithm {
public:
  BasicAlgorithm(Function f, Triangle seed_triangle, numeric e_size,
                 realsymbol x, realsymbol y, realsymbol z, BoundingBox b)
      : F(f), active_edges(), checked_edges(), my_mesh(seed_triangle),
        e_size(e_size), x(x), y(y), z(z), bounding_box(b) {
    active_edges.push_back(seed_triangle.AB());
    active_edges.push_back(seed_triangle.BC());
    active_edges.push_back(seed_triangle.CA());
  }

  Mesh calculate();
  int update_border(const Edge &new_edge1, const Edge &new_edge2);
  bool fix_prev_next(const Edge &working_edge, const bool is_prev);
  bool fix_overlap(const Edge &working_edge, Point overlap_point);
  bool fix_proj(const Edge &working_edge);
  bool step(const Edge &working_edge);
  void add_marks();
  int fix_holes(const Edge &working_edge);
  void ending();
  void starting();
  bool is_active(const Edge &edge) const;
  bool is_checked(const Edge &edge) const;
  bool is_bounding(const Edge &edge) const;
  bool is_border(const Edge &edge) const;
  bool is_border_point(Point P) const;
  void delete_from_active(const Edge &edge);
  void delete_from_checked(const Edge &edge);
  void delete_from_bounding(const Edge &edge);
  void push_edge_to_active(const Edge &edge);
  void push_edge_to_checked(const Edge &edge);
  void push_edge_to_bounding(const Edge &edge);
  bool good_edges(const Edge &working_edge, const Point &P);
  std::optional<Point> get_closest_point(const Edge &working_edge, const Triangle &N);
  std::optional<pair<Edge, numeric>> get_closest_edge(const Point &P, const Triangle &N);
  pair<std::optional<Point>, std::optional<Point>> find_closest_prev_next(const Edge &working_edge,
  const vector<Point> &prev, const vector<Point> &next);
  pair<Point, Point> find_prev_next(const Edge &working_edge);
  bool is_vertex_good_possibility(const Point candidate, const Point prev,
                                const Point next, const Edge &working_edge,
                                const Triangle &neighbour_triangle);

  Point get_projected(const Edge &working_edge);
  std::optional<vector<Point>>
  empty_surrounding(Point P, const Edge working_edge, const Triangle N) const;
  







private:
  Function F;
  vector<Edge> active_edges;
  vector<Edge> checked_edges;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
  BoundingBox bounding_box;
};

#endif