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
  void starting();
  void ending();
  bool step(const Edge &working_edge);
  bool fix_proj(const Edge &working_edge);
  bool fix_prev_next(const Edge &working_edge, const bool is_prev);
  bool fix_overlap(const Edge &working_edge, Point overlap_point);
  int fix_holes(const Edge &working_edge);
  int update_border(const Edge &new_edge1, const Edge &new_edge2);
  bool basic_triangle(const Edge &working_edge,
                      const Triangle &neighbour_triangle);
  bool create_triangle(const Edge &working_edge, const Point &P);
  bool Delaunay_conditions(const Edge &working_edge, const Point &P,
                           const Triangle &neighbour_triangle);
  Point get_projected(const Edge &working_edge);
  bool check_overlap_normal(const Point candidate, const Point prev,
                            const Point next, const Edge &working_edge);

  void add_marks();

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