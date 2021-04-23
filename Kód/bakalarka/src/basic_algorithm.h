#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H

#include "function.h"
#include "mesh.h"
#include "triangle.h"

#include <vector>

class BasicAlgorithm {
public:
  BasicAlgorithm(Function f, Triangle seed_triangle, numeric e_size,
                 realsymbol x, realsymbol y, realsymbol z)
      : F(f), active_edges(), checked_edges(), my_mesh(seed_triangle),
        e_size(e_size), x(x), y(y), z(z) {
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
  pair<std::optional<Point>, std::optional<Point>>
  find_closest_prev_next(const Edge &working_edge,
                       const vector<Point> &prev, const vector<Point> &next);
  pair<Point, Point> find_prev_next(const Edge &working_edge);
  Point get_projected(const Edge &working_edge);

private:
  Function F;
  vector<Edge> active_edges;
  vector<Edge> checked_edges;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
};


#endif