#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H

#include "function.h"
#include "mesh.h"
#include "triangle.h"
#include "bounding_box.h"

#include <vector>

class BasicAlgorithm {
public:
  BasicAlgorithm(Function f, Triangle seed_triangle, numeric e_size,
                 realsymbol x, realsymbol y, realsymbol z, BoundingBox bounding_box)
      : F(f), active_edges(), checked_edges(), my_mesh(seed_triangle),
        e_size(e_size), x(x), y(y), z(z), bounding_box(bounding_box) {
    if(bounding_box.new_bounding_edge(seed_triangle.AB()))
      bounding_edges.push_back(seed_triangle.AB());
    else
      active_edges.push_back(seed_triangle.AB());
    
    if(bounding_box.new_bounding_edge(seed_triangle.BC()))
      bounding_edges.push_back(seed_triangle.BC());
    else
      active_edges.push_back(seed_triangle.BC());

    if(bounding_box.new_bounding_edge(seed_triangle.CA()))
      bounding_edges.push_back(seed_triangle.CA());
    else
      active_edges.push_back(seed_triangle.CA());
    
  }

  Mesh calculate();
  void starting();
  void ending();
  bool step(const Edge &working_edge);
  bool fix_proj(const Edge &working_edge, const Triangle &neighbour_triangle);
  bool fix_prev_next(const Edge &working_edge,
                     const Triangle &neighbour_triangle, const bool is_prev);
  bool fix_overlap(const Edge &working_edge, const Triangle &neighbour_triangle,
                   Point overlap_point);
  int fix_holes(const Edge &working_edge, const Triangle &neighbour_triangle);
  pair<std::optional<Point>, std::optional<Point>>
  find_closest_prev_next(const Edge &working_edge,
                         const Triangle &neighbour_triangle,
                         const vector<Point> &prev, const vector<Point> &next);
  pair<Point, Point> find_prev_next(const Edge &working_edge,
                                    const Triangle &neighbour_triangle);
  Point get_projected(const Edge &working_edge,
                      const Triangle &neighbour_triangle);
  bool overlap_normals_check(const Point candidate, const Point prev,
                             const Point next, const Edge &working_edge,
                             const Triangle &neighbour_triangle);
  void add_marks();
  int update_border(const Edge &new_edge1, const Edge &new_edge2);
  bool Delaunay_conditions(const Edge &working_edge, const Point &P,
                           const Triangle &neighbour_triangle);

  void create_triangle(const Edge &working_edge, const Point &P);
  bool good_edges(const Edge &working_edge, const Point &P);
  bool basic_triangle(const Edge &working_edge,
                      const Triangle &neighbour_triangle, const Point &prev,
                      const Point &next);

  bool is_active(const Edge &edges);
  bool is_checked(const Edge &edge);
  bool is_bounding(const Edge &edge);
  bool is_border(const Edge &edge);
  bool is_border_point(Point P);
  void delete_from_active(const Edge &edge);
  void delete_from_checked(const Edge &edge);
  void push_edge_to_active(const Edge &edge);
  void push_edge_to_checked(const Edge &edge);

  std::optional<Point> get_closest_point(const Edge &working_edge,
                                         const Triangle &N);
  std::optional<pair<Edge, numeric>> get_closest_edge(const Point &P,
                                                      const Triangle &N);

private:
  Function F;
  vector<Edge> active_edges;
  vector<Edge> checked_edges;
  vector<Edge> bounding_edges;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
  BoundingBox bounding_box;
};

#endif