#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "algorithms.h"
#include "edge.h"
#include "function.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>
#include <optional>
#include <vector>

using namespace GiNaC;
using namespace std;

class BoundingBox {
public:
  BoundingBox(numeric _min_x, numeric _max_x, numeric _min_y, numeric _max_y,
              numeric _min_z, numeric _max_z);

  //vector<Edge> bounding_edges;

  numeric min_x() const;
  numeric max_x() const;
  numeric min_y() const;
  numeric max_y() const;
  numeric min_z() const;
  numeric max_z() const;

  /*
  Vector normal_x() const;
  Vector normal_y() const;
  Vector normal_z() const;
  */

  bool in_interval_x(const numeric x) const;
  bool in_interval_y(const numeric y) const;
  bool in_interval_z(const numeric z) const;

  bool is_inside(const Point P) const;
  bool is_on(const Point P) const;
  bool new_bounding_edge(const Edge &e) const;
  int faces(const Point &P) const;

  std::optional<Point> project_on_min_x(const Point &midpoint, const Point &P) const;
  std::optional<Point> project_on_max_x(const Point &midpoint, const Point &P) const;
  std::optional<Point> project_on_min_y(const Point &midpoint, const Point &P) const;
  std::optional<Point> project_on_max_y(const Point &midpoint, const Point &P) const;
  std::optional<Point> project_on_min_z(const Point &midpoint, const Point &P) const;
  std::optional<Point> project_on_max_z(const Point &midpoint, const Point &P) const;

  Point project_on_box(const Point &midpoint, const Point &P) const;
  Point crop_to_box(const Point &midpoint, const Point &P,
                    const numeric &e_size, const Function &F) const;

private:
  numeric _min_x;
  numeric _max_x;
  numeric _min_y;
  numeric _max_y;
  numeric _min_z;
  numeric _max_z;
};

#endif
