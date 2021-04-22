#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "vector.h"
#include "edge.h"
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

  vector<Edge> bounding_edges;

  numeric min_x() const;
  numeric max_x() const;
  numeric min_y() const;
  numeric max_y() const;
  numeric min_z() const;
  numeric max_z() const;

  Vector normal_x() const;
  Vector normal_y() const;
  Vector normal_z() const;

  bool in_interval_x(const numeric x) const;
  bool in_interval_y(const numeric y) const;
  bool in_interval_z(const numeric z) const;

  bool is_inside(const Point P) const;
  bool is_on(const Point P) const;

  // returns set of close walls, 1 for min_x wall, 2 for max_x, 3 for
  // min_y, 4 for max_y, 5 for min_z, 6 for max_z
  std::set<int> close_walls(const Point P, numeric e_size) const;
  Point crop_to_box(Point P, const std::set<int> & close_walls);

private:
  numeric _min_x;
  numeric _max_x;
  numeric _min_y;
  numeric _max_y;
  numeric _min_z;
  numeric _max_z;
};



#endif