#include "bounding_box.h"

BoundingBox::BoundingBox(numeric _min_x, numeric _max_x, numeric _min_y,
                         numeric _max_y, numeric _min_z, numeric _max_z)
    : bounding_edges(), _min_x(_min_x), _max_x(_max_x), _min_y(_min_y),
      _max_y(_max_y), _min_z(_min_z), _max_z(_max_z){};

numeric BoundingBox::min_x() const { return _min_x; }
numeric BoundingBox::max_x() const { return _max_x; }
numeric BoundingBox::min_y() const { return _min_y; }
numeric BoundingBox::max_y() const { return _max_y; }
numeric BoundingBox::min_z() const { return _min_z; }
numeric BoundingBox::max_z() const { return _max_z; }

Vector BoundingBox::normal_x() const { return Vector(1, 0, 0); }
Vector BoundingBox::normal_y() const { return Vector(0, 1, 0); }
Vector BoundingBox::normal_z() const { return Vector(0, 0, 1); }

bool BoundingBox::in_interval_x(const numeric x) const {
  return (x >= _min_x && x <= _max_x);
}
bool BoundingBox::in_interval_y(const numeric y) const {
  return (y >= _min_y && y <= _max_y);
}
bool BoundingBox::in_interval_z(const numeric z) const {
  return (z >= _min_z && z <= _max_z);
}

bool BoundingBox::is_inside(const Point P) const {
  return in_interval_x(P.x()) && in_interval_y(P.y()) && in_interval_z(P.z());
}

bool BoundingBox::is_on(const Point P) const {
  bool x_wall =
      ((abs(P.x() - _min_x) < 10e-10 || abs(P.x() - _max_x) < 10e-10) &&
       in_interval_y(P.y()) && in_interval_z(P.z()));
  bool y_wall =
      ((abs(P.y() - _min_y) < 10e-10 || abs(P.y() - _max_y) < 10e-10) &&
       in_interval_x(P.x()) && in_interval_z(P.z()));
  bool z_wall =
      ((abs(P.z() - _min_z) < 10e-10 || abs(P.z() - _max_z) < 10e-10) &&
       in_interval_y(P.y()) && in_interval_x(P.x()));
  return (x_wall || y_wall || z_wall);
}

std::set<int> BoundingBox::close_walls(const Point P, numeric e_size) const {
  std::set<int> result;
  bool x_min_wall = (abs(P.x() - _min_x) < e_size / 3);
  bool x_max_wall = (abs(P.x() - _max_x) < e_size / 3);
  bool y_min_wall = (abs(P.y() - _min_y) < e_size / 3);
  bool y_max_wall = (abs(P.y() - _max_y) < e_size / 3);
  bool z_min_wall = (abs(P.z() - _min_z) < e_size / 3);
  bool z_max_wall = (abs(P.z() - _max_z) < e_size / 3);

  if (x_min_wall)
    result.insert(1);
  if (x_max_wall)
    result.insert(2);
  if (y_min_wall)
    result.insert(3);
  if (y_max_wall)
    result.insert(4);
  if (z_min_wall)
    result.insert(5);
  if (z_max_wall)
    result.insert(6);
  assertm(result.size() <= 3,
          "Wrong output of close function! Try bigger bounding box!");
  return result;
}

Point BoundingBox::crop_to_box(Point P, const std::set<int> &close_walls) {
  Point new_point = P;
  if (close_walls.find(1) != close_walls.end()) {
    assertm(close_walls.find(2) == close_walls.end(),
            "Wrong output from close_walls function!");
    new_point = Point(min_x(), new_point.y(), new_point.z());
    std::cout << "Close wall 1!" << endl;
  } else if (close_walls.find(2) != close_walls.end()) {
    new_point = Point(max_x(), new_point.y(), new_point.z());
    std::cout << "Close wall 2!" << endl;
  }
  if (close_walls.find(3) != close_walls.end()) {
    assertm(close_walls.find(4) == close_walls.end(),
            "Wrong output from close_walls function!");
    new_point = Point(new_point.x(), min_y(), new_point.z());
    std::cout << "Close wall 3!" << endl;
  } else if (close_walls.find(4) != close_walls.end()) {
    new_point = Point(new_point.x(), max_y(), new_point.z());
    std::cout << "Close wall 4!" << endl;
  }
  if (close_walls.find(5) != close_walls.end()) {
    assertm(close_walls.find(6) == close_walls.end(),
            "Wrong output from close_walls function!");
    new_point = Point(new_point.x(), new_point.y(), min_z());
    std::cout << "Close wall 5!" << endl;
  } else if (close_walls.find(6) != close_walls.end()) {
    new_point = Point(new_point.x(), new_point.y(), max_z());
    std::cout << "Close wall 6!" << endl;
  }
}