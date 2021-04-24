#include "bounding_box.h"

BoundingBox::BoundingBox(numeric _min_x, numeric _max_x, numeric _min_y,
                         numeric _max_y, numeric _min_z, numeric _max_z)
    : _min_x(_min_x), _max_x(_max_x), _min_y(_min_y),
      _max_y(_max_y), _min_z(_min_z), _max_z(_max_z){};

numeric BoundingBox::min_x() const { return _min_x; }
numeric BoundingBox::max_x() const { return _max_x; }
numeric BoundingBox::min_y() const { return _min_y; }
numeric BoundingBox::max_y() const { return _max_y; }
numeric BoundingBox::min_z() const { return _min_z; }
numeric BoundingBox::max_z() const { return _max_z; }

/*
Vector BoundingBox::normal_x() const { return Vector(1, 0, 0); }
Vector BoundingBox::normal_y() const { return Vector(0, 1, 0); }
Vector BoundingBox::normal_z() const { return Vector(0, 0, 1); }
*/

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

bool BoundingBox::new_bounding_edge(const Edge &e) const{
  return is_on(e.A()) && is_on(e.B());
}

std::optional<Point> BoundingBox::project_on_min_x(const Point &midpoint, const Point &P) const {
  Vector v = Vector(P, midpoint);
  //Vector v = Vector(1, 0, 0);

  numeric dist = Vector(midpoint, P).get_length();

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_x - _min_x;
  vector<ex> d_input = {1, 0, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.x() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) || Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_max_x(const Point &midpoint, const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  //Vector v = Vector(1, 0, 0);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_x - _max_x;
  vector<ex> d_input = {1, 0, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);


  std::optional<Point> projected = std::nullopt;
  if (v.x() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) || Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;

}
std::optional<Point> BoundingBox::project_on_min_y(const Point &midpoint, const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  //Vector v = (0, 1, 0);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_y - _min_y;
  vector<ex> d_input = {0, 1, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.y() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) || Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_max_y(const Point &midpoint, const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  //Vector v = (0, 1, 0);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_y - _max_y;
  vector<ex> d_input = {0, 1, 0};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.y() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) || Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_min_z(const Point &midpoint, const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();

  //Vector v = Vector(0, 0, 1);

  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_z - _min_z;
  vector<ex> d_input = {0, 0, 1};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.z() != 0) {
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) || Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}
std::optional<Point> BoundingBox::project_on_max_z(const Point &midpoint, const Point &P) const {
  Vector v = Vector(P, midpoint);
  numeric dist = Vector(midpoint, P).get_length();
  
  //Vector v = Vector(0, 0, 1);
  
  assertm(v.x() != 0 || v.y() != 0 || v.z() != 0,
          "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");

  ex input = b_z - _max_z;
  vector<ex> d_input = {0, 0, 1};
  Function F = Function(b_x, b_y, b_z, input, d_input);

  std::optional<Point> projected = std::nullopt;
  if (v.z() != 0){
    projected = project(P, v, F, std::nullopt);
    if (!is_on(projected.value()) || Vector(P, projected.value()).get_length() > dist)
      projected = std::nullopt;
  }
  return projected;
}


Point BoundingBox::project_on_box(const Point &midpoint, const Point &P) const {
  
  // projecting point on all of the walls
  vector<std::optional<Point>> points = {std::nullopt, std::nullopt,
                                         std::nullopt, std::nullopt,
                                         std::nullopt, std::nullopt};
  points[0] = project_on_min_x(midpoint, P);
  points[1] = project_on_max_x(midpoint, P);
  points[2] = project_on_min_y(midpoint, P);
  points[3] = project_on_max_y(midpoint, P);
  points[4] = project_on_min_z(midpoint, P);
  points[5] = project_on_max_z(midpoint, P);

  std::optional<numeric> min_dist = std::nullopt;
  std::optional<int> min_dist_index = std::nullopt;
  for (int i = 0; i < 6; ++i) {
    if (points[i].has_value()) {
      numeric dist = Vector(P, points[i].value()).get_length();
      if (!min_dist.has_value()) {
        min_dist = dist;
        min_dist_index = i;
      } else if (dist < min_dist.value()) {
        min_dist = dist;
        min_dist_index = i;
      }
    }
  }
  assertm(min_dist_index.has_value(), "Index without value!");
  assertm(points[min_dist_index.value()].has_value(), "Projected point without value!");
  return points[min_dist_index.value()].value();
}

Point BoundingBox::crop_to_box(const Point &midpoint, const Point &P,
                               const numeric &e_size) const {
  numeric precision = e_size/3;
  Point projected = P;

  if (!is_inside(P)) {
    projected = project_on_box(midpoint, P);
  } else if (abs(P.x() - _min_x) < precision && project_on_min_x(midpoint, P).has_value()) {
    cout<<"Close to min_x"<<endl;
    projected =  project_on_min_x(midpoint, P).value();
  } else if (abs(P.x() - _max_x) < precision && project_on_max_x(midpoint, P).has_value()) {
    cout<<"Close to max_x"<<endl;
    projected =  project_on_max_x(midpoint, P).value();
  } else if (abs(P.y() - _min_y) < precision && project_on_min_y(midpoint, P).has_value()) {
    cout<<"Close to min_y"<<endl;
    projected =  project_on_min_y(midpoint, P).value();
  } else if (abs(P.y() - _max_y) < precision && project_on_max_y(midpoint, P).has_value()) {
    cout<<"Close to max_y"<<endl;
    projected =  project_on_max_y(midpoint, P).value();
  } else if (abs(P.z() - _min_z) < precision && project_on_min_z(midpoint, P).has_value()) {
    cout<<"Close to min_z"<<endl;
    projected =  project_on_min_z(midpoint, P).value();
  } else if (abs(P.z() - _max_z) < precision && project_on_max_z(midpoint, P).has_value()) {
    cout<<"Close to max_z"<<endl;
    projected =  project_on_max_z(midpoint, P).value();
  } else {
    projected = P;
  }
  return projected;
}
