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


Point BoundingBox::project_on_box(const Edge &working_edge, const Point &P){
  Vector v = Vector(P, working_edge.get_midpoint());

  assertm(v.x() != 0 && v.y() != 0 && v.z() != 0, "Projecting in direction of zero vector!");

  realsymbol b_x("b_x"), b_y("b_y"), b_z("b_z");
  
  //implicit functions of all 6 bounding_box walls
  ex input_min_x = b_x - _min_x;
  vector<ex>d_input_min_x = {1, 0, 0};
  Function F_min_x = Function(b_x, b_y, b_z, input_min_x, d_input_min_x);

  ex input_max_x = b_x - _max_x;
  vector<ex>d_input_max_x = {1, 0, 0};
  Function F_max_x = Function(b_x, b_y, b_z, input_max_x, d_input_max_x);

  ex input_min_y = b_y - _min_y;
  vector<ex>d_input_min_y = {0, 1, 0};
  Function F_min_y = Function(b_x, b_y, b_z, input_min_y, d_input_min_y);

  ex input_max_y = b_y - _max_y;
  vector<ex>d_input_max_y = {0, 1, 0};
  Function F_max_y = Function(b_x, b_y, b_z, input_max_y, d_input_max_y);

  ex input_min_z = b_z - _min_z;
  vector<ex>d_input_min_z = {0, 0, 1};
  Function F_min_z = Function(b_x, b_y, b_z, input_min_z, d_input_min_z);

  ex input_max_z = b_z - _max_z;
  vector<ex>d_input_max_z = {0, 0, 1};
  Function F_max_z = Function(b_x, b_y, b_z, input_max_z, d_input_max_z);

  // projecting point on all of the walls
  vector<std::optional<Point>> points = {std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt, std::nullopt};
  if(v.x() != 0){
    points[0] = project(P, v, F_min_x, std::nullopt);
    points[1] = project(P, v, F_max_x, std::nullopt);
    if(!is_on(points[0].value())) points[0] = std::nullopt;
    if(!is_on(points[1].value())) points[1] = std::nullopt;
  }
  if(v.y() != 0){
    points[2] = project(P, v, F_min_y, std::nullopt);
    points[3] = project(P, v, F_max_y, std::nullopt);
    if(!is_on(points[2].value())) points[2] = std::nullopt;
    if(!is_on(points[3].value())) points[3] = std::nullopt;
  }
  if(v.z() != 0){
    points[4] = project(P, v, F_min_z, std::nullopt);
    points[5] = project(P, v, F_max_z, std::nullopt);
    if(!is_on(points[4].value())) points[4] = std::nullopt;
    if(!is_on(points[5].value())) points[5] = std::nullopt;
  }
  
  std::optional<numeric> min_dist = std::nullopt;
  std::optional<int> min_dist_index = std::nullopt;
  for (int i = 0; i<6; ++i){
    if(points[i].has_value()){
      numeric dist = Vector(P, points[i].value()).get_length();
      if(!min_dist.has_value()){
        min_dist = dist;
        min_dist_index = i;
      }
      else if(dist < min_dist.value()){
        min_dist = dist;
        min_dist_index = i;
      }
    }
  }
  assertm(min_dist_index.has_value(), "Index without value!");
  return points[min_dist_index.value()].value();

}


Point BoundingBox::crop_to_box(const Edge &working_edge, const Point& P, const numeric &e_size) {
  numeric precision = e_size/4;
  Point projected = P;
  if(!is_inside(P)){
    projected = project_on_box(working_edge, P);
  } 
}


/*
Point BoundingBox::crop_to_box(Point &P, const Vector &v) 
{ 
  std::optional<numeric> tmin = std::nullopt, tmax = std::nullopt;
  if(v.x() != 0)
  {
    tmin = (_min_x - P.x()) / v.x(); 
    tmax = (_max_x - P.x()) / v.x(); 
      
    if (tmin.value() > tmax.value()) swap(tmin, tmax);
  }

  std::optional<numeric> tymin = std::nullopt, tymax = std::nullopt;
  if(v.y() != 0)
  {
    tymin = (_min_y - P.y()) / v.y(); 
    tymax = (_max_y - P.y()) / v.y(); 
 
    if (tymin.value() > tymax.value()) swap(tymin, tymax);
  }
  
  if(tmin.has_value() && tymin.has_value()){
    if ((tmin.value() > tymax.value()) || (tymin.value() > tmax.value())) 
        //return false; 
 
    if (tymin.value() > tmin.value()) 
        tmin = tymin.value(); 
    
    if (tymax.value() < tmax.value()) 
        tmax = tymax.value(); 
  }
   
  std::optional<numeric> tzmin = std::nullopt, tzmax = std::nullopt;
    
  if(v.z() != 0)
  {
    tzmin = (_min_z - P.z()) / v.z(); 
    tzmax = (_max_z - P.z()) / v.z(); 
 
    if (tzmin.value() > tzmax.value()) swap(tzmin, tzmax); 
  }
  
  if(tmin.has_value() && tzmin.has_value()){
    if ((tmin.value() > tzmax.value()) || (tzmin.value() > tmax.value())) 
        //return false; 
 
    if (tzmin.value() > tmin.value()) 
        tmin = tzmin.value(); 
 
    if (tzmax.value() < tmax.value()) 
        tmax = tzmax.value(); 
 
    return true;
} 

*/

