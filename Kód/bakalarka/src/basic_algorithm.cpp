#include "basic_algorithm.h"
#include "algorithms.h"
#include "assertm.h"
#include <iostream>

using std::cout;

void BasicAlgorithm::add_marks() {
  auto border = connect_edges(active_edges, checked_edges);
  for (auto edge : border) {
    auto dir =
        e_size / 5 *
        find_direction(edge, my_mesh.find_triangle_with_edge(edge));
    assertm(!Vector(edge.A(), edge.B()).is_zero(), "Zero edge vector!");
    Vector edge_dir = e_size / 5 * Vector(edge.A(), edge.B()).unit();
    Edge new_e(Point(edge.get_midpoint(), edge_dir.vector_inverse() / 2),
               Point(edge.get_midpoint(), edge_dir / 2));
    auto new_p = Point(edge.get_midpoint(), dir);
    my_mesh.add_triangle(new_e, new_p);
  }
}


// returns vector of points closer than 0.4*e_size to point P sorted from closest to working edge
// if there are no points returns std::nullopt
std::optional<vector<Point>> BasicAlgorithm::
find_close_points(Point P, const Edge &working_edge, const Triangle &neighbour_triangle) const {
  numeric min_dist = 0.4 * e_size;
  vector<Point> close_points;

  vector<Edge>edges = connect_edges(connect_edges(active_edges, checked_edges), bounding_edges);
  for(auto edge : edges){
    if (edge.A() != P && edge.A() != working_edge.A() && edge.A() != working_edge.B()) {
      numeric dist = Vector(edge.A(), P).get_length();
      if(dist<min_dist){
        bool found = false;
        for(auto point : close_points){
          if(point == edge.A()){
            found = true;
          }
        }
        if(!found){
          close_points.push_back(edge.A());
        }
      }
    }
    if (edge.B() != P && edge.B() != working_edge.A() && edge.B() != working_edge.B()) {
      numeric dist = Vector(edge.B(), P).get_length();
      if(dist<min_dist){
        bool found = false;
        for(auto point : close_points){
          if(point == edge.B()){
            found = true;
          }
        }
        if(!found){
          close_points.push_back(edge.B());
        }
      }
    }
  }
  if (close_points.empty())
    return std::nullopt;

  sort(close_points.begin(), close_points.end(),
    [&working_edge, &neighbour_triangle](auto i, auto j) {
              return line_point_dist(working_edge, i, neighbour_triangle) <
                     line_point_dist(working_edge, j, neighbour_triangle);
            });

  return close_points;
}


// fixes corners of bounded triangulation
void BasicAlgorithm::fix_corners() {
  realsymbol my_x("my_x"), my_y("my_y"), my_z("my_z");
  for (auto edge : bounding_edges){
    // binary number, ones at place of side at which point is lying on
    int faces_values_A = bounding_box.faces(edge.A());
    int faces_values_B = bounding_box.faces(edge.B());
    // logical and
    int common_faces = faces_values_A & faces_values_B;
    int index_A;
    int index_B;
    if(common_faces != 0) continue;
      for (int i = 0; i<6; ++i){
        // value of ith bit
        if(faces_values_A & (1<<i)){
          index_A = i;
        }
        if(faces_values_B & (1<<i)){
          index_B = i;
        }
      }
      // direciton vector of intersection line of two faces on which the edge lies
      Vector v(1, 1, 1);
      std::optional<numeric> Px;
      std::optional<numeric> Py;
      std::optional<numeric> Pz;
      if(index_A == 0 || index_B == 0 || index_A == 1 || index_B == 1)
      {
        v = v - Vector(1, 0, 0);
        if(index_A == 0 || index_B == 0){
          Px = bounding_box.min_x();
        }
        else if(index_A == 1 || index_B == 1){
          Px = bounding_box.max_x();
        }
        else{
          assertm(false, "Point of an edge lying on x_min and x_max sides!");
        }
      }
      if(index_A == 2 || index_B == 2 || index_A == 3 || index_B == 3)
      {
        v = v - Vector(0, 1, 0);
        if(index_A == 2 || index_B == 2){
          Py = bounding_box.min_y();
        }
        else if(index_A == 3 || index_B == 3){
          Py = bounding_box.max_y();
        }
        else{
          assertm(false, "Point of an edge lying on y_min and y_max sides!");
        }
      }
      if(index_A == 4 || index_B == 4 || index_A == 5 || index_B == 5)
      {
        v = v - Vector(0, 0, 1);
        if(index_A == 4 || index_B == 4){
          Pz = bounding_box.min_z();
        }
        else if(index_A == 5 || index_B == 5){
          Pz = bounding_box.max_z();
        }
        else{
          assertm(false, "Point of an edge lying on z_min and z_max sides!");
        }
      }
      if(!Px.has_value()) Px = edge.get_midpoint().x();
      if(!Py.has_value()) Py = edge.get_midpoint().y();
      if(!Pz.has_value()) Pz = edge.get_midpoint().z();

      std::optional<Point> projected = std::nullopt;
      Point P(Px.value(), Py.value(), Pz.value());
      if(v.is_zero()){
        projected = P;
      }
      else{
        projected = project(P, v, F, e_size);
      }
      assertm(projected.has_value(), "Point without value!");
      if(Vector(projected.value(), edge.get_midpoint()).get_length() < edge.get_length()){
        if(Triangle(edge.A(), edge.B(), projected.value()).is_triangle()/* && good_orientation(edge, projected, my_mesh.find_triangle_with_edge(edge))*/){
          create_triangle(edge, projected.value());
        }
      }
  }
  return;
}

// checks if conditions required in the first part of the algorithm are
// satisfied
bool BasicAlgorithm::Delaunay_conditions(const Edge &working_edge,
                                         const Point &P,
                                         const Triangle &neighbour_triangle) const {
  if(working_edge.A() == P || working_edge.B() == P) return false;
  Triangle T = Triangle(working_edge.A(), working_edge.B(), P);
  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);
  // all the conditions for a new triangle in the first part of algorithm
  return (T.is_triangle() &&
          good_orientation(working_edge, P, neighbour_triangle) &&
          my_mesh.check_Delaunay(T, working_edge, neighbour_triangle) && good_edges(working_edge, P));
}

// checks if conditions required in the first part of the algorithm are
// satisfied
bool BasicAlgorithm::non_Delaunay_conditions(const Edge &working_edge,
                                         const Point &P,
                                         const Triangle &neighbour_triangle) const{
  if(working_edge.A() == P || working_edge.B() == P) return false;
  Triangle T = Triangle(working_edge.A(), working_edge.B(), P);
  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);
  // all the conditions for a new triangle in the second part of algorithm
  return (T.is_triangle() &&
          good_orientation(working_edge, P, neighbour_triangle) &&
          good_edges(working_edge, P));
}

// creates new triangle and adds it to mesh
void BasicAlgorithm::create_triangle(const Edge &working_edge, const Point &P) {

  Edge new_edge1(working_edge.A(), P);
  Edge new_edge2(working_edge.B(), P);

  assertm(new_edge1 != new_edge2, "Same edges!");

  my_mesh.add_triangle(working_edge, P);
  update_border(new_edge1, new_edge2);
  return;
}

// updates active and checked edges and returns number of new edges
int BasicAlgorithm::update_border(const Edge &new_edge1,
                                  const Edge &new_edge2) {

  assertm(new_edge1 != new_edge2, "Same edges while updating border!");
  int new_edges = 0;

  if (is_border(new_edge1)) {
    delete_from_active(new_edge1);
    delete_from_checked(new_edge1);
  } else {
    push_edge_to_active(new_edge1);
    new_edges++;
  }

  if (is_border(new_edge2)) {
    delete_from_active(new_edge2);
    delete_from_checked(new_edge2);
  } else {
    push_edge_to_active(new_edge2);
    new_edges++;
  }
  
  if(bounding_box.new_bounding_edge(new_edge1)){
    delete_from_active(new_edge1);
    delete_from_checked(new_edge1);
    if(!is_bounding(new_edge1))
      bounding_edges.push_back(new_edge1);
  }
  if(bounding_box.new_bounding_edge(new_edge2)){
    delete_from_active(new_edge2);
    delete_from_checked(new_edge2);
    if(!is_bounding(new_edge2))
      bounding_edges.push_back(new_edge2);
  }

  return new_edges;
}

// makes triangle if prev == next
// TODO basic triangle na okraji!!
bool BasicAlgorithm::basic_triangle(const Edge &working_edge,
                                    const Triangle &neighbour_triangle,
                                    const Point &prev, const Point &next) {

  // determine other point on the neighbour triangle
  Point neighbour_triangle_point = neighbour_triangle.A();
  if (neighbour_triangle.A() != working_edge.A() &&
      neighbour_triangle.A() != working_edge.B())
    neighbour_triangle_point = neighbour_triangle.A();
  else if (neighbour_triangle.B() != working_edge.A() &&
           neighbour_triangle.B() != working_edge.B())
    neighbour_triangle_point = neighbour_triangle.B();
  else
    neighbour_triangle_point = neighbour_triangle.C();

  // if prev and next are the and its not only one triangle same create triangle
  if (prev == next && neighbour_triangle_point != prev && Triangle(working_edge.A(), working_edge.B(), prev).is_triangle()) {
    create_triangle(working_edge, prev);
    return true;
  }
  return false;
}

// true if edge is active
bool BasicAlgorithm::is_active(const Edge &edge) const{
  int counter = 0;
  for (auto my_edge : active_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in active_edges!");
  return (counter == 1);
}

// true if edge is checked
bool BasicAlgorithm::is_checked(const Edge &edge) const{
  int counter = 0;
  for (auto my_edge : checked_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in active_edges!");
  return (counter == 1);
}

bool BasicAlgorithm::is_bounding(const Edge &edge) const {
  int counter = 0;
  for (auto my_edge : bounding_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in bounding_edges!");
  return (counter == 1);
}

// true if edge is active or checked
bool BasicAlgorithm::is_border(const Edge &edge) const {
  return (is_active(edge) || is_checked(edge) || is_bounding(edge));
}

// true if point is on border of mesh
bool BasicAlgorithm::is_border_point(Point P) const {
  for (auto edge : active_edges) {
    if (edge.A() == P || edge.B() == P)
    {
      return true;
    }
  }

  for (auto edge : checked_edges) {
    if (edge.A() == P || edge.B() == P)
    {
      return true;
    }
  }

  for (auto edge : bounding_edges) {
    if (edge.A() == P || edge.B() == P)
    {
      return true;
    }
  }
  
  return false;
}

// throws error if it is found more than once
void BasicAlgorithm::delete_from_active(const Edge &edge) {
  size_t counter = 0;
  std::optional<size_t> index = std::nullopt;
  for (size_t i = 0; i < active_edges.size(); ++i) {
    if (active_edges[i] == edge) {
      index = i;
      counter++;
    }
  }
  assertm(counter == 1 || counter == 0,
          "More than one edge found while deleting!");
  if (counter == 0)
    return;
  else {
    std::swap(active_edges[index.value()], active_edges.back());
    active_edges.pop_back();
  }
  return;
}

// throws error if it is found more than once
void BasicAlgorithm::delete_from_checked(const Edge &edge) {
  int counter = 0;
  std::optional<size_t> index = std::nullopt;
  for (size_t i = 0; i < checked_edges.size(); ++i) {
    if (checked_edges[i] == edge) {
      index = i;
      counter++;
    }
  }
  assertm(counter == 1 || counter == 0,
          "More than one edge found while deleting!");
  if (counter == 0)
    return;
  else {
    std::swap(checked_edges[index.value()], checked_edges.back());
    checked_edges.pop_back();
  }
  return;
}

// throws error if it is already there
void BasicAlgorithm::push_edge_to_active(const Edge &edge) {

  assertm(!is_active(edge), "Edge already in active edges!");
  active_edges.push_back(edge);
  return;
}

// throws error if it is already there
void BasicAlgorithm::push_edge_to_checked(const Edge &edge) {

  assertm(!is_checked(edge), "Edge already in checked_edges!");
  checked_edges.push_back(edge);
  return;
}

// finds closest border point to edge
std::optional<Point> BasicAlgorithm::get_closest_point(const Edge &working_edge,
                                                       const Triangle &neighbour_triangle) const {
  auto border_edges = connect_edges(active_edges, checked_edges);
  std::optional<Point> closest_point = std::nullopt;
  std::optional<numeric> min_dist;
  for (auto edge : border_edges) {
      if (!closest_point.has_value()) {
        if (edge.A()!= working_edge.A() && edge.A() != working_edge.B() &&
        non_Delaunay_conditions(working_edge, edge.A(), neighbour_triangle))
        {
          closest_point = edge.A();
          min_dist = line_point_dist(working_edge, closest_point.value(), neighbour_triangle);
        }
        if (!closest_point.has_value()
        && edge.B()!= working_edge.A() && edge.B() != working_edge.B()
        && non_Delaunay_conditions(working_edge, edge.B(), neighbour_triangle))
        {
          closest_point = edge.B();
          min_dist = line_point_dist(working_edge, closest_point.value(), neighbour_triangle);
        }
      }
      if (closest_point.has_value()) {
        if (line_point_dist(working_edge, edge.A(), neighbour_triangle) < min_dist 
                && edge.A()!= working_edge.A() && edge.A() != working_edge.B() && 
                non_Delaunay_conditions(working_edge, edge.A(), neighbour_triangle))
                {
                  closest_point = edge.A();
                  min_dist = line_point_dist(working_edge, closest_point.value(), neighbour_triangle);
                }
        if (line_point_dist(working_edge, edge.B(), neighbour_triangle) < min_dist 
            && edge.B()!= working_edge.A() && edge.B() != working_edge.B() &&
            non_Delaunay_conditions(working_edge, edge.B(), neighbour_triangle))
            {
              closest_point = edge.B();
              min_dist = line_point_dist(working_edge, closest_point.value(), neighbour_triangle);
            }
      }
    }
    

  if (closest_point.has_value() && min_dist.has_value() && min_dist.value() < 1.5*e_size) {
    return closest_point.value();
  }
  return std::nullopt;
}

// checks if edges of new triangle are active or are not im mesh
bool BasicAlgorithm::good_edges(const Edge &working_edge, const Point &P) const {
  
  Edge new_edge1(working_edge.A(), P);
  Edge new_edge2(working_edge.B(), P);

  return !((my_mesh.is_in_mesh(new_edge1) && !is_border(new_edge1)) ||
           (my_mesh.is_in_mesh(new_edge2) && !is_border(new_edge2)));
}

// finds closest border edge to point P
std::optional<pair<Edge, numeric>>
BasicAlgorithm::get_closest_edge(const Point &P, const Triangle &N) const {
  auto border = connect_edges(active_edges, checked_edges);
  std::optional<pair<Edge, numeric>> closest_edge = std::nullopt;
  numeric dist = 0;
  for (auto edge : border) {
    dist = line_point_dist(edge, P, N);
    if (!closest_edge.has_value()) {
      closest_edge = pair(edge, dist);
    } else if (dist < closest_edge.value().second) {
      closest_edge = pair(edge, dist);
    }
  }
  assertm(closest_edge.has_value(), "Edge without value!");
  return closest_edge; // TODO: return non-optional
}

// finds neighbour of prev/next which has the smallest angle with the working
// edge
pair<std::optional<Point>, std::optional<Point>>
BasicAlgorithm::find_closest_prev_next(const Edge &working_edge,
                                       const Triangle &neighbour_triangle,
                                       const vector<Point> &prev,
                                       const vector<Point> &next) const {

  std::optional<numeric> min_prev_angle = std::nullopt;
  std::optional<Point> min_prev_point = std::nullopt;

  std::optional<numeric> my_angle = std::nullopt;

  for (auto prev_point : prev) {
    // we will use angle function to find smallest angle
    my_angle = angle(working_edge, prev_point, neighbour_triangle);
    assertm(my_angle.has_value(), "Angle without value!");
    if (my_angle.value() < 0)
      my_angle = my_angle.value() + ex_to<numeric>(2 * Pi.evalf());
    assertm(my_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                my_angle.value() >= 0,
            "Wrong angle interval!");

    if (min_prev_angle.has_value()) {
      if (my_angle.value() < min_prev_angle.value()) {
        min_prev_angle = my_angle.value();
        min_prev_point = prev_point;
      }
    } else {
      min_prev_angle = my_angle.value();
      min_prev_point = prev_point;
    }
    assertm(min_prev_angle.has_value(), "Angle without value!");
    assertm(min_prev_point.has_value(), "Point without value!");
    assertm(min_prev_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                min_prev_angle.value() >= 0,
            "Angle not in right interval!");
  }

  std::optional<numeric> min_next_angle = std::nullopt;
  std::optional<Point> min_next_point = std::nullopt;

  my_angle = std::nullopt;

  for (auto next_point : next) {
    my_angle = angle(Edge(working_edge.B(), working_edge.A()), next_point,
                     neighbour_triangle);
    assertm(my_angle.has_value(), "Angle without value!");

    if (my_angle.value() < 0)
      my_angle = my_angle.value() + ex_to<numeric>(2 * Pi.evalf());
    assertm(my_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                my_angle.value() >= 0,
            "Wrong angle interval!");

    if (min_next_angle.has_value()) {
      if (my_angle.value() < min_next_angle.value()) {
        min_next_angle = my_angle.value();
        min_next_point = next_point;
      }
    } else {
      min_next_angle = my_angle.value();
      min_next_point = next_point;
    }

    assertm(min_next_angle.has_value(), "Angle without value!");
    assertm(min_next_point.has_value(), "Point without value!");
    assertm(min_next_angle.value() < ex_to<numeric>(2 * Pi.evalf()) &&
                min_next_angle.value() >= 0,
            "Angle not in right interval!");
  }

  return pair(min_prev_point, min_next_point);
}

pair<Point, Point>
BasicAlgorithm::find_prev_next(const Edge &working_edge,
                               const Triangle &neighbour_triangle) const {
  vector<Point> prev;
  vector<Point> next;

  vector<Edge> border_edges = connect_edges(connect_edges(active_edges, checked_edges), bounding_edges);

  for (auto curr_edge : border_edges) {
    if (curr_edge.A() == working_edge.A()) {
      prev.push_back(curr_edge.B());
    }
    if (curr_edge.B() == working_edge.A()) {
      prev.push_back(curr_edge.A());
    }
    if (curr_edge.A() == working_edge.B()) {
      next.push_back(curr_edge.B());
    }
    if (curr_edge.B() == working_edge.B()) {
      next.push_back(curr_edge.A());
    }
  }

  assertm(!prev.empty() && !next.empty(),
          "Neighbour edge not found in border edges!");

  for (auto prev_point : prev) {
    assertm(prev_point != working_edge.B(),
            "Working edge found in border edges!");
  }
  for (auto next_point : next) {
    assertm(next_point != working_edge.A(),
            "Working edge found in border edges!");
  }

  if (prev.size() == 1 && next.size() == 1) {
    return pair(prev[0], next[0]);
  }

  auto [closest_prev, closest_next] =
      find_closest_prev_next(working_edge, neighbour_triangle, prev, next);

  assertm(closest_prev.has_value() && closest_next.has_value(),
          "Closest prev or next vithout value!");

  return pair(closest_prev.value(), closest_next.value());
}

bool BasicAlgorithm::overlap_normals_check(const Point candidate,
                                           const Edge &working_edge) const {
  /*
  if (candidate == prev || candidate == next || candidate == working_edge.A() ||
      candidate == working_edge.B())
    return false;
    */
   if (candidate == working_edge.A() ||
      candidate == working_edge.B())
    return false;

  Triangle my_triangle(working_edge.A(), working_edge.B(), candidate);

  if (my_triangle.is_triangle() && good_edges(working_edge, candidate)) {
    Vector my_normal = F.outside_normal(my_triangle, e_size);

    for (auto edge : active_edges) {

      // vertex found in active_edges
      if (edge.A() == candidate || edge.B() == candidate) {

        // normal of overlap triangle
        Triangle overlap_triangle = my_mesh.find_triangle_with_edge(edge);
        Vector overlap_normal = F.outside_normal(overlap_triangle, e_size);

        // if normals have the same orientation
        if (overlap_normal * my_normal > 0) {
          return true;
        }
        return false;
      }
    }

    for (auto edge : checked_edges) {
      // vertex found in checked_edges
      if (edge.A() == candidate || edge.B() == candidate) {

        // normal of overlap triangle
        Triangle overlap_triangle = my_mesh.find_triangle_with_edge(edge);
        Vector overlap_normal = F.outside_normal(overlap_triangle, e_size);

        // if normals have the same orientation
        if (overlap_normal * my_normal > 0) {
          return true;
        }
        return false;
      }
    }
  }

  return false;
}

Point BasicAlgorithm::get_projected(const Edge &working_edge,
                                    const Triangle &neighbour_triangle) const {
  Point center = working_edge.get_midpoint();
  assertm(Vector(working_edge.A(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6 &&
              Vector(working_edge.B(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6,
          "Wrong get_midpoint function!");

  
  auto [neighbour1, neighbour2] =
      find_prev_next(working_edge, neighbour_triangle);
  
  // height of equilateral triangle based on e_size
   numeric basic_height = e_size * sqrt(numeric(3)) / 2;
  
  numeric average =
      (1 / numeric(3)) * (Edge(working_edge.A(), neighbour1).get_length() +
                          Edge(working_edge.A(), neighbour1).get_length() +
                          working_edge.get_length());

  // non adaptive height
  numeric height = basic_height;

  // height of equilateral triangle based on working_edge size
  // numeric height = working_edge.get_length() * sqrt(numeric(3)) / 2;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;
  
  // height of equilateral triangle based on neighbour edges size
  // numeric height = average * sqrt(numeric(3)) / 2;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;
  
  // height of equilateral triangle based on neighbour edges size with influence of e_size
  // numeric height = 0.75*(average * sqrt(numeric(3)) / 2) + 0.25*basic_height;
  // if(height<e_size/3) height = 0.25*(e_size/3) + 0.75*height;
  // if(height>2*e_size) height = 2*e_size;
  

  Vector direction =
      height * find_direction(working_edge, neighbour_triangle);
  assertm(direction * neighbour_triangle.get_normal() < 10e-10,
          "Wrong direction!");
  assertm(direction * Vector(working_edge.A(), working_edge.B()) < 10e-10,
          "Wrong direction!");

  Point P(center, direction);
  assertm(Vector(center, P).get_length() - height < 10e-10,
          "Wrong point to project!");

  /*
  Vector n_A = F.get_gradient_at_point(working_edge.A()).unit();
  Vector n_B = F.get_gradient_at_point(working_edge.B()).unit();

  Vector normal = (n_A + n_B) / 2;
  */
  assertm(!F.get_gradient_at_point(P).is_zero(), "Zero gradient!");
  Vector normal = F.get_gradient_at_point(P).unit();
  Point projected = project(P, normal, F, {e_size});

  /*
  assertm((Vector(working_edge.A(), projected).get_length() < 4 * e_size) &&
              (Vector(working_edge.B(), projected).get_length() < 4 * e_size),
          "Projected point too far!");
  */
  return projected;
}

// returns true if prev/next triangle is added to mesh else returns false
bool BasicAlgorithm::fix_prev_next(const Edge &working_edge,
                                   const Triangle &neighbour_triangle,
                                   const bool is_prev, const bool Delaunay) {

  // find previous and next points
  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);
  
  /*
  assertm(Vector(working_edge.A(), prev).get_length() < 3 * e_size,
          "Wrong prev point!");
  assertm(Vector(working_edge.B(), next).get_length() < 3 * e_size,
          "Wrong prev point!");
  */
  assertm(is_border(Edge(working_edge.A(), prev)), "Prev edge not in border!");
  assertm(is_border(Edge(working_edge.B(), next)), "Next edge not in border!");

  Point vertex = prev;
  Edge edge = working_edge;
  // fix so that vertex is the new vertex and edge.A() is the adjacent vertex
  if (!is_prev) {
    edge = Edge(working_edge.B(), working_edge.A());
    vertex = next;
  }

  // potentialy new triangle
  Triangle maybe_new_T(edge.A(), edge.B(), vertex);

  if (maybe_new_T.AB().get_length() > 2 * e_size ||
      maybe_new_T.CA().get_length() > 2 * e_size ||
      maybe_new_T.BC().get_length() > 2 * e_size) {
    return false;
  }

  // checks if the potential triangle has good orientation and angle near A is
  // less than 90 degrees and checks Delaunay
  if ((Delaunay*Delaunay_conditions(edge, vertex, neighbour_triangle)
  || !Delaunay*non_Delaunay_conditions(edge, vertex, neighbour_triangle))) {

    create_triangle(edge, vertex);
    return true;
  }
  return false;
}

// returns true if overlap triangle is added to mesh else returns false
bool BasicAlgorithm::fix_overlap(const Edge &working_edge,
                                 const Triangle &neighbour_triangle,
                                 Point overlap_point, const bool Delaunay) {

  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);

  // checks if overlap point is not neighbour or working edge point also if
  // overlap is on the border and if overlap triangle has good orientation of
  // normal
  if (overlap_normals_check(overlap_point, working_edge))
  {
    Triangle maybe_new_T =
        Triangle(working_edge.A(), working_edge.B(), overlap_point);

    // if Delaunay constraint is satisfied add the triangle to
    // triangulation and end
    if (Delaunay*Delaunay_conditions(working_edge, overlap_point, neighbour_triangle) 
    || !Delaunay*non_Delaunay_conditions(working_edge, overlap_point, neighbour_triangle)) {

      /*
      assertm(Vector(working_edge.A(), overlap_point).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of overlap point!");
      assertm(Vector(working_edge.B(), overlap_point).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of overlap point!");
      */
      create_triangle(working_edge, overlap_point);

      return true;
    }
  }
  return false;
}

bool BasicAlgorithm::fix_proj(const Edge &working_edge, const Point &projected,
                              const Triangle &neighbour_triangle, const Point &prev, const Point &next) {

  Triangle maybe_new_T(working_edge.A(), working_edge.B(), projected);
  assertm(maybe_new_T.is_triangle(), "Projected triangle not valid!");

  // checks if there are some points very close to projected point
  // as surrounding points were taken working_edge points
  if (auto surrounding_points = find_close_points(
          projected, working_edge, neighbour_triangle);
      surrounding_points.has_value()) {

    // points closer to projected point than 0.4*e_size sorted from closest
    vector<Point> close_points = surrounding_points.value();
    for (auto close_point : close_points) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), close_point);
      if (Delaunay_conditions(working_edge, close_point, neighbour_triangle)) {

        // if close point is prev we want to try fix prev
        if (close_point == prev) {
          if (fix_prev_next(working_edge, neighbour_triangle, true, true)) {
            return true;
          }
        }
        // if close point is next we want to try fix next
        else if (close_point == next) {
          if (fix_prev_next(working_edge, neighbour_triangle, false, true)) {
            return true;
          }
        }

        // if close point is overlap we want to try fix overlap
        else {
          if (fix_overlap(working_edge, neighbour_triangle, close_point, true)) {
            return true;
          }
        }
      }
    }
    // if there are close points but nothing worked we want to try to
    // construct original triangle
  }
  auto close_edge = get_closest_edge(projected, neighbour_triangle).value();
  assertm(is_border(close_edge.first), "Closest edge not in border edges!");
  if (close_edge.second < e_size / 3) {
    Edge closest_edge = close_edge.first;
    Point P1 = closest_edge.A();
    Point P2 = closest_edge.B();

    /*
    Vector n_A = F.get_gradient_at_point(closest_edge.A()).unit();
    Vector n_B = F.get_gradient_at_point(closest_edge.B()).unit();
    Vector normal = (n_A + n_B) / 2;
    */

    assertm(!F.get_gradient_at_point(closest_edge.get_midpoint()).is_zero(), "Zero gradient!");
    Vector normal = F.get_gradient_at_point(closest_edge.get_midpoint()).unit();
    
    
    if (P1 != working_edge.A() && P1 != working_edge.B()) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P1);
      if (Delaunay_conditions(working_edge, P1, neighbour_triangle)) {
        create_triangle(working_edge, P1);
        { return true; }
      }
    }

    if (P2 != working_edge.A() && P2 != working_edge.B()) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P2);
      if (Delaunay_conditions(working_edge, P2, neighbour_triangle)) {
        create_triangle(working_edge, P2);
        { return true; }
      }
    }
    // TODO opraviť divide_triangle_by_point
    
    // point in the middle of side
    Point P3 = project(closest_edge.get_midpoint(), normal, F, {e_size});
    //Point P3 = closest_edge.get_midpoint();
    if (P3 != working_edge.A() && P3 != working_edge.B()) {

      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P3);
      if (Delaunay_conditions(working_edge, P3, neighbour_triangle)) {

        // point P3 is not a vertex, so we need to subdivide triagle with this
        // edge and add subdivided triangles to mesh
        
        Edge new_edge1(working_edge.A(), P3);
        Edge new_edge2(working_edge.B(), P3);

        delete_from_active(closest_edge);
        delete_from_checked(closest_edge); 
        my_mesh.divide_triangle_by_point(closest_edge, P3);
        //update_border(new_edge1, new_edge2);
        //if(working_edge.A() != closest_edge.A() && working_edge.B() != closest_edge.A())
          push_edge_to_active(Edge(closest_edge.A(), P3));
        //if(working_edge.A() != closest_edge.B() && working_edge.B() != closest_edge.B())
          push_edge_to_active(Edge(closest_edge.B(), P3));

        
        my_mesh.add_triangle(working_edge, P3);

        update_border(new_edge1, new_edge2);
        return true;
      }
    }


  }
  if (Delaunay_conditions(working_edge, projected, neighbour_triangle)
   /* my_mesh.check_Delaunay(maybe_new_T) &&
      good_edges(working_edge, projected)*/) {
    Point clipped = bounding_box.crop_to_box(working_edge.get_midpoint(), projected, e_size, F);
    if(non_Delaunay_conditions(working_edge, clipped, neighbour_triangle))
    {
      create_triangle(working_edge, clipped);
      return true;
    }
  }
  // if nothing worked
  return false;
}

bool BasicAlgorithm::fix_breakers(const Edge &working_edge, const Point &projected, const Triangle &neighbour_triangle, const bool Delaunay){
  
  Triangle proj_T(working_edge.A(), working_edge.B(), projected);
  
  // points that break Delaunay constraint
  vector<Point> breakers =
      my_mesh.get_breakers(proj_T, active_edges, checked_edges);

  // sort from closest to working_edge
  std::sort(breakers.begin(), breakers.end(),
            [&working_edge, &neighbour_triangle](auto i, auto j) {
              return line_point_dist(working_edge, i, neighbour_triangle) <
                     line_point_dist(working_edge, j, neighbour_triangle);
            });

  // try create triangle with breakers
  for (auto point : breakers) {
    if (is_border_point(point) && fix_overlap(working_edge, neighbour_triangle, point, Delaunay)) {
        return true;
    }
  }
  return false;
}
bool BasicAlgorithm::fix_same_points(const Edge &working_edge, const Point &projected, const Triangle &neighbour_triangle){
  if(is_border_point(projected) && Delaunay_conditions(working_edge, projected, neighbour_triangle)){
    create_triangle(working_edge, projected);
    return true;
  }
  return false;
}

// one step of the algorithm
bool BasicAlgorithm::step(const Edge &working_edge) {

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);
  // find candidate point for working_edge
  Point projected = get_projected(working_edge, neighbour_triangle);

  assertm(!is_border(working_edge), "Working edge found in border edges!");

  // if point projects on already existing border point
  if(fix_same_points(working_edge, projected, neighbour_triangle)){
    return true;
  }

  // if there is a hole in triangle shape, fill it with triangle
  if (basic_triangle(working_edge, neighbour_triangle, prev, next)) {
    return true;
  }

  if(fix_breakers(working_edge, projected, neighbour_triangle, true)){
    return true;
  }

  if (fix_proj(working_edge, projected, neighbour_triangle, prev, next)) {
    return true;
  }
  // tries to add triangle with prev point, true for prev
  if (fix_prev_next(working_edge, neighbour_triangle, true, true)) {
    return true;
  }
  // tries to add triangle with next point, false for next
  if (fix_prev_next(working_edge, neighbour_triangle, false, true)) {
    return true;
  }

  assertm(!is_border(working_edge), "Working edge found in border edges!");
  push_edge_to_checked(working_edge);
  return false;
}

bool BasicAlgorithm::fix_holes2(const Edge &working_edge){
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  auto [prev, next] = find_prev_next(working_edge, neighbour_triangle);

  // if there is a hole in triangle shape, fill it with triangle
  if (basic_triangle(working_edge, neighbour_triangle, prev, next)) {
    return true;
  }
  // find candidate point for working_edge
  Point projected = get_projected(working_edge, neighbour_triangle);
  std::optional<Point> closest_point = get_closest_point(working_edge, neighbour_triangle);
  if(closest_point.has_value()){
    create_triangle(working_edge, closest_point.value());
    return true;
  }

  if(fix_breakers(working_edge, projected, neighbour_triangle, false)){
    return true;
  }
  if(fix_prev_next(working_edge, neighbour_triangle, true, false)){
    return true;
  }
  if(fix_prev_next(working_edge, neighbour_triangle, false, false)){
    return true;
  }

  push_edge_to_checked(working_edge);
  return false;
}

void BasicAlgorithm::starting() {
  int round = 0;
  while (!active_edges.empty()) {
    round++;
    std::optional<Edge> working_edge = std::nullopt;

    // std::random_shuffle(active_edges.begin(), active_edges.end());
    reverse(active_edges.begin(), active_edges.end());
    working_edge = active_edges.back();
    active_edges.pop_back();
    reverse(active_edges.begin(), active_edges.end());
    // step returns true if new triangle is created
    if (step(working_edge.value()))
      assertm(!is_border(working_edge.value()),
              "Working edge found in border!");
    else
      assertm(is_checked(working_edge.value()), "Checked edge not checked");

    // once in every 50 steps prints number of triangles and number of active
    // edges and updates output file
    if (round % 50 == 0) {
      my_mesh.obj_format(name);
      my_mesh.cout_triangles_number();
      cout << "Number of active edges: " << active_edges.size() << endl << endl;
    }

    if(round % 50 == 0){
      //my_mesh.obj_format(name);
    }
  }

    // output
    my_mesh.obj_format(name);

  // this means there are no holes so we finished
  if (checked_edges.empty()) {
    cout << "No holes!" << endl;
    return;
  }
  return;
}

void BasicAlgorithm::ending() {
  assertm(active_edges.empty(), "Called ending with non-empty active edges!");

  // else push checked edges to active, clear checked and call ending
  
  int round0 = 0;
  while (!checked_edges.empty() && round0 < 2){
    round0++;
    active_edges = checked_edges;
    checked_edges.clear();
  
    int round = 0;
    cout<<"Fixing holes!"<<endl;
    while (!active_edges.empty()) {

      round++;
      std::optional<Edge> working_edge = std::nullopt;

      working_edge = active_edges.back();
      active_edges.pop_back();

      assertm(working_edge.has_value(), "No working edge!");

      Triangle neighbour_triangle =
          my_mesh.find_triangle_with_edge(working_edge.value());
      fix_holes2(working_edge.value());
      if (round % 50 == 0) {
        my_mesh.cout_triangles_number();
        cout << "Number of active edges: " << active_edges.size() << endl << endl;
      }
    }
  }
  if(!checked_edges.empty()) 
  {
    cout<<"Some holes stayed!"<<endl;
    for (auto edge : checked_edges)
    {
      cout<<edge<<endl;
    }
  }
  return;
}

Mesh BasicAlgorithm::calculate() {

  starting();
  ending();
  fix_corners();
  //my_mesh.measure(bounding_edges, F, name, e_size);
  //my_mesh.adaptive(0.005, F, e_size);

  my_mesh.obj_format(name);
  return my_mesh;
}
