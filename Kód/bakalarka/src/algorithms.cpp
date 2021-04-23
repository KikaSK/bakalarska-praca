#include "algorithms.h"
#include "assertm.h"
#include "function.h"

// N-R method for root finding, not necessarily the closest root
numeric Newton_Raphson(const realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point) {

  numeric precision = 1e-6;

  numeric iter = starting_point;
  int iterations = 0;
  while (abs(f.subs(my_x == iter).evalf()) > precision && iterations < 1'000) {
    iterations++;
    assertm(abs(df.subs(my_x == iter).evalf()) > 10e-6,
            "Division by 0 in N-R method!");
    iter -= ex_to<numeric>(f.subs(my_x == iter).evalf() /
                           df.subs(my_x == iter).evalf());
  }
  return iter;
}

// bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter) {
  if (f.subs(my_x == point1).evalf() < 10e-6)
    return point1;
  if (f.subs(my_x == point2).evalf() < 10e-6)
    return point2;
  assertm(iter < 1000, "Too much iterations in bisection method!");
  assertm(f.subs(my_x == point1).evalf() * f.subs(my_x == point2).evalf() <= 0,
          "Wrong call for Bisect function!");
  numeric new_point = (point1 + point2) / 2;
  std::optional<numeric> projected = std::nullopt;
  if (f.subs(my_x == point1).evalf() * f.subs(my_x == new_point).evalf() <= 0) {
    projected = Bisect(my_x, f, point1, new_point, ++iter);
  } else if (f.subs(my_x == point2).evalf() *
                 f.subs(my_x == new_point).evalf() <=
             0) {
    projected = Bisect(my_x, f, new_point, point2, ++iter);
  } else {
    assertm(false, "Wrong value in Bisect!");
  }
  assertm(projected.has_value(), "Projected point withou value!");
  return projected.value();
}

// Bisection is called when N-R proejcts to distant point
// finds 2 points on opposite sides of surface and returns result of bisection
// on these two points
numeric Bisection(const realsymbol my_x, const ex &f, numeric starting_point,
                  numeric e_size) {
  numeric dx = e_size / 10;
  numeric new_point1 = starting_point, last_point1 = starting_point;
  numeric new_point2 = starting_point, last_point2 = starting_point;
  int iterations = 0;
  while ((f.subs(my_x == new_point1).evalf() *
                  f.subs(my_x == last_point1).evalf() >
              0 &&
          f.subs(my_x == new_point2).evalf() *
                  f.subs(my_x == last_point2).evalf() >
              0) &&
         iterations < 1000) {
    iterations++;
    last_point1 = new_point1;
    new_point1 += dx;
    last_point2 = new_point2;
    new_point2 -= dx;
  }

  assertm(iterations < 1000,
          "Bisection method failed, try smaller triangle edge size");

  std::optional<numeric> projected = std::nullopt;

  if (f.subs(my_x == new_point1).evalf() *
          f.subs(my_x == last_point1).evalf() <=
      0)
    projected = Bisect(my_x, f, new_point1, last_point1, 0);
  else if (f.subs(my_x == new_point2).evalf() *
               f.subs(my_x == last_point2).evalf() <=
           0)
    projected = Bisect(my_x, f, new_point2, last_point2, 0);
  else
    assertm(false, "Wrong call in Bisection function!");

  assertm(projected.has_value(), "Projected point without value!");
  return projected.value();
}

// returns projected point in the direction of normal
Point project(Point point_to_project, Vector normal, const Function &F,
              const numeric e_size) {

  realsymbol my_x("my_x");

  numeric starting_point;

  Point P = point_to_project;
  Vector n = normal;

  ex param_x, param_y, param_z;

  // parametric equations of line given by P and n
  // expressing parameter and substituing to other equations

  if (n.x() != 0) {
    starting_point = P.x();
    param_x = my_x;
    param_y = P.y() - n.y() * P.x() / n.x() + (n.y() / n.x()) * my_x;
    param_z = P.z() - n.z() * P.x() / n.x() + (n.z() / n.x()) * my_x;
  } else if (n.y() != 0) {
    starting_point = P.y();
    param_x = P.x();
    param_y = my_x;
    param_z = P.z() - n.z() * P.y() / n.y() + (n.z() / n.y()) * my_x;
  } else if (n.z() != 0) {
    starting_point = P.z();
    param_x = P.x();
    param_y = P.y();
    param_z = my_x;
  } else {
    assertm(false, "Normal is a zero vector!");
  }

  // after this in param_x, param_y and param_z are equations of only one
  // variable my_x substitute to F to get function for Newton-Raphson method

  ex f = F.get_function().subs(
      lst{F.get_x() == param_x, F.get_y() == param_y, F.get_z() == param_z});
  ex df = f.diff(my_x, 1);

  std::optional<Point> projected = std::nullopt;

  numeric root = Newton_Raphson(my_x, f, df, starting_point);
  numeric projected_x = ex_to<numeric>(param_x.subs(my_x == root).evalf());
  numeric projected_y = ex_to<numeric>(param_y.subs(my_x == root).evalf());
  numeric projected_z = ex_to<numeric>(param_z.subs(my_x == root).evalf());

  /*assertm(F.substitute(lst{F.get_x() == projected_x, F.get_y() == projected_y,
                           F.get_z() == projected_z}) < 10e-6,
          "Wrong return from NR method!");*/

  projected = Point(projected_x, projected_y, projected_z);
  if (Vector(point_to_project, projected.value()).get_length() > 4 * e_size) {

    root = Bisection(my_x, f, starting_point, e_size);

    numeric projected_x = ex_to<numeric>(param_x.subs(my_x == root).evalf());
    numeric projected_y = ex_to<numeric>(param_y.subs(my_x == root).evalf());
    numeric projected_z = ex_to<numeric>(param_z.subs(my_x == root).evalf());
    assertm(F.substitute(lst{F.get_x() == projected_x, F.get_y() == projected_y,
                             F.get_z() == projected_z}) < 10e-6,
            "Wrong Bisection calculation!");
  }
  assertm(projected.has_value(), "Not found projected point!");
  assertm(Vector(point_to_project, projected.value()).get_length() < 2 * e_size,
          "Wrong calculation in project function!");
  return projected.value();
}

// connects two vectors of edges
vector<Edge> connect_edges(const vector<Edge> &v1, const vector<Edge> &v2) {
  vector<Edge> connected = v1;
  for (auto edge : v2) {
    connected.push_back(edge);
  }
  return connected;
}

// connects two vectors of points
vector<Point> connect_points(const vector<Point> &v1, const vector<Point> &v2) {
  vector<Point> connected = v1;
  for (auto point : v2) {
    connected.push_back(point);
  }
  return connected;
}

// angle BAP in range (-Pi, Pi) with respect to neighbour triangle
numeric angle(const Edge &working_edge, const Point P, const Triangle &N) {
  if (P == working_edge.A() || P == working_edge.B())
    return false;
  Triangle T = Triangle(working_edge.A(), working_edge.B(), P);
  if (!T.is_triangle())
    return false;
  assertm(N.is_triangle(), "Invalid neighbour triangle!");

  // third point in neighbour triangle

  std::optional<Point> NPoint = std::nullopt;

  if (N.A() != working_edge.A() && N.A() != working_edge.B())
    NPoint = N.A();

  else if (N.B() != working_edge.A() && N.B() != working_edge.B())
    NPoint = N.B();

  else if (N.C() != working_edge.A() && N.C() != working_edge.B())
    NPoint = N.C();

  assertm(NPoint.has_value(), "Not found neighbour point!");

  Vector AB = Vector(working_edge.A(), working_edge.B());
  Vector AP = Vector(working_edge.A(), P);
  Vector AN = Vector(working_edge.A(), NPoint.value());

  numeric cos = (AB.unit() * AP.unit());
  assertm(abs(cos) <= numeric(1), "Wrong value of cosine!");

  numeric angle = acos(cos);

  assertm(angle >= 0, "Angle in the wrong range!");

  Vector ABcrossAN = (AB ^ AN);
  Vector ABcrossAP = (AB ^ AP);

  // if same_direction is > 0 the normals are pointing in aprox same direction
  // thus the points lie on the same side

  numeric same_direction = ABcrossAP * ABcrossAN;

  if (same_direction > 0) {
    angle = -angle;
  }
  assertm(abs(angle) <= ex_to<numeric>(Pi.evalf()),
          "Angle in the wrong range!");

  return angle;
}

// true if angle is between 0 and 9*pi/10 with respect to neighbour triangle
bool good_orientation(const Edge &working_edge, const Point P,
                      const Triangle &N) {
  // Edge reversed_edge = Edge(working_edge.B(), working_edge.A());
  return angle(working_edge, P, N) > ex_to<numeric>(Pi.evalf()) / 10 &&
         angle(working_edge, P, N) <
             numeric(9) * ex_to<numeric>(Pi.evalf()) / 10;
}

// https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d

// returns ditance between point and line segment given by working edge
numeric line_point_dist(const Edge &working_edge, const Point P,
                        const Triangle &neighbour_triangle) {

  Vector AB_unit = Vector(working_edge.A(), working_edge.B()).unit();
  Vector AP = Vector(working_edge.A(), P);
  numeric t = AP * AB_unit;
  Point p = Point(working_edge.A(), t * AB_unit);

  numeric angle1 = angle(working_edge, P, neighbour_triangle);
  numeric angle2 =
      angle(Edge(working_edge.B(), working_edge.A()), P, neighbour_triangle);

  if (abs(angle1) <= ex_to<numeric>(Pi.evalf()) / 2 &&
      abs(angle2) <= ex_to<numeric>(Pi.evalf()) / 2) {
    return Vector(P, p).get_length();
  } else if (abs(angle1) > ex_to<numeric>(Pi.evalf()) / 2) {
    return Vector(working_edge.A(), P).get_length();
  } else {
    return Vector(working_edge.B(), P).get_length();
  }
  assertm(false, "Wrong line point distance!");
  return 1000;
}

// true if edge is active
bool is_active(const Edge &edge, const vector<Edge> &active_edges) {
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
bool is_checked(const Edge &edge, const vector<Edge> &checked_edges) {
  int counter = 0;
  for (auto my_edge : checked_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in active_edges!");
  return (counter == 1);
}

bool is_bounding(const Edge &edge, const BoundingBox &bounding_box) {
  int counter = 0;
  for (auto my_edge : bounding_box.bounding_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in bounding_edges!");
  return (counter == 1);
}

// true if edge is active or checked
bool is_border(const Edge &edge, const vector<Edge> &active_edges,
               const vector<Edge> &checked_edges,
               const BoundingBox &bounding_box) {
  return (is_active(edge, active_edges) || is_checked(edge, checked_edges) ||
          is_bounding(edge, bounding_box));
}

// true if point is on border of mesh
bool is_border_point(Point P, const vector<Edge> &active_edges,
                     const vector<Edge> &checked_edges,
                     const vector<Edge> &bounding_edges) {
  for (auto edge : active_edges) {
    if (edge.A() == P || edge.B() == P)
      return true;
  }

  for (auto edge : checked_edges) {
    if (edge.A() == P || edge.B() == P)
      return true;
  }

  for (auto edge : bounding_edges) {
    if (edge.A() == P || edge.B() == P)
      return true;
  }
  return false;
}

// throws error if it is found more than once
void delete_from_active(const Edge &edge, vector<Edge> &active_edges) {
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
void delete_from_checked(const Edge &edge, vector<Edge> &checked_edges) {
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

// throws error if it is found more than once
void delete_from_bounding(const Edge &edge, BoundingBox &bounding_box) {
  size_t counter = 0;
  std::optional<size_t> index = std::nullopt;
  for (size_t i = 0; i < bounding_box.bounding_edges.size(); ++i) {
    if (bounding_box.bounding_edges[i] == edge) {
      index = i;
      counter++;
    }
  }
  assertm(counter == 1 || counter == 0,
          "More than one edge found while deleting!");
  if (counter == 0)
    return;
  else {
    std::swap(bounding_box.bounding_edges[index.value()],
              bounding_box.bounding_edges.back());
    bounding_box.bounding_edges.pop_back();
  }
  return;
}

// throws error if it is already there
void push_edge_to_active(const Edge &edge, vector<Edge> &active_edges) {

  assertm(!is_active(edge, active_edges), "Edge already in active edges!");
  active_edges.push_back(edge);
  return;
}

// throws error if it is already there
void push_edge_to_checked(const Edge &edge, vector<Edge> &checked_edges) {

  assertm(!is_checked(edge, checked_edges), "Edge already in checked_edges!");
  checked_edges.push_back(edge);
  return;
}

void push_edge_to_bounding(const Edge &edge, BoundingBox &bounding_box) {

  assertm(!is_bounding(edge, bounding_box), "Edge already in bounding_edges!");
  bounding_box.bounding_edges.push_back(edge);
  return;
}

// checks if edges of new triangle are active or are not im mesh
bool good_edges(const Mesh &my_mesh, const vector<Edge> &active_edges,
                const vector<Edge> &checked_edges, const Edge &working_edge,
                const Point &P, const BoundingBox &bounding_box) {
  const Triangle T = Triangle(working_edge.A(), working_edge.B(), P);
  Edge new_edge1(working_edge.A(), P);
  Edge new_edge2(working_edge.B(), P);
  // if(!bounding_box.is_inside(P)) return false;
  bool intersections_off = true;
  return (!((my_mesh.is_in_mesh(new_edge1) &&
             !is_border(new_edge1, active_edges, checked_edges,
                        bounding_box)) ||
            (my_mesh.is_in_mesh(new_edge2) &&
             !is_border(new_edge2, active_edges, checked_edges,
                        bounding_box)))) &&
         (intersections_off || (!my_mesh.intersections(T)));
}

// finds closest border point to edge
std::optional<Point> get_closest_point(const Mesh &my_mesh,
                                       const vector<Edge> &active_edges,
                                       const vector<Edge> &checked_edges,
                                       const Edge &working_edge,
                                       const Triangle &N, const numeric &e_size,
                                       const BoundingBox &bounding_box) {
  // auto border_edges = connect_edges(connect_edges(active_edges,
  // checked_edges), bounding_box.bounding_edges);
  auto border_edges = connect_edges(active_edges, checked_edges);

  std::optional<pair<Point, numeric>> closest_point = std::nullopt;
  for (auto edge : border_edges) {
    if (!closest_point.has_value()) {
      if (good_orientation(working_edge, edge.A(), N) &&
          Triangle(working_edge.A(), working_edge.B(), edge.A())
              .is_triangle() &&
          good_edges(my_mesh, active_edges, checked_edges, working_edge,
                     edge.A(), bounding_box))
        closest_point =
            pair(edge.A(), line_point_dist(working_edge, edge.A(), N));
      if (!closest_point.has_value() &&
          good_orientation(working_edge, edge.B(), N) &&
          Triangle(working_edge.A(), working_edge.B(), edge.B())
              .is_triangle() &&
          good_edges(my_mesh, active_edges, checked_edges, working_edge,
                     edge.B(), bounding_box))
        closest_point =
            pair(edge.B(), line_point_dist(working_edge, edge.B(), N));
    }

    if (closest_point.has_value()) {
      if (line_point_dist(working_edge, edge.A(), N) <
              closest_point.value().second &&
          good_orientation(working_edge, edge.A(), N) &&
          Triangle(working_edge.A(), working_edge.B(), edge.A())
              .is_triangle() &&
          good_edges(my_mesh, active_edges, checked_edges, working_edge,
                     edge.A(), bounding_box))
        closest_point =
            pair(edge.A(), line_point_dist(working_edge, edge.A(), N));
      if (line_point_dist(working_edge, edge.B(), N) <
              closest_point.value().second &&
          good_orientation(working_edge, edge.B(), N) &&
          Triangle(working_edge.A(), working_edge.B(), edge.B())
              .is_triangle() &&
          good_edges(my_mesh, active_edges, checked_edges, working_edge,
                     edge.B(), bounding_box))
        closest_point =
            pair(edge.B(), line_point_dist(working_edge, edge.B(), N));
    }
  }

  if (closest_point.has_value() && closest_point.value().second < 2 * e_size) {
    return closest_point.value().first;
  }
  return std::nullopt;
}


// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(Edge e, const Triangle &T, numeric e_size) {
  assertm(T.is_triangle(), "Getting normal of non-valid triangle!");
  Vector normal = T.get_normal();
  Vector edge_vector(e.A(), e.B());

  Vector direction = (normal ^ edge_vector).unit();

  numeric min_side_length = std::min(
      T.AB().get_length(), std::min(T.BC().get_length(), T.CA().get_length()));
  numeric delta = min_side_length / 20;

  Point P1 = Point(e.get_midpoint(), delta * direction);

  if (T.is_in_triangle(P1)) {
    if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction)))
      direction = numeric(-1) * direction;
    else
      assertm(false, "Both points in triangle!");
  } else
    assertm(T.is_in_triangle(Point(e.get_midpoint(), -delta * direction)),
            "No points in triangle!");

  return direction;
}



