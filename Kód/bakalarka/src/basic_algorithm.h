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

#include "algorithms.h"
#include "assertm.h"
#include <iostream>

using std::cout;

Point project(Point point_to_project, Vector normal, const Function &F,
              const numeric e_size) {

  realsymbol my_x("my_x");

  numeric starting_point;

  // for better orientation
  Point P = point_to_project;
  // std::cout<<"Point to project: " << point_to_project<<endl;
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

  assertm(F.substitute(lst{F.get_x() == projected_x, F.get_y() == projected_y,
                           F.get_z() == projected_z}) < 10e-6,
          "Wrong return from NR method!");

  projected = Point(projected_x, projected_y, projected_z);
  // cout<< "Projected point afetr NR: " << projected.value() <<endl;
  // cout<< "Distance from point to project: " << Vector(point_to_project,
  // projected.value()).get_length() << endl; cout<< "Edge size is: " << e_size
  // <<endl;
  if (Vector(point_to_project, projected.value()).get_length() > 4 * e_size) {

    // cout<<"Starting point in bisection: " << starting_point << endl;
    root = Bisection(my_x, f, starting_point, e_size);
    // cout<< "Fount root of bisection:" << root << endl;
    numeric projected_x = ex_to<numeric>(param_x.subs(my_x == root).evalf());
    numeric projected_y = ex_to<numeric>(param_y.subs(my_x == root).evalf());
    numeric projected_z = ex_to<numeric>(param_z.subs(my_x == root).evalf());
    assertm(F.substitute(lst{F.get_x() == projected_x, F.get_y() == projected_y,
                             F.get_z() == projected_z}) < 10e-6,
            "Wrong Bisection calculation!");
  }
  assertm(projected.has_value(), "Not found projected point!");
  // cout<< "Projected point afetr Bisection: " << projected.value() <<endl;
  // cout<< "Distance from point to project: " << Vector(point_to_project,
  // projected.value()).get_length() << endl; cout<< "Edge size is: " << e_size
  // <<endl;
  assertm(Vector(point_to_project, projected.value()).get_length() < 4 * e_size,
          "Wrong calculation in project function!");
  return projected.value();
}

// a
// b
// c

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(Edge e, Triangle &T, numeric e_size) {
  assertm(T.is_triangle(), "Getting normal of non-valid triangle!");
  Vector normal = T.get_normal();
  Vector edge_vector(e.A(), e.B());

  Vector direction = (normal ^ edge_vector).unit();

  // cout<< "Inside function, my directon: " << direction <<endl;
  numeric min_side_length = std::min(
      T.AB().get_length(), std::min(T.BC().get_length(), T.CA().get_length()));
  numeric delta = min_side_length / 20;

  Point P1 = Point(e.get_midpoint(), delta * direction);

  if (T.is_in_triangle(P1)) {
    if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction))) {
      direction = numeric(-1) * direction;
    } else {
      /*cout<<"Triangle: " << T << endl;
      cout<<"Point1: " << P1 << endl;
      cout<<"Point2: " << Point(e.get_midpoint(), -delta*direction) << endl;
      */
      assertm(false, "Both points in triangle!");
    }
  } else {
    if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction)))
      cout << "Would die!" << endl;
    /*if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction))) {
      assertm(false, "No points in triangle!");
    }*/
  }

  return direction;
}

// a
// b
// c

bool is_active(const Edge &edge, vector<Edge> &active_edges) {
  int counter = 0;
  for (auto my_edge : active_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in active_edges!");
  return (counter == 1);
}

bool is_checked(const Edge &edge, vector<Edge> &checked_edges) {
  int counter = 0;
  for (auto my_edge : checked_edges) {
    if (my_edge == edge)
      counter++;
  }
  assertm(counter == 0 || counter == 1,
          "More than one same edges in active_edges!");
  return (counter == 1);
}

bool is_border(const Edge &edge, vector<Edge> &active_edges,
               vector<Edge> &checked_edges) {
  return (is_active(edge, active_edges) || is_checked(edge, checked_edges));
}

bool is_border_point(Point P, vector<Edge> &active_edges,
                     vector<Edge> &checked_edges) {
  for (auto edge : active_edges) {
    if (edge.A() == P || edge.B() == P)
      return true;
  }

  for (auto edge : checked_edges) {
    if (edge.A() == P || edge.B() == P)
      return true;
  }
  return false;
}

// a
// b
// c

void delete_from_active(const Edge &edge, vector<Edge> &active_edges) {
  size_t counter = 0;
  std::optional<size_t> index = std::nullopt;
  for (size_t i = 0; i < active_edges.size(); ++i) {
    if (active_edges[i] == edge) {
      index = i;
      counter++;
    }
  }
  // cout<< counter<<endl;
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

// a
// b
// c

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

// a
// b
// c

void push_edge_to_active(const Edge &edge, vector<Edge> &active_edges) {

  assertm(!is_active(edge, active_edges), "Edge already in active edges!");
  active_edges.push_back(edge);
  return;
}

// a
// b
// c

void push_edge_to_checked(const Edge &edge, vector<Edge> &checked_edges) {

  assertm(!is_checked(edge, checked_edges), "Edge already in checked_edges!");
  checked_edges.push_back(edge);
  return;
}

// a
// b
// c

// angle BAP in range (-Pi, Pi)
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
  // std::cout<< "NPoint: " << NPoint.value()<< endl;

  Vector AB = Vector(working_edge.A(), working_edge.B());
  Vector AP = Vector(working_edge.A(), P);
  Vector AN = Vector(working_edge.A(), NPoint.value());

  numeric cos = (AB.unit() * AP.unit());
  assertm(abs(cos) <= numeric(1), "Wrong value of cosine!");

  numeric angle = acos(cos);

  assertm(angle >= 0, "Angle in the wrong range!");

  Vector ABcrossAN = (AB ^ AN);
  Vector ABcrossAP = (AB ^ AP);

  // std::cout << "AB cross AN: " << ABcrossAN << endl;
  // std::cout << "AB cross AP: " << ABcrossAP << endl;

  // if same_direction is > 0 the normals are pointing in aprox same direction
  // thus the points lie on the same side

  numeric same_direction = ABcrossAP * ABcrossAN;

  // std::cout << "Dir: " << same_direction << endl;

  if (same_direction > 0) {
    angle = -angle;
  }
  assertm(abs(angle) <= ex_to<numeric>(Pi.evalf()),
          "Angle in the wrong range!");

  return angle;
}

bool good_orientation(const Edge &working_edge, const Point P,
                      const Triangle &N) {
  return angle(working_edge, P, N) > 0 &&
         angle(working_edge, P, N) < 3 * ex_to<numeric>(Pi.evalf()) / 4;
}

// https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

// returns ditance between point and line given by working edge
numeric line_point_dist(const Edge &working_edge, const Point P,
                        const Triangle &neighbour_triangle) {
  Vector AP(working_edge.A(), P);
  Vector BP(working_edge.B(), P);

  numeric angle1 = angle(working_edge, P, neighbour_triangle);
  numeric angle2 =
      angle(Edge(working_edge.B(), working_edge.A()), P, neighbour_triangle);

  if (abs(angle1) < ex_to<numeric>(Pi.evalf()) &&
      abs(angle2) < ex_to<numeric>(Pi.evalf())) {
    return (AP ^ BP).get_length() / working_edge.get_length();
  } else if (abs(angle1) > ex_to<numeric>(Pi.evalf())) {
    return Vector(working_edge.A(), P).get_length();
  } else {
    return Vector(working_edge.B(), P).get_length();
  }
  assertm(false, "Wrong line point distance!");
  return 1000;
}

// finds point breaking Delaunay closest to line given by working edge
std::optional<Point> closest_point(const Edge &working_edge,
                                   const vector<Point> &breakers,
                                   const Triangle &neighbour_triangle) {
  std::optional<numeric> min_dist = std::nullopt;
  std::optional<Point> min_dist_point = std::nullopt;

  std::optional<numeric> my_dist = std::nullopt;

  for (auto point : breakers) {
    if (point == working_edge.A() || point == working_edge.B())
      continue;
    my_dist = line_point_dist(working_edge, point, neighbour_triangle);
    assertm(my_dist.has_value(), "Point without value!");

    if (!min_dist.has_value()) {
      min_dist = my_dist.value();
      min_dist_point = point;
    } else if (my_dist.value() < min_dist.value()) {
      min_dist = my_dist.value();
      min_dist_point = point;
    }
  }

  return min_dist_point;
}

// a
// b
// c

// finds neighbour of prev/next which has the smallest angle with the working
// edge
pair<std::optional<Point>, std::optional<Point>>
find_closest_prev_next(const Mesh &my_mesh, const Edge &working_edge,
                       const vector<Point> &prev, const vector<Point> &next) {

  std::optional<numeric> min_prev_angle = std::nullopt;
  std::optional<Point> min_prev_point = std::nullopt;

  const Triangle neighbour_triangle =
      my_mesh.find_triangle_with_edge(working_edge);

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

pair<Point, Point> find_prev_next(const Mesh &my_mesh, const Edge &working_edge,
                                  const vector<Edge> &active_edges,
                                  const vector<Edge> &checked_edges) {
  vector<Point> prev;
  vector<Point> next;

  int counter = 0;

  for (auto curr_edge : active_edges) {
    if (curr_edge.A() == working_edge.A()) {
      prev.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == working_edge.A()) {
      prev.push_back(curr_edge.A());
      ++counter;
    }
    if (curr_edge.A() == working_edge.B()) {
      next.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == working_edge.B()) {
      next.push_back(curr_edge.A());
      ++counter;
    }
  }

  for (auto curr_edge : checked_edges) {
    if (curr_edge.A() == working_edge.A()) {
      prev.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == working_edge.A()) {
      prev.push_back(curr_edge.A());
      ++counter;
    }
    if (curr_edge.A() == working_edge.B()) {
      next.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == working_edge.B()) {
      next.push_back(curr_edge.A());
      ++counter;
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
  } else {
    std::optional<Point> my_prev = std::nullopt;
    std::optional<Point> my_next = std::nullopt;

    auto [closest_prev, closest_next] =
        find_closest_prev_next(my_mesh, working_edge, prev, next);

    assertm(closest_prev.has_value() && closest_next.has_value(),
            "Closest prev or next vithout value!");
    my_prev = closest_prev.value();
    my_next = closest_next.value();
    /*if (closest_prev.has_value()) {
      my_prev = closest_prev.value();
    } else {
      // it's going to be checked and will fail
      my_prev = prev[0];
    }

    if (closest_next.has_value()) {
      my_next = closest_next.value();
    } else {
      // it's going to be checked and will fail
      my_next = next[0];
    }*/

    assertm(my_prev.has_value() && my_next.has_value(),
            "Prev or next point without value!");

    return pair(my_prev.value(), my_next.value());
  }
}

// a
// b
// c

bool is_vertex_good_possibility(const Point candidate, const Point prev,
                                const Point next, const Edge &working_edge,
                                const Triangle &neighbour_triangle,
                                const vector<Edge> &active_edges,
                                const vector<Edge> &checked_edges,
                                const Mesh &my_mesh, const Function &F) {
  if (candidate == prev || candidate == next || candidate == working_edge.A() ||
      candidate == working_edge.B())
    return false;

  Triangle my_triangle(working_edge.A(), working_edge.B(), candidate);

  if (my_triangle.is_triangle()) {
    Vector my_normal = F.outside_normal(my_triangle);

    for (auto edge : active_edges) {

      // vertex found in active_edges
      if (edge.A() == candidate || edge.B() == candidate) {

        // normal of overlap triangle
        Triangle overlap_triangle = my_mesh.find_triangle_with_edge(edge);
        Vector overlap_normal = F.outside_normal(overlap_triangle);

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
        Vector overlap_normal = F.outside_normal(overlap_triangle);

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

// a
// b
// c

Point get_projected(const Edge &working_edge, const numeric e_size,
                    const Mesh &my_mesh, const Function F) {
  Point center = working_edge.get_midpoint();
  assertm(Vector(working_edge.A(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6 &&
              Vector(working_edge.B(), center).get_length() -
                      working_edge.get_length() / 2 <
                  10e-6,
          "Wrong get_midpoint function!");
  // cout<<"Center: " << center<<endl;
  // numeric height = working_edge.get_length() * sqrt(numeric(3)) / 2;
  numeric height = e_size * sqrt(numeric(3)) / 2;
  // assertm(abs(height - e_size * sqrt(numeric(3)) / 2) < e_size / 3,
  //        "Weird size of height!");

  // TODO: checkovat ci je triangle iba jeden
  cout<<"OK1"<<endl;
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  cout<<"OK2"<<endl;
  assertm(neighbour_triangle.is_triangle(), "Neighbour triangle not valid!");
  // cout<<"Neighbour Triangle: " << endl << neighbour_triangle <<endl;
  Vector direction =
      height * find_direction(working_edge, neighbour_triangle, e_size);
  // cout<<"Direction: " << direction << endl;
  assertm(direction * neighbour_triangle.get_normal() < 10e-8,
          "Wrong direction!");
  assertm(direction * Vector(working_edge.A(), working_edge.B()) < 10e-8,
          "Wrong direction!");
  // TODO: kontrolovat ci je P v rovine neighbour triangle
  Point P(center, direction);
  assertm(Vector(center, P).get_length() - height < 10e-8,
          "Wrong point to project!");
  // TODO: checknut gradient
  Vector n_A = F.get_gradient_at_point(working_edge.A()).unit();
  Vector n_B = F.get_gradient_at_point(working_edge.B()).unit();

  Vector normal = (n_A + n_B) / 2;
  Point projected = project(P, normal, F, e_size);

  assertm((Vector(working_edge.A(), projected).get_length() < 4 * e_size) &&
              (Vector(working_edge.B(), projected).get_length() < 4 * e_size),
          "Projected point too far!");

  return projected;
}

// a
// b
// c

// returns true if prev/next triangle is added to mesh else returns false
bool fix_prev_next(Mesh &my_mesh, vector<Edge> &active_edges,
                   vector<Edge> &checked_edges, const Edge &working_edge,
                   const bool is_prev, const numeric e_size) {

  // find previous and next points
  auto [prev, next] =
      find_prev_next(my_mesh, working_edge, active_edges, checked_edges);
  assertm(Vector(working_edge.A(), prev).get_length() < 3 * e_size,
          "Wrong prev point!");
  assertm(Vector(working_edge.B(), next).get_length() < 3 * e_size,
          "Wrong prev point!");
  assertm(is_border(Edge(working_edge.A(), prev), active_edges, checked_edges),
          "Prev edge not in border!");
  assertm(is_border(Edge(working_edge.B(), next), active_edges, checked_edges),
          "Next edge not in border!");

  Point vertex = prev;
  Edge edge = working_edge;

  // fix so that vertex is the new vertex and edge.A() is the adjacent vertex
  if (!is_prev) {
    edge = Edge(working_edge.B(), working_edge.A());
    vertex = next;
  }

  // potentialy new triangle
  Triangle maybe_new_T(edge.A(), edge.B(), vertex);

  // already existing neighbour triangle of working edge
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(edge);

  // checks if the potential triangle has good orientation and angle near A is
  // less than 90 degrees and checks Delaunay
  if (maybe_new_T.is_triangle() &&
      good_orientation(edge, vertex, neighbour_triangle) &&
      my_mesh.check_Delaunay(maybe_new_T)) {

    // cout << "Found prev triangle!" << endl;

    assertm(Vector(working_edge.A(), vertex).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev or next point!");
    assertm(Vector(working_edge.B(), vertex).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev or next point!");

    my_mesh.add_triangle(edge, vertex);

    assertm(is_border(Edge(edge.A(), vertex), active_edges, checked_edges),
            "Neighbour edge not in border!");
    assertm(!(is_border(Edge(edge.B(), vertex), active_edges, checked_edges)) ||
                (prev == next),
            "New edges already in border!");

    // delete_from_active(Edge(edge.B(), vertex), active_edges);
    // delete_from_checked(Edge(edge.B(), vertex), checked_edges);

    delete_from_active(Edge(edge.A(), vertex), active_edges);
    delete_from_checked(Edge(edge.A(), vertex), checked_edges);

    if (prev == next) {
      delete_from_active(Edge(edge.B(), vertex), active_edges);
      delete_from_checked(Edge(edge.B(), vertex), checked_edges);
    } else
      push_edge_to_active(Edge(edge.B(), vertex), active_edges);
    return true;
  }
  return false;
}

// a
// b
// c

// returns true if overlap triangle is added to mesh else returns false
bool fix_overlap(Mesh &my_mesh, const Edge &working_edge,
                 vector<Edge> &active_edges, vector<Edge> &checked_edges,
                 Point overlap_point, const Function &F) {
  // assertm(false, "In check overlap!");

  auto [prev, next] =
      find_prev_next(my_mesh, working_edge, active_edges, checked_edges);

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  assertm(neighbour_triangle.is_triangle(), "Neighbour triangle not valid!");

  // checks if overlap point is not neighbour or working edge point also if
  // overlap is on the border and if overlap triangle has good orientation of
  // normal
  if (is_vertex_good_possibility(overlap_point, prev, next, working_edge,
                                 neighbour_triangle, active_edges,
                                 checked_edges, my_mesh, F))

  {
    Triangle maybe_new_T =
        Triangle(working_edge.A(), working_edge.B(), overlap_point);

    // if Delaunay constraint is satisfied add the triangle to
    // triangulation and end
    if (maybe_new_T.is_triangle() &&
        good_orientation(working_edge, overlap_point, neighbour_triangle) &&
        my_mesh.check_Delaunay(maybe_new_T)) {

      // assertm(false, "Found overlap triangle, hooray!");
      // cout << "Found overlap triangle" << endl;

      assertm(Vector(working_edge.A(), overlap_point).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of overlap point!");
      assertm(Vector(working_edge.B(), overlap_point).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of overlap point!");
      my_mesh.add_triangle(working_edge, overlap_point);

      Edge new_edge1 = Edge(working_edge.A(), overlap_point);
      Edge new_edge2 = Edge(working_edge.B(), overlap_point);

      assertm(new_edge1 != new_edge2, "Same edges!");

      if (!is_border(new_edge1, active_edges, checked_edges)) {
        push_edge_to_active(new_edge1, active_edges);
      } else {
        delete_from_active(new_edge1, active_edges);
        delete_from_checked(new_edge1, checked_edges);
      }
      if (!is_border(new_edge2, active_edges, checked_edges)) {
        push_edge_to_active(new_edge2, active_edges);
      } else {
        delete_from_active(new_edge2, active_edges);
        delete_from_checked(new_edge2, checked_edges);
      }
      return true;
    }
  }
  return false;
}

// a
// b
// c

bool fix_proj(Mesh &my_mesh, vector<Edge> &active_edges,
              vector<Edge> &checked_edges, const Edge &working_edge,
              numeric e_size, const Function &F) {

  Point projected = get_projected(working_edge, e_size, my_mesh, F);

  Triangle maybe_new_T(working_edge.A(), working_edge.B(), projected);
  assertm(maybe_new_T.is_triangle(), "Proj triangle not valid!");

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  // checks if there are some points very close to projected point
  if (my_mesh.empty_surrounding(projected, e_size, active_edges, checked_edges)
          .has_value()) {

    // points closer to projected point than 0.3*e_size sorted from closest
    // TODO sort not working!!!
    vector<Point> close_points =
        my_mesh
            .empty_surrounding(projected, e_size, active_edges, checked_edges)
            .value();

    for (auto close_point : close_points) {

      Triangle maybe_new_T(working_edge.A(), working_edge.B(), close_point);

      if (maybe_new_T.is_triangle() &&
          good_orientation(working_edge, close_point, neighbour_triangle)) {

        if (my_mesh.check_Delaunay(maybe_new_T)) {

          auto [prev, next] = find_prev_next(my_mesh, working_edge,
                                             active_edges, checked_edges);

          // if close point is prev we want to try fix prev
          if (close_point == prev) {
            if (fix_prev_next(my_mesh, active_edges, checked_edges,
                              working_edge, true, e_size))
              return true;
          }
          // if close point is next we want to try fix next
          else if (close_point == next) {
            if (fix_prev_next(my_mesh, active_edges, checked_edges,
                              working_edge, false, e_size))
              return true;
          }

          // if close point is overlap we want to try fix overlap
          else {
            if (fix_overlap(my_mesh, working_edge, active_edges, checked_edges,
                            close_point, F))
              return true;
          }
        }
      }
    }
    // if there are close points but nothing worked we want to try to
    // construct original triangle
    // cout << "There are close points but nothing worked!" << endl;
  }

  else if (my_mesh.check_Delaunay(maybe_new_T)) {

    // cout << "Found new triangle" << endl;

    assertm(Vector(working_edge.A(), projected).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of projected point!");
    assertm(Vector(working_edge.B(), projected).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of projected point!");

    my_mesh.add_triangle(working_edge, projected);

    /*assertm(!is_border(Edge(working_edge.A(), projected), active_edges,
                       checked_edges) &&
                !is_border(Edge(working_edge.B(), projected), active_edges,
                           checked_edges),
            "New edges already in border!");*/
    assertm(Edge(working_edge.A(), projected) !=
                Edge(working_edge.B(), projected),
            "Same edges!");
    if (is_border(Edge(working_edge.A(), projected), active_edges,
                  checked_edges)) {
      delete_from_active(Edge(working_edge.A(), projected), active_edges);
    } else {
      push_edge_to_active(Edge(working_edge.A(), projected), active_edges);
    }
    if (is_border(Edge(working_edge.B(), projected), active_edges,
                  checked_edges)) {
      delete_from_active(Edge(working_edge.B(), projected), active_edges);
    } else {
      push_edge_to_active(Edge(working_edge.B(), projected), active_edges);
    }
    return true;
  }
  // if nothing worked we move on to another step
  return false;
}

// a
// b
// c

// one step of the algorithm
void step(Mesh &my_mesh, vector<Edge> &active_edges,
          vector<Edge> &checked_edges, const Edge &working_edge,
          const numeric e_size, const Function &F) {

  my_mesh.obj_format();

  assertm(!is_border(working_edge, active_edges, checked_edges),
          "Workind edge found in border edges!");

  if (fix_proj(my_mesh, active_edges, checked_edges, working_edge, e_size, F)) {
    return;
  } else if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                           true, e_size))
    return;
  else if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                         false, e_size))
    return;
  else {
    Point projected = get_projected(working_edge, e_size, my_mesh, F);
    Triangle proj_T(working_edge.A(), working_edge.B(), projected);
    vector<Point> breakers = my_mesh.get_breakers(proj_T, active_edges, checked_edges);

    for (auto point : breakers) {
      Triangle neighbour_T = my_mesh.find_triangle_with_edge(working_edge);
      if (good_orientation(working_edge, point, neighbour_T)) {
        if (is_border_point(point, active_edges, checked_edges)) {
          assertm(Vector(point, working_edge.A()).get_length() < 3 * e_size,
                  "Too big distance of break point!");
          if (fix_overlap(my_mesh, working_edge, active_edges, checked_edges,
                          point, F))
            return;
        }
      }
    }

    // cout << "No triangle found" << endl;
    // if no triangle was created
    assertm(!is_border(working_edge, active_edges, checked_edges),
            "Something very wrong!");
    // something is wtong with find_prev_next
    /*
        auto [prev, next] =
            find_prev_next(my_mesh, working_edge, active_edges, checked_edges);

        Triangle neighbour_T = my_mesh.find_triangle_with_edge(working_edge);

        if (prev == next) {
          Triangle new_T(working_edge.A(), working_edge.B(), prev);
          if (new_T.is_triangle() && !(neighbour_T == new_T) &&
              good_orientation(working_edge, prev, neighbour_T)) {
            my_mesh.add_triangle(working_edge, prev);
            delete_from_active(Edge(working_edge.A(), prev), active_edges);
            delete_from_active(Edge(working_edge.B(), prev), active_edges);
            delete_from_checked(Edge(working_edge.A(), prev), checked_edges);
            delete_from_checked(Edge(working_edge.B(), prev), checked_edges);
            return;
          }
        } else {
          std::optional<Point> min_dist_point = closest_point(working_edge,
       breakers); if(!min_dist_point.has_value()){ Triangle
       new_prev_T(working_edge.A(), working_edge.B(), prev); Triangle
       new_next_T(working_edge.A(), working_edge.B(), next);
            if(new_prev_T.is_triangle() && !(neighbour_T == new_prev_T) &&
              good_orientation(working_edge, prev, neighbour_T)){
                delete_from_active(Edge(working_edge.A(), prev), active_edges);
                delete_from_checked(Edge(working_edge.A(), prev),
       checked_edges); push_edge_to_active(Edge(working_edge.B(), prev),
       active_edges);
              }
            else if(new_next_T.is_triangle() && !(neighbour_T == new_next_T) &&
              good_orientation(working_edge, next, neighbour_T)){
                delete_from_active(Edge(working_edge.B(), next), active_edges);
                delete_from_checked(Edge(working_edge.B(), next),
       checked_edges); push_edge_to_active(Edge(working_edge.A(), next),
       active_edges);
            }
            else{
              push_edge_to_checked(working_edge, checked_edges);
            }
            return;
          }
          assertm(min_dist_point.value() != working_edge.A() &&
                      min_dist_point.value() != working_edge.B(),
                  "Closest_point function is wrong!");

          Triangle closest_point_T(working_edge.A(), working_edge.B(),
                                   min_dist_point.value());
          if (closest_point_T.is_triangle() &&
              good_orientation(working_edge, min_dist_point.value(),
       neighbour_T)) { my_mesh.add_triangle(working_edge,
       min_dist_point.value());

            Edge new_edge1(working_edge.A(), min_dist_point.value());
            Edge new_edge2(working_edge.B(), min_dist_point.value());

            assertm(new_edge1 != new_edge2, "Same edges!");
            assertm(prev != next, "Prev should not be next!");

            if(is_border(new_edge1, active_edges, checked_edges)){
              delete_from_active(new_edge1, active_edges);
              delete_from_checked(new_edge1, checked_edges);
            }
            else{
              push_edge_to_active(new_edge1, active_edges);
            }

            if(is_border(new_edge2, active_edges, checked_edges)){
              delete_from_active(new_edge2, active_edges);
              delete_from_checked(new_edge2, checked_edges);
            }
            else{
              push_edge_to_active(new_edge2, active_edges);
            }
            return;
          }
          else {*/
    push_edge_to_checked(working_edge, checked_edges);
    //}
    //}
  }
  // push_edge_to_checked(working_edge, checked_edges);
  return;
}

int fix_holes(Mesh &my_mesh, const Function &F, const Edge &working_edge,
               vector<Edge> &active_edges, vector<Edge> &checked_edges,
               const numeric e_size) {
  
  int number_of_new_edges = 0;

  cout << "In fix_holes!" << endl;

  Point projected = get_projected(working_edge, e_size, my_mesh, F);

  Triangle maybe_new_T(working_edge.A(), working_edge.B(), projected);
  assertm(maybe_new_T.is_triangle(), "Proj triangle not valid!");

  cout<<"HERE"<<endl;
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  cout<<"HERE2"<<endl;
  auto breakers = my_mesh.get_breakers(maybe_new_T, active_edges, checked_edges);
  auto close_points =
      my_mesh.empty_surrounding(projected, e_size, active_edges, checked_edges);
  if (close_points.has_value()) {
    if (breakers.empty()) {
      breakers = close_points.value();
    }
    breakers.insert(breakers.end(), close_points.value().begin(),
                    close_points.value().end());
  }
  if (breakers.empty()) {
    if (fix_proj(my_mesh, active_edges, checked_edges, working_edge, e_size, F))
      return true;
    else
      push_edge_to_checked(working_edge, checked_edges);
    return false;
  }
  assertm(!breakers.empty(), "No close points!");

  std::optional<pair<Point, numeric>> closest_point = std::nullopt;

  for (auto point : breakers) {
    if (!closest_point.has_value()) {
      closest_point =
          pair(point, line_point_dist(working_edge, point, neighbour_triangle));
    } else if (Triangle(working_edge.A(), working_edge.B(), point)
                   .is_triangle() &&
               line_point_dist(working_edge, point, neighbour_triangle) <
                   closest_point.value().second &&
               good_orientation(working_edge, point, neighbour_triangle)) {
      closest_point =
          pair(point, line_point_dist(working_edge, point, neighbour_triangle));
    }
  }
  if (closest_point.has_value()) {
    Triangle new_T(working_edge.A(), working_edge.B(),
                   closest_point.value().first);
    my_mesh.add_triangle(working_edge, closest_point.value().first);
    my_mesh.obj_format();
    cout << "New triangle1!" << endl;
    Edge new_edge1(working_edge.A(), closest_point.value().first);
    Edge new_edge2(working_edge.B(), closest_point.value().first);
    bool new_edge = false;
    if (is_border(new_edge1, active_edges, checked_edges)) {
      delete_from_checked(new_edge1, checked_edges);
      delete_from_active(new_edge1, active_edges);
      cout << "deleted edge" << endl;
    } else {
      push_edge_to_active(new_edge1, active_edges);
      number_of_new_edges++;
    }
    if (is_border(new_edge2, active_edges, checked_edges)) {
      delete_from_checked(new_edge2, checked_edges);
      delete_from_active(new_edge2, active_edges);
    } else {
      push_edge_to_active(new_edge2, active_edges);
      number_of_new_edges++;
    }
    return number_of_new_edges;
  }

  auto [prev, next] =
      find_prev_next(my_mesh, working_edge, active_edges, checked_edges);
  Triangle T = Triangle(working_edge.A(), working_edge.B(), prev);
  Triangle NeighbourT = my_mesh.find_triangle_with_edge(working_edge);
  if (T.is_triangle() && good_orientation(working_edge, prev, NeighbourT)) {
    assertm(Vector(working_edge.A(), prev).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev point!");
    assertm(Vector(working_edge.B(), prev).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev point!");

    my_mesh.add_triangle(working_edge, prev);
    my_mesh.obj_format();
    cout << "Find new triangle2!" << endl;
    delete_from_active(Edge(working_edge.A(), prev), active_edges);
    delete_from_checked(Edge(working_edge.A(), prev), checked_edges);
    if (is_border(Edge(working_edge.B(), prev), active_edges, checked_edges)) {
      delete_from_active(Edge(working_edge.B(), prev), active_edges);
      delete_from_checked(Edge(working_edge.B(), prev), checked_edges);
    } else {
      push_edge_to_active(Edge(working_edge.B(), prev), active_edges);
      number_of_new_edges++;
    }
  } else if (!is_border(Edge(working_edge.B(), prev), active_edges,
                        checked_edges)) {
    T = Triangle(working_edge.A(), working_edge.B(), next);
    if (T.is_triangle() && good_orientation(working_edge, next, NeighbourT)) {

      assertm(Vector(working_edge.A(), next).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of next point!");
      assertm(Vector(working_edge.B(), next).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of next point!");

      my_mesh.add_triangle(working_edge, next);
      my_mesh.obj_format();
      cout << "New triangle3!" << endl;
      delete_from_active(Edge(working_edge.B(), next), active_edges);
      if (is_border(Edge(working_edge.A(), next), active_edges,
                    checked_edges)) {
        delete_from_active(Edge(working_edge.A(), next), active_edges);
        delete_from_checked(Edge(working_edge.A(), next), checked_edges);
      } else {
        push_edge_to_active(Edge(working_edge.A(), next), active_edges);
        number_of_new_edges++;
      }
    }
  } else {
    push_edge_to_checked(working_edge, checked_edges);
  }
  return number_of_new_edges;
}

void ending(Mesh my_mesh, vector<Edge> & active_edges, vector<Edge> & checked_edges, const Function & F, numeric e_size);

void starting(Mesh my_mesh, vector<Edge> & active_edges, vector<Edge> & checked_edges, const Function & F, numeric e_size)
{
  int round = 0;
  while (!active_edges.empty()) {
    round++;
    std::optional<Edge> working_edge = std::nullopt;

    // std::random_shuffle(active_edges.begin(), active_edges.end());

    working_edge = active_edges.back();

    assertm(working_edge.has_value(), "No working edge!");

    delete_from_active(working_edge.value(), active_edges);

    // active_edges.pop_back();

    // cout << "Current working edge: " << endl << working_edge.value() <<
    // endl;

    step(my_mesh, active_edges, checked_edges, working_edge.value(), e_size, F);
    // my_mesh.cout_triangles();
    if (round % 100 == 0) {
      my_mesh.cout_triangles_number();
      cout << "Number of edges in active_edges: " << active_edges.size()
           << endl;
      cout << endl;
    }
    my_mesh.obj_format();
  }
  
  return;
}
void ending(Mesh my_mesh, vector<Edge> & active_edges, vector<Edge> & checked_edges, const Function & F, numeric e_size){
  int round = 0;
  while (!active_edges.empty()) {

    round++;
    std::optional<Edge> working_edge = std::nullopt;

    // std::random_shuffle(active_edges.begin(), active_edges.end());

    working_edge = active_edges.back();

    assertm(working_edge.has_value(), "No working edge!");

    // delete_from_active(working_edge.value(), active_edges);

    active_edges.pop_back();

    // cout << "Current working edge: " << endl << working_edge.value() <<
    // endl;

    //fix_holes(my_mesh, F, working_edge.value(), active_edges, checked_edges,
    //          e_size);
    int new_edges = fix_holes(my_mesh,  F, working_edge.value(), active_edges,
    checked_edges, e_size);
    assertm(new_edges == 0 || new_edges == 1 || new_edges == 2, "Wrong number of new edges!");
    
    if(new_edges>0){
      vector<Edge>new_active_edges = {active_edges.back()};
      active_edges.pop_back();
      if(new_edges == 2){
        new_active_edges.push_back(active_edges.back());
        active_edges.pop_back();
      }
      vector<Edge>new_checked_edges = checked_edges;
      new_checked_edges.insert(new_checked_edges.end(), active_edges.begin(), active_edges.end());
      
      starting(my_mesh, new_active_edges, new_checked_edges, F, e_size);
    }
    // my_mesh.cout_triangles();
    if (round % 10 == 0) {
      my_mesh.cout_triangles_number();
      cout << "Number of edges in active_edges: " << active_edges.size()
           << endl;
      cout << endl;
    }
    my_mesh.obj_format();
  }

}

// a
// b
// c

Mesh BasicAlgorithm::calculate() {
  

  /*
    Point T11(numeric(0), numeric(0), numeric(0)), T12(numeric(1), numeric(0),
    numeric(0)), T13(numeric(0), numeric(-1), numeric(0)); Triangle T1(T11, T12,
    T13); Point P1(numeric(1), numeric(-1), numeric(0)); Point P2(numeric(0),
    numeric(-1), numeric(0)); Point P3(numeric(-1), numeric(-1), numeric(0));
    Point P4(numeric(-1), numeric(0), numeric(0));
    Point P5(numeric(-1), numeric(1), numeric(0));
    Point P6(numeric(0), numeric(1), numeric(0));
    Point P7(numeric(1), numeric(1), numeric(0));

    cout<< "Angle1: " << angle(Edge(T11, T12), P1, T1) <<endl;
    cout<< "Angle2: " << angle(Edge(T11, T12), P2, T1) <<endl;
    cout<< "Angle3: " << angle(Edge(T11, T12), P3, T1) <<endl;
    cout<< "Angle4: " << angle(Edge(T11, T12), P4, T1) <<endl;
    cout<< "Angle5: " << angle(Edge(T11, T12), P5, T1) <<endl;
    cout<< "Angle6: " << angle(Edge(T11, T12), P6, T1) <<endl;
    cout<< "Angle7: " << angle(Edge(T11, T12), P7, T1) <<endl;
  */
  starting(my_mesh, active_edges, checked_edges, F, e_size);

  active_edges.clear();
  active_edges = checked_edges;
  checked_edges.clear();

  cout << "Fixing holes!" << endl;
  cout << "Number of active edges: " << active_edges.size() << endl;

  ending(my_mesh, active_edges, checked_edges, F, e_size);

  my_mesh.obj_format();
  // my_mesh.output();
  return my_mesh;
}

#endif