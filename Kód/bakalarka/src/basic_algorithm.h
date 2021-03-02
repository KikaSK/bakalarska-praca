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

Point project(Point point_to_project, Vector normal, const Function &F) {

  realsymbol my_x("my_x");

  numeric starting_point;

  // for better orientation
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

  numeric projected_point = Newton_Raphson(my_x, f, df, starting_point);
  numeric projected_x =
      ex_to<numeric>(param_x.subs(my_x == projected_point).evalf());
  numeric projected_y =
      ex_to<numeric>(param_y.subs(my_x == projected_point).evalf());
  numeric projected_z =
      ex_to<numeric>(param_z.subs(my_x == projected_point).evalf());

  assertm(F.substitute(lst{F.get_x() == projected_x, F.get_y() == projected_y,
                           F.get_z() == projected_z}) < 10e-6,
          "Projected point not on surface!");

  return Point(projected_x, projected_y, projected_z);
}

// a
// b
// c

Vector find_direction(Edge e, Triangle &T, numeric e_size) {
  Vector normal = T.get_normal();
  Vector edge_vector(e.A(), e.B());

  Vector direction = (normal ^ edge_vector).unit();

  // cout<< "Inside function, my directon: " << direction <<endl;
  numeric delta = e_size / 4;

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
    if (!T.is_in_triangle(Point(e.get_midpoint(), -delta * direction))) {
      assertm(false, "No points in triangle!");
    }
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
  return false;
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

// a
// b
// c

void delete_from_active(const Edge &edge, vector<Edge> &active_edges) {
  int counter = 0;
  int index = -1;
  for (int i = 0; i < active_edges.size(); ++i) {
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
    std::swap(active_edges[index], active_edges.back());
    active_edges.pop_back();
  }
  return;
}

// a
// b
// c

void delete_from_checked(const Edge &edge, vector<Edge> &checked_edges) {
  int counter = 0;
  int index = -1;
  for (int i = 0; i < checked_edges.size(); ++i) {
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
    std::swap(checked_edges[index], checked_edges.back());
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

/*
thought process:
we want to make a prev triangle if the oriented angle (BAP) is less than 90°
the first part ensures interval (-90°, 90°)
the second part ensures interval (0, 180)
if boths conditions are satisfied we have interval
*/

bool good_orientation(const Edge &working_edge, const Point P,
                      const Triangle &N) {

  // first part
  Triangle T = Triangle(working_edge.A(), working_edge.B(), P);
  Vector direction = (-1) * working_edge.get_length() / 2 *
                     find_direction(working_edge, T, working_edge.get_length());
  Point test_point(working_edge.get_midpoint(), direction);

  // second part

  // the third point in neighbour triangle
  std::optional<Point> NPoint = std::nullopt;

  if (T.A() != working_edge.A() && T.A() != working_edge.B())
    NPoint = T.A();

  else if (T.B() != working_edge.A() && T.B() != working_edge.B())
    NPoint = T.B();

  else if (T.C() != working_edge.A() && T.C() != working_edge.B())
    NPoint = T.C();

  assertm(NPoint.has_value(), "Not found neighbour point!");

  // if cross products points in the opposite direction all is good, else we
  // refuse this triangle
  Vector AB = Vector(working_edge.A(), working_edge.B());
  Vector AP = Vector(working_edge.A(), P);
  Vector AN = Vector(working_edge.A(), NPoint.value());

  Vector cross1 = AB ^ AP;
  Vector cross2 = AB ^ AN;

  numeric dot = cross1 * cross2;

  return (T.is_in_triangle(test_point) && dot < 0);
}

// a
// b 
// c

pair<std::optional<Point>, std::optional<Point>> find_closest_prev_next(const Mesh & my_mesh, const Edge & working_edge, const vector<Point> &prev, const vector<Point> &next)
{
  vector<std::optional<numeric> > prev_angles;
  vector<std::optional<numeric> > next_angles;

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  for (auto prev_point : prev){
    // we will use the good oreintation function to ensure the 0, 90° interval and then find smallest angle
    if(good_orientation(working_edge, prev_point, neighbour_triangle)){
      prev_angles.push_back(90);
    }
    else{
      prev_angles.push_back(std::nullopt);
    }
  }

}

pair<Point, Point> find_prev_next(const vector<Edge> &active_edges,
                                  const vector<Edge> &checked_edges,
                                  Edge edge) {

  vector<Point> prev;
  vector<Point> next;

  int counter = 0;
  
  for (auto curr_edge : active_edges) {
    if (curr_edge.A() == edge.A()) {
      prev.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == edge.A()) {
      prev.push_back(curr_edge.A());
      ++counter;
    }
    if (curr_edge.A() == edge.B()) {
      next.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == edge.B()) {
      next.push_back(curr_edge.A());
      ++counter;
    }
  }

  for (auto curr_edge : checked_edges) {
    if (curr_edge.A() == edge.A()) {
      prev.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == edge.A()) {
      prev.push_back(curr_edge.A());
      ++counter;
    }
    if (curr_edge.A() == edge.B()) {
      next.push_back(curr_edge.B());
      ++counter;
    }
    if (curr_edge.B() == edge.B()) {
      next.push_back(curr_edge.A());
      ++counter;
    }
  }
  assertm(!prev.empty() && !next.empty(),
          "Neighbour edge not found in border edges!");
  for (auto prev_point : prev){
    assertm(prev_point != edge.B(),
          "Working edge found in border edges!");
  }
  for (auto next_point : next){
    assertm(next_point != edge.A(), "Working edge found in border edges!");
  }
  
  
  if (prev.size() == 1 && next.size() == 1)
    return pair(prev[0], next[0]);
  else 
    return find_closest_prev_next(prev, next);
  
}

// a
// b
// c

bool is_vertex_good_possibility(const Point candidate, const Point prev,
                                const Point next, const Edge &working_edge,
                                const Triangle &neighbour_triangle,
                                const vector<Edge> &active_edges,
                                const vector<Edge> &checked_edges,
                                const Mesh &my_mesh, Function F) {
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

  numeric height = working_edge.get_length() * sqrt(numeric(3)) / 2;
  assertm(abs(height - e_size * sqrt(numeric(3)) / 2) < e_size / 3,
          "Weird size of height!");

  // TODO: checkovat ci je triangle iba jeden
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  assertm(neighbour_triangle.is_triangle(), "Neighbour triangle not valid!");

  Vector direction =
      height * find_direction(working_edge, neighbour_triangle, e_size);
  // TODO: kontrolovat ci je P v rovine neighbour triangle
  Point P(center, direction);
  // TODO: checknut gradient
  Vector n_A = F.get_gradient_at_point(working_edge.A()).unit();
  Vector n_B = F.get_gradient_at_point(working_edge.B()).unit();

  Vector normal = (n_A + n_B) / 2;
  Point projected = project(P, normal, F);

  return projected;
}


// a
// b
// c

// returns true if prev/next triangle is added to mesh else returns false
bool fix_prev_next(Mesh &my_mesh, vector<Edge> &active_edges,
                   vector<Edge> &checked_edges, const Edge &working_edge,
                   const bool is_prev) {

  // find previous and next points
  auto [prev, next] = find_prev_next(active_edges, checked_edges, working_edge);
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

    cout << "Found prev triangle!" << endl;

    my_mesh.add_triangle(edge, vertex);

    assertm(is_border(Edge(edge.A(), vertex), active_edges, checked_edges),
            "Neighbour edge not in border!");
    assertm(!is_border(Edge(edge.B(), vertex), active_edges, checked_edges),
            "New edges already in border!");

    cout << "OK3" << endl;
    delete_from_active(Edge(edge.A(), vertex), active_edges);
    delete_from_checked(Edge(edge.A(), vertex), checked_edges);

    push_edge_to_active(Edge(edge.B(), vertex), active_edges);
    cout << "OK4" << endl;
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
  //assertm(false, "In check overlap!");

  auto [prev, next] = find_prev_next(active_edges, checked_edges, working_edge);

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  assertm(neighbour_triangle.is_triangle(), "Neighbour triangle not valid!");

  // checks if overlap point is not neighbour or working edge point also if
  // overlap is on the border and if overlap triangle has good orientation of
  // normal
  if (is_vertex_good_possibility(overlap_point, prev, next, working_edge,
                                 neighbour_triangle, active_edges,
                                 checked_edges, my_mesh, F))

  {
    cout << "OK1eeeeeeeeeee" << endl;
    cout << "tatatatatata" << endl;
    Triangle maybe_new_T =
        Triangle(working_edge.A(), working_edge.B(), overlap_point);
    cout << "OK2eeeeeee" << endl;

    // if Delaunay constraint is satisfied add the triangle to
    // triangulation and end
    if (maybe_new_T.is_triangle() && my_mesh.check_Delaunay(maybe_new_T)) {

      //assertm(false, "Found overlap triangle, hooray!");
      cout << "Found overlap triangle" << endl;
      my_mesh.add_triangle(working_edge, overlap_point);

      Edge new_edge1 = Edge(working_edge.A(), overlap_point);
      Edge new_edge2 = Edge(working_edge.B(), overlap_point);
      if (!is_border(new_edge1, active_edges, checked_edges)) {
        push_edge_to_active(new_edge1, active_edges);
      } else {
        delete_from_active(new_edge1, active_edges);
        delete_from_checked(new_edge1, checked_edges);
      }
      if (!is_border(new_edge2, active_edges, checked_edges)) {
        push_edge_to_active(new_edge1, active_edges);
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

  // checks if there are some points very close to projected point
  if (my_mesh.empty_surrounding(projected, e_size).has_value()) {

    // points closer to projected point than 0.3*e_size sorted from closest
    vector<Point> close_points =
        my_mesh.empty_surrounding(projected, e_size).value();

    for (auto close_point : close_points) {
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), close_point);

      if (maybe_new_T.is_triangle()) {

        if (my_mesh.check_Delaunay(maybe_new_T)) {

          auto [prev, next] =
              find_prev_next(active_edges, checked_edges, working_edge);

          // if close point is prev we want to try fix prev
          if (close_point == prev) {
            if (fix_prev_next(my_mesh, active_edges, checked_edges,
                              working_edge, true))
              return true;
          }
          // if close point is next we want to try fix next
          else if (close_point == next) {
            if (fix_prev_next(my_mesh, active_edges, checked_edges,
                              working_edge, false))
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
    // if there are close points but nothing worked we want to try to construct
    // original triangle
    cout << "There are close points but nothing worked!" << endl;
  } else
    cout << "No close points!" << endl;

  if (my_mesh.check_Delaunay(maybe_new_T)) {

    cout << "Found new triangle" << endl;
    my_mesh.add_triangle(working_edge, projected);

    assertm(!is_border(Edge(working_edge.A(), projected), active_edges,
                       checked_edges) &&
                !is_border(Edge(working_edge.B(), projected), active_edges,
                           checked_edges),
            "New edges already in border!");

    cout << "OK1" << endl;
    push_edge_to_active(Edge(working_edge.A(), projected), active_edges);
    push_edge_to_active(Edge(working_edge.B(), projected), active_edges);
    cout << "OK2" << endl;
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
          const numeric e_size, const Function &F, realsymbol x, realsymbol y,
          realsymbol z) {

  my_mesh.obj_format();
  assertm(!is_border(working_edge, active_edges, checked_edges),
          "Workind edge found in border edges!");

  if (fix_proj(my_mesh, active_edges, checked_edges, working_edge, e_size, F))
    return;
  else if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                         true))
    return;
  else if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                         false))
    return;
  else {
    cout << "No triangle found" << endl;
    // if no triangle was created
    assertm(!is_border(working_edge, active_edges, checked_edges),
            "Something very wrong!");
    push_edge_to_checked(working_edge, checked_edges);
    return;
  }
}
// a
// b
// c

Mesh BasicAlgorithm::calculate() {
  while (!active_edges.empty()) {
    std::optional<Edge> working_edge = std::nullopt;

    std::random_shuffle(active_edges.begin(), active_edges.end());

    working_edge = active_edges.back();

    assertm(working_edge.has_value(), "No working edge!");
    active_edges.pop_back();

    // cout << "Current working edge: " << endl << working_edge.value() <<
    // endl;

    step(my_mesh, active_edges, checked_edges, working_edge.value(), e_size, F,
         x, y, z);
    // my_mesh.cout_triangles();
    my_mesh.cout_triangles_number();
    cout << "Number of edges in active_edges: " << active_edges.size() << endl;
    cout << endl;
  }
  my_mesh.output();
  return my_mesh;
}

#endif