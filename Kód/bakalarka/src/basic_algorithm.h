#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H

#include "function.h"
#include "mesh.h"
#include "triangle.h"

#include <vector>

class BasicAlgorithm {
public:
  BasicAlgorithm(Function f, Triangle seed_triangle, numeric e_size, realsymbol x, realsymbol y, realsymbol z)
      : F(f), active_edges(), my_mesh(seed_triangle), e_size(e_size), x(x), y(y), z(z) {
    active_edges.push_back(pair(seed_triangle.AB(), false));
    active_edges.push_back(pair(seed_triangle.BC(), false));
    active_edges.push_back(pair(seed_triangle.CA(), false));
  }

  Mesh calculate();

private:
  Function F;
  vector<pair<Edge, bool>> active_edges;
  Mesh my_mesh;
  numeric e_size = 0.1;
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
  numeric prejcted_x =
      ex_to<numeric>(param_x.subs(my_x == projected_point).evalf());
  numeric prejcted_y =
      ex_to<numeric>(param_y.subs(my_x == projected_point).evalf());
  numeric prejcted_z =
      ex_to<numeric>(param_z.subs(my_x == projected_point).evalf());

  return Point(prejcted_x, prejcted_y, prejcted_z);
}

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

void push_edge_to_active(const pair<Edge, bool> &edge,
                         vector<pair<Edge, bool>> &active_edges) {

  for (auto my_edge : active_edges) {
    assertm(my_edge.first != edge.first, "Edge already in active_edges!");
  }
  active_edges.push_back(edge);
  return;
}

void delete_from_active(const Edge &edge,
                        vector<pair<Edge, bool>> &active_edges) {
  int counter = 0;
  int index = -1;
  for (auto i = 0; i < active_edges.size(); ++i) {
    if (active_edges[i].first == edge) {
      index = i;
      counter++;
    }
  }
  assertm(counter == 1 || counter == 0,
          "More than one edge found while deleting!");
  if (counter == 0)
    return;
  else {
    swap(active_edges[index], active_edges.back());
    active_edges.pop_back();
  }
  return;
}

pair<Point, Point> find_prev_next(const vector<pair<Edge, bool>> &active_edges,
                                  Edge edge) {

  numeric zero = numeric(0);
  Point prev(zero, zero, zero), next(zero, zero, zero);

  for (auto curr_edge : active_edges) {
    if (curr_edge.first.A() == edge.A()) {
      prev = curr_edge.first.B();
    }
    if (curr_edge.first.B() == edge.A()) {
      prev = curr_edge.first.A();
    }
    if (curr_edge.first.A() == edge.B()) {
      next = curr_edge.first.B();
    }
    if (curr_edge.first.B() == edge.B()) {
      next = curr_edge.first.A();
    }
  }

  assertm(prev != Point(zero, zero, zero) && next != Point(zero, zero, zero),
          "Neighbour edge not found in active_edges!");
  assertm(prev != edge.B() && next != edge.A(),
          "Working edge found in active_edges!");

  return pair(prev, next);
}

bool is_vertex_good_possibility(const Point candidate, const Point prev,
                                const Point next,
                                const Edge & working_edge,
                                const Triangle &neighbour_triangle,
                                const vector<pair<Edge, bool>> &active_edges,
                                const Mesh &my_mesh, Function F) {
  if (candidate == prev || candidate == next || candidate == working_edge.A() || candidate == working_edge.B())
    return false;

  for (auto edge_bool : active_edges) {
    Edge edge = edge_bool.first;

    // vertex found in active_edges
    if (edge.A() == candidate || edge.B() == candidate) {
      Triangle overlap_triangle = my_mesh.find_triangle_with_edge(edge);
      Vector overlap_normal = F.outside_normal(overlap_triangle);
      Vector my_normal = F.outside_normal(neighbour_triangle);
      if (overlap_normal * my_normal > 0) {
        return true;
      }
    }
  }
  return false;
}

// one step of the algorithm
void step(Mesh &my_mesh, vector<pair<Edge, bool>> &active_edges,
          const Edge &working_edge, const numeric e_size, const Function &F, realsymbol x, realsymbol y, realsymbol z) {

  Point center = working_edge.get_midpoint();
  numeric height = working_edge.get_length() * sqrt(numeric(3)) / 2;

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  Vector direction =
      height * find_direction(working_edge, neighbour_triangle, e_size);
  Point P(center, direction);
  Vector n_A = F.get_gradient_at_point(working_edge.A()).unit();
  Vector n_B = F.get_gradient_at_point(working_edge.B()).unit();

  Vector normal = (n_A + n_B) / 2;
  Point projected = project(P, normal, F);

  assertm(F.substitute(lst{x == projected.x(), y==projected.y(), z==projected.z()})<10e-6, "Projected point not on surface!");
  
  Triangle maybe_new_T(working_edge.A(), working_edge.B(), projected);

  // not satisfied Delaunay for new triangle
  if (!my_mesh.check_Delaunay(maybe_new_T)) {

    // finding prev and next edges
    pair<Point, Point> prev_next = find_prev_next(active_edges, working_edge);
    Point prev = prev_next.first;
    Point next = prev_next.second;

    maybe_new_T = Triangle(working_edge.A(), working_edge.B(), prev);

    // not satisfied Delaunay for prev triangle
    if (!my_mesh.check_Delaunay(maybe_new_T)) {
      maybe_new_T = Triangle(working_edge.A(), working_edge.B(), next);

      // not satisfied Delaunay for next triangle
      if (!my_mesh.check_Delaunay(maybe_new_T)) {

        // find points that break Delaunay
        vector<Point> breakers = my_mesh.get_breakers(maybe_new_T);
        assertm(!breakers.empty(),
                "No vertices that break Delaunay constraint!");
        // through the vertices that break Delaunay constraint
        for (auto vertex : breakers) {
          
          // try each "good" vertex
          if (is_vertex_good_possibility(vertex, prev, next, working_edge, neighbour_triangle,
                                         active_edges, my_mesh, F))
            
            {cout<< "OK1eeeeeeeeeee"<<endl;
            cout<<"tatatatatata"<<endl;
            maybe_new_T = Triangle(working_edge.A(), working_edge.B(), vertex);
            cout<< "OK2eeeeeee"<<endl;}

          // if Delaunay constraint is satisfied add the triangle to
          // triangulation and end
          if (my_mesh.check_Delaunay(maybe_new_T)) {

            cout<< "Found overlap triangle" <<endl;
            my_mesh.add_triangle(working_edge, vertex);
            push_edge_to_active(pair(Edge(working_edge.A(), vertex), false),
                                active_edges);
            push_edge_to_active(pair(Edge(working_edge.B(), vertex), false),
                                active_edges);
            return;
          }
        }


        cout<< "No triangle found" <<endl;
        // if no triangle was created
        push_edge_to_active(pair(working_edge, true), active_edges);
        return;
      }

      // the triangle next satisfies Delaunay constraint
      else {

        cout<< "Found next triangle" <<endl;
        my_mesh.add_triangle(working_edge, next);
        delete_from_active(Edge(working_edge.B(), next), active_edges);
        if (prev != next)
          push_edge_to_active(pair(Edge(working_edge.A(), next), false),
                              active_edges);
        else {
          assertm(false, "You should not get here!");
          // delete_from_active(Edge(working_edge.A(), next), active_edges);
        }
        return;
      }
    }
    // the triangle prev satisfies Delaunay constraint
    else {

     cout<< "Found prev triangle" <<endl;
      my_mesh.add_triangle(working_edge, prev);
      delete_from_active(Edge(working_edge.A(), prev), active_edges);
      if (prev != next)
        push_edge_to_active(pair(Edge(working_edge.B(), prev), false),
                            active_edges);
      else
        delete_from_active(Edge(working_edge.B(), prev), active_edges);
      return;
    }
  }
  // the found triangle satisfies Delaunay constraint
  else {
    cout<< "Found new triangle" <<endl;
    my_mesh.add_triangle(working_edge, projected);
    push_edge_to_active(pair(Edge(working_edge.A(), projected), false),
                        active_edges);
    push_edge_to_active(pair(Edge(working_edge.B(), projected), false),
                        active_edges);
    return;
  }

  push_edge_to_active(pair(working_edge, true), active_edges);
  return;
}

Mesh BasicAlgorithm::calculate() {
  while (!active_edges.empty()) {
    bool stopped = false;
    std::optional<Edge> working_edge = std::nullopt;

    if (!active_edges.back().second) {
      working_edge = active_edges.back().first;
    } else {
      for (int i = active_edges.size() - 1; i >= -1; --i) {
        if (i == -1) {
          assertm(false, "I shoud stop here!");
          active_edges.clear();
          stopped = true;
          break;
        }
        if (!active_edges[i].second) {
          std::swap(active_edges[i], active_edges.back());
          working_edge = {active_edges.back().first};
          break;
        }
      }
    }

    if (stopped)
      break;
    
    assertm(working_edge.has_value(), "No working edge!");

    assertm(!active_edges.back().second, "Something is odd!");

    assertm( working_edge == active_edges.back().first, "Wrong swapping!");

    //cout << "Current working edge: " << endl << working_edge.value() << endl;
    active_edges.pop_back();
    for (int i = 0; i < active_edges.size(); ++i) {
      if (active_edges[i].first == working_edge.value()) {
        assertm(false, "Still duplicate edges!");
      }
    }

    step(my_mesh, active_edges, working_edge.value(), e_size, F, x, y, z);
    // my_mesh.cout_triangles();
    my_mesh.cout_triangles_number();
    cout<< "Number of edges in active_edges: " <<active_edges.size() <<endl;
    if(active_edges.size() == 50)
    {
        my_mesh.output();
        return my_mesh;
    }
    cout << endl;
  }
  my_mesh.output();
  return my_mesh;
}

#endif