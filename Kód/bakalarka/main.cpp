#include <cassert>
#include <cmath>
#include <ginac/ginac.h>
#include <iostream>
#include <ostream>
#include <queue>
#include <vector>

#include "assertm.h"
#include "edge.h"
#include "function.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"
#include "mesh.h"

using std::cout;
using std::endl;
using std::queue;
using std::vector;
using namespace GiNaC;

numeric Newton_Raphson(realsymbol my_x, const ex & f, const ex & df, numeric starting_point) {

  numeric precision = 1e-6;

  numeric iter = starting_point;
  int iterations = 0;
  while (abs(f.subs(my_x == iter).evalf()) > precision && iterations < 1'000) {
    iterations++;
    iter -= ex_to<numeric>(f.subs(my_x == iter).evalf() /
                           df.subs(my_x == iter).evalf());
  }
  cout << "Number of iterations: " << iterations << endl;

  return iter;
}

Point project(Point point_to_project, Vector normal, const Function & F) {

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

Edge get_seed_edge(Point seed_point, const Function & F, numeric edge_size) {

  // point to project
  Vector edge_size_tangent =
      edge_size * (F.get_tangent_at_point(seed_point).unit());

  Point point_to_project(seed_point, edge_size_tangent);

  // direction of projection
  Vector normal = F.get_gradient_at_point(seed_point).unit();

  Point projected_point = project(point_to_project, normal, F);

  return Edge(seed_point, projected_point);
}

Point get_seed_triangle(const Edge &e, numeric edge_size, const Function & F) {

  Point center = e.get_midpoint();

  // normals at endpoints of the edge
  Vector n_A = F.get_gradient_at_point(e.A()).unit();
  Vector n_B = F.get_gradient_at_point(e.B()).unit();

  // average of endpoints normals
  Vector center_normal((n_A + n_B) / 2);

  Vector edge_vector(e.A(), e.B());
  Vector center_tangent = center_normal ^ edge_vector;

  assertm(abs(center_normal * center_tangent) < 1e-6, "Not perpendicular!");

  // height of equilateral triangle with side edge_size
  numeric height = sqrt(numeric(3)) / 2 * edge_size;

  Point point_to_project(center, height * center_tangent.unit());

  Point projected = project(point_to_project, center_normal, F);

  return projected;
}

Vector find_direction(Edge e, Triangle & T){
  Vector normal = T.get_normal();
  Vector edge_vector(e.A(), e.B());

  Vector direction = (normal^edge_vector).unit();

  //naprogramovat ci je bod v trojuholniku
  
}

void step(Mesh & my_mesh, const Edge & e, const Function &F) { 
  
  Point center = e.get_midpoint();

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(e);


  return; 
}

int main() {
  realsymbol x("x"), y("y"), z("z");
  Digits = 15;

  ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
  vector<ex> input_dF;
  input_dF.push_back(2 * x);
  input_dF.push_back(2 * y);
  input_dF.push_back(2 * z);

  queue<Edge> active_edges;

  Function F(x, y, z, input_F, input_dF);

  numeric e_size = 0.05;
  Point seed(1, 0, 0);
  // Point seed(sqrt(numeric(119) / 144), numeric(1) / 4, numeric(1) / 3);

  // cout<< "Seed point: " << seed <<endl;

  // cout<< "Gradient at seed point: " << F.get_gradient_at_point(seed)<<endl;
  // cout<< "Is seed outside: " << F.is_outside(seed) << endl;
  // cout<< F.get_tangent_at_point(seed)<<endl;

  Edge seed_edge = get_seed_edge(seed, F, e_size);
  Point projected = seed_edge.B();

  // cout<< "Projected point: " << projected << endl;

  ex my_func = F.get_function();

  // functional value of projected
  // cout<< "Value at projected point: " << ex_to<numeric>(my_func.subs(lst{x ==
  // projected.x(), y == projected.y(), z == projected.z()}).evalf()) <<endl;
  // cout<<endl;

  Vector dist(seed, projected);

  // cout<< "Distance from seed: " << dist.get_length();

  cout << endl;
  active_edges.push(seed_edge);

  Point Q = get_seed_triangle(seed_edge, e_size, F);

  cout << "New point: " << Q << endl;

  Vector dist1(seed, Q), dist2(projected, Q);

  cout << "Distances from original points: " << dist1.get_length() << " and "
       << dist2.get_length() << endl;
  cout << "Value at new point: "
       << ex_to<numeric>(
              my_func.subs(lst{x == Q.x(), y == Q.y(), z == Q.z()}).evalf())
       << endl;

  Edge e1(seed, Q), e2(projected, Q);
  active_edges.push(e1);
  active_edges.push(e2);

  Mesh my_mesh(Triangle(seed, Q, projected));

  while(!active_edges.empty())
  {
    Edge working_edge = active_edges.front();
    active_edges.pop();
    step(my_mesh, working_edge, F);
  }
}
