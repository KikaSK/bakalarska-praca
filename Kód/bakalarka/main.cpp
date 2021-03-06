#include <cassert>
#include <cmath>
#include <ginac/ginac.h>
#include <iostream>
#include <ostream>
#include <queue>
#include <vector>

#include "assertm.h"
#include "basic_algorithm.h"
#include "edge.h"
#include "function.h"
#include "mesh.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"

using std::cout;
using std::endl;
using std::queue;
using std::vector;
using namespace GiNaC;

numeric substitute(const Function F, GiNaC::ex il) {
  return ex_to<numeric>(F.get_function().subs(il).evalf());
}

Edge get_seed_edge(Point seed_point, const Function &F, numeric edge_size) {

  // point to project
  Vector edge_size_tangent =
      edge_size * (F.get_tangent_at_point(seed_point).unit());

  Point point_to_project(seed_point, edge_size_tangent);

  // direction of projection
  Vector normal = F.get_gradient_at_point(seed_point).unit();

  Point projected_point = project(point_to_project, normal, F);

  assertm(seed_point != projected_point, "Error in get_seed_edge");

  return Edge(seed_point, projected_point);
}

Point get_seed_triangle(const Edge &e, numeric edge_size, const Function &F) {

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

Triangle find_seed_triangle(const Function &F, Point seed, numeric e_size) {
  Edge seed_edge = get_seed_edge(seed, F, e_size);
  Point projected = seed_edge.B();
  ex my_func = F.get_function();
  Point Q = get_seed_triangle(seed_edge, e_size, F);

  return Triangle(seed_edge.A(), seed_edge.B(), Q);
}

void test_find_seed_triangle();

int main() {
  Digits = 15;

  // test_find_seed_triangle();

  realsymbol x("x"), y("y"), z("z");


  numeric e_size = 0.3;
  
// sphere
  ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
  vector<ex> input_dF;

    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));
  //input_dF.push_back(2 * x);
  //input_dF.push_back(2 * y);
  //input_dF.push_back(2 * z);

  Function F(x, y, z, input_F, input_dF);

  Point seed(1, 0, 0);

  
/*
  
//egg

  numeric e_size = 0.35;
  ex input_F =
        pow((x - 1) / 2, 2) + pow((y - 1) / 3, 2) + pow((z - 1), 2) - 1;
    vector<ex> input_dF;
    input_dF.push_back((x - 1) / 2);
    input_dF.push_back(2 * (y - 1) / 3);
    input_dF.push_back(2 * (z - 1));

    Function F(x, y, z, input_F, input_dF);
    Point seed(1, 1, 2);
  

*/
/*

// torus

  numeric e_size = 10;
  ex input_F =
        pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + 40*40 - 15*15, 2) - 4*40*40*(pow(x, 2) + pow(y, 2));
    vector<ex> input_dF;
    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));

    Function F(x, y, z, input_F, input_dF);
    Point seed(55, 0, 0);
  

*/

/*
// genus

  numeric e_size =0.1;
  ex input_F =
  2*y*(y*y - 3*x*x)*(1-z*z)+pow((x*x+y*y), 2)-(9*z*z-1)*(1-z*z);
       
    vector<ex> input_dF;
    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));

    Function F(x, y, z, input_F, input_dF);
    Point seed(0, 0, 1);

*/

  Triangle seed_triangle = find_seed_triangle(F, seed, e_size);

  cout << "Side lenghts of seed triangle: " << endl
       << seed_triangle.AB().get_length() << " "
       << seed_triangle.BC().get_length() << " "
       << seed_triangle.CA().get_length() << endl;

  assertm(seed_triangle.AB() != seed_triangle.BC() &&
              seed_triangle.AB() != seed_triangle.CA() &&
              seed_triangle.BC() != seed_triangle.CA(),
          "Seed triangle contains duplicit edges!");

  BasicAlgorithm alg(F, seed_triangle, e_size, x, y, z);

  alg.calculate();
}

void test_find_seed_triangle() {
  realsymbol x("x"), y("y"), z("z");
  numeric e_size = 0.01;

  {
    ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
    vector<ex> input_dF;
    input_dF.push_back(2 * x);
    input_dF.push_back(2 * y);
    input_dF.push_back(2 * z);

    Function F(x, y, z, input_F, input_dF);
    Point seed(1, 0, 0);

    Triangle t = find_seed_triangle(F, seed, e_size);

    auto total_length =
        t.AB().get_length() + t.BC().get_length() + t.CA().get_length();
    assertm(abs(total_length - e_size * 3) < 10e-2,
            "Triangle not of specified length.");

    numeric value1 =
        F.substitute(lst{x == t.A().x(), y == t.A().y(), z == t.A().z()});
    numeric value2 =
        F.substitute(lst{x == t.B().x(), y == t.B().y(), z == t.B().z()});
    numeric value3 =
        F.substitute(lst{x == t.C().x(), y == t.C().y(), z == t.C().z()});

    assertm(value1 < 10e-6 && value2 < 10e-6 && value3 < 10e-6, "Bad test!");
  }
  {
    ex input_G =
        pow((x - 1) / 2, 2) + pow((y - 1) / 3, 2) + pow((z - 1), 2) - 1;
    vector<ex> input_dG;
    input_dG.push_back((x - 1) / 2);
    input_dG.push_back(2 * (y - 1) / 3);
    input_dG.push_back(2 * (z - 1));

    Function G(x, y, z, input_G, input_dG);
    Point seed(1, 1, 2);

    Triangle t = find_seed_triangle(G, seed, e_size);
    auto total_length =
        t.AB().get_length() + t.BC().get_length() + t.CA().get_length();
    assertm(abs(total_length - e_size * 3) < 10e-2,
            "Triangle not of specified length.");

    numeric value1 =
        G.substitute(lst{x == t.A().x(), y == t.A().y(), z == t.A().z()});
    numeric value2 =
        G.substitute(lst{x == t.B().x(), y == t.B().y(), z == t.B().z()});
    numeric value3 =
        G.substitute(lst{x == t.C().x(), y == t.C().y(), z == t.C().z()});

    assertm(value1 < 10e-6 && value2 < 10e-6 && value3 < 10e-6, "Bad test!");
  }
}
