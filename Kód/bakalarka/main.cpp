#include <cassert>
#include <cmath>
#include <ginac/ginac.h>
#include <iostream>
#include <ostream>
#include <queue>
#include <vector>

#include "algorithms.h"
#include "assertm.h"
#include "basic_algorithm.h"
#include "edge.h"
#include "function.h"
#include "mesh.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"
#include "bounding_box.h"

using std::cout;
using std::endl;
using std::queue;
using std::vector;
using namespace GiNaC;

// substitutes to function and returns numeric
numeric substitute(const Function F, GiNaC::ex il) {
  return ex_to<numeric>(F.get_function().subs(il).evalf());
}

// finds first edge from seed point
Edge get_seed_edge(Point seed_point, const Function &F, numeric edge_size, const BoundingBox & bounding_box) {

  Vector edge_size_tangent =
      edge_size * (F.get_tangent_at_point(seed_point).unit());

  Point point_to_project(seed_point, edge_size_tangent);

  // direction of projection
  Vector direction = F.get_gradient_at_point(point_to_project).unit();

  Point projected_point = project(point_to_project, direction, F, {edge_size});
  projected_point = bounding_box.crop_to_box(seed_point, projected_point, edge_size);

  assertm(seed_point != projected_point, "Error in get_seed_edge");

  return Edge(seed_point, projected_point);
}

// finds third point in first triangle from seed edge
Point get_seed_triangle(const Edge &e, numeric edge_size, const Function &F, const BoundingBox &bounding_box) {

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

  Vector normal = F.get_gradient_at_point(point_to_project).unit();

  Point projected = project(point_to_project, normal, F, {edge_size});
  projected = bounding_box.crop_to_box(e.get_midpoint(), projected, edge_size);
  return projected;
}

// returns first triangle
Triangle find_seed_triangle(const Function &F, Point seed, numeric e_size, BoundingBox bounding_box) {

  Vector normal = F.get_gradient_at_point(seed).unit();
  // project point on surface just to be sure it is lying on the surface with
  // enough precision
  seed = project(seed, normal, F, {e_size});
  assertm(bounding_box.is_inside(seed), "Seed point outside of bounding box!");
  // gets seed edge
  Edge seed_edge = get_seed_edge(seed, F, e_size, bounding_box);

  // gets third point in seed triangle
  Point Q = get_seed_triangle(seed_edge, e_size, F, bounding_box);

  // return seed triangle
  return Triangle(seed_edge.A(), seed_edge.B(), Q);
}

void test_find_seed_triangle();

int main() {
  Digits = 15;

  // test_find_seed_triangle();

  realsymbol x("x"), y("y"), z("z");

  // sphere
  // OK: 0.2, 0.4, 0.6
  /*
  BoundingBox my_bounding_box(numeric(-0.5), numeric(1.5), numeric(-1.5),
                                  numeric(0.9), numeric(-0.9), numeric(0.1));
  numeric e_size = 0.2;
  ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
  vector<ex> input_dF;

  input_dF.push_back(diff(input_F, x));
  input_dF.push_back(diff(input_F, y));
  input_dF.push_back(diff(input_F, z));
  // input_dF.push_back(2 * x);
  // input_dF.push_back(2 * y);
  // input_dF.push_back(2 * z);

  Function F(x, y, z, input_F, input_dF);

  Point seed(1, 0, 0);
  
*/
  
    //egg
    //OK: 0.3, 0.6, 0.8
    /*
  BoundingBox my_bounding_box(numeric(-10), numeric(2), numeric(-1.5),
                                  numeric(5), numeric(-4), numeric(3));
      numeric e_size = 0.4;
      ex input_F =
            pow((x - 1) / 2, 2) + pow((y - 1) / 3, 2) + pow((z - 1), 2) - 1;
        vector<ex> input_dF;
        input_dF.push_back((x - 1) / 2);
        input_dF.push_back(2 * (y - 1) / 3);
        input_dF.push_back(2 * (z - 1));

        Function F(x, y, z, input_F, input_dF);
        Point seed(1, 1, 2);
*/
  
    // torus
    // OK: 5 10 15
    //max e_size = 20
    /*
    BoundingBox my_bounding_box(numeric(-60), numeric(60), numeric(-60),
                              numeric(60), numeric(-60), numeric(5));

    numeric e_size = 10;
    ex input_F = pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + 40 * 40 - 15 * 15, 2) -
                 4 * 40 * 40 * (pow(x, 2) + pow(y, 2));
    vector<ex> input_dF;
    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));

    Function F(x, y, z, input_F, input_dF);
    Point seed(55, 0, 0);
    */

   /*
    BoundingBox my_bounding_box(numeric(-10), numeric(10), numeric(-10),
                              numeric(10), numeric(-5), numeric(5));
      // plane
      numeric e_size = 2.6;

      ex input_F = z - 1;
      vector<ex> input_dF;

      input_dF.push_back(diff(input_F, x));
      input_dF.push_back(diff(input_F, y));
      input_dF.push_back(diff(input_F, z));

      Function F(x, y, z, input_F, input_dF);

      Point seed(0, 0, 1);
    */

  // genus
  /*
      BoundingBox my_bounding_box(numeric(-0.5), numeric(3), numeric(-0.6),
                              numeric(2), numeric(-10), numeric(50));
      //max size not falling: 0.1 takes very long
      numeric e_size =0.2;
      ex input_F =
      2*y*(y*y - 3*x*x)*(1-z*z)+pow((x*x+y*y), 2)-(9*z*z-1)*(1-z*z);

        vector<ex> input_dF;
        input_dF.push_back(diff(input_F, x));
        input_dF.push_back(diff(input_F, y));
        input_dF.push_back(diff(input_F, z));

        Function F(x, y, z, input_F, input_dF);
        Point seed(0, 0, 1);

  */

    // blobby

      //OK: 0.08 0.1 0.12 0.15 0.17 0.2 0.21 0.22 0.23 0.235
      //max size: 0.25
      /*
      BoundingBox my_bounding_box(numeric(-2), numeric(2), numeric(-2),
                              numeric(2), numeric(-0.09), numeric(50));

      numeric e_size = 0.14;
      ex input_F =
      sqrt((x-1)*(x-1)+y*y+z*z)*sqrt((x+1)*(x+1)+y*y+z*z)*sqrt(x*x+(y-1)*(y-1)+z*z)*sqrt(x*x+(y+1)*(y+1)+z*z)-1.1;

        vector<ex> input_dF;
        input_dF.push_back(diff(input_F, x));
        input_dF.push_back(diff(input_F, y));
        input_dF.push_back(diff(input_F, z));

        Function F(x, y, z, input_F, input_dF);
        Point seed(1.2038, 0, 0);
        seed = project(seed, F.get_gradient_at_point(seed), F, e_size);
      */

      // cyllinder
      
      BoundingBox my_bounding_box(numeric(-2), numeric(2), numeric(-2),
                              numeric(2), numeric(-5), numeric(5));

      numeric e_size = 0.5;
      ex input_F = x*x + y*y - 1;
        vector<ex> input_dF;
        input_dF.push_back(diff(input_F, x));
        input_dF.push_back(diff(input_F, y));
        input_dF.push_back(diff(input_F, z));

        Function F(x, y, z, input_F, input_dF);
        Point seed(1, 0, 0);
        seed = project(seed, F.get_gradient_at_point(seed), F, e_size);
      
  /*
  // cubedsphere
    //OK: 0.1 0.2 0.3 0.4 0.5
    //max: 1.5

      numeric e_size = 0.3;
      ex input_F = x*x*x*x + y*y*y*y + z*z*z*z - 1;

        vector<ex> input_dF;
        input_dF.push_back(diff(input_F, x));
        input_dF.push_back(diff(input_F, y));
        input_dF.push_back(diff(input_F, z));

        Function F(x, y, z, input_F, input_dF);
        Point seed(1, 0, 0);
  */

  // tetrahedron by ajko
  /*
      numeric e_size = 0.5;
      ex input_F = x*x*x*x + 2*x*x*y*y + 2*x*x*z*z + y*y*y*y + 2*y*y*z*z+z*z*z*z
     + 8*x*y*z - 10*x*x - 10*y*y - 10*z*z + 20; vector<ex> input_dF;
      input_dF.push_back(diff(input_F, x));
      input_dF.push_back(diff(input_F, y));
      input_dF.push_back(diff(input_F, z));

      Function F(x, y, z, input_F, input_dF);
      Point seed(1.66250775, 0, 0);

  */

  // akokeby 4 preepojene gule
  /*
   numeric e_size = 0.3;
      ex input_F = x*x*x*x-5*x*x + y*y*y*y-5*y*y+z*z*z*z-5*z*z + 11.8;

        vector<ex> input_dF;
        input_dF.push_back(diff(input_F, x));
        input_dF.push_back(diff(input_F, y));
        input_dF.push_back(diff(input_F, z));

        Function F(x, y, z, input_F, input_dF);
        Point seed(numeric(-2.26634),-numeric(1.58114), -numeric(1.58114));
  */
  /*
  cout << "Side lenghts of seed triangle: " << endl
       << seed_triangle.AB().get_length() << " "
       << seed_triangle.BC().get_length() << " "
       << seed_triangle.CA().get_length() << endl;
  */

  // making seed triangle from seed point lying on the surface
  Triangle seed_triangle = find_seed_triangle(F, seed, e_size, my_bounding_box);

  assertm(seed_triangle.AB() != seed_triangle.BC() &&
              seed_triangle.AB() != seed_triangle.CA() &&
              seed_triangle.BC() != seed_triangle.CA(),
          "Seed triangle contains duplicit edges!");

  BasicAlgorithm alg(F, seed_triangle, e_size, x, y, z, my_bounding_box);

  alg.calculate();
}

void test_find_seed_triangle() {
  realsymbol x("x"), y("y"), z("z");
  numeric e_size = 0.01;

  {
    BoundingBox test_bounding_box = BoundingBox(-1.5, 1.5, -1.5, 1.5, -1.5, 1.5);
    ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
    vector<ex> input_dF;
    input_dF.push_back(2 * x);
    input_dF.push_back(2 * y);
    input_dF.push_back(2 * z);

    Function F(x, y, z, input_F, input_dF);
    Point seed(1, 0, 0);

    Triangle t = find_seed_triangle(F, seed, e_size, test_bounding_box);

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
    BoundingBox test_bounding_box = BoundingBox(-1.5, 1.5, -1.5, 1.5, -1.5, 1.5);
    ex input_G =
        pow((x - 1) / 2, 2) + pow((y - 1) / 3, 2) + pow((z - 1), 2) - 1;
    vector<ex> input_dG;
    input_dG.push_back((x - 1) / 2);
    input_dG.push_back(2 * (y - 1) / 3);
    input_dG.push_back(2 * (z - 1));

    Function G(x, y, z, input_G, input_dG);
    Point seed(1, 1, 2);

    Triangle t = find_seed_triangle(G, seed, e_size, test_bounding_box);
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
