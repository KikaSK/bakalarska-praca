#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "edge.h"
#include "point.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>

using namespace GiNaC;
using std::endl;

class Triangle {
private:
  Point _A, _B, _C;
  Edge _e1, _e2, _e3;

public:
  Triangle(Point A, Point B, Point C);
  Triangle();

  Point get_gravity_center() const;
  bool is_triangle() const;
  Point get_circumcenter() const;
  Vector get_normal() const;
  friend std::ostream &operator<<(std::ostream &os, const Triangle &T) {
    os << "A: " << T._A << endl
       << "B: " << T._B << endl
       << "C: " << T._C << endl;
    return os;
  }
};

#endif