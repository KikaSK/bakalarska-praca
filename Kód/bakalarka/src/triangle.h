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
  Triangle() = delete;

  Point A() const;
  Point B() const;
  Point C() const;

  Edge AB() const;
  Edge BC() const;
  Edge CA() const;

  Point get_gravity_center() const;
  bool is_triangle() const;
  Point get_circumcenter() const;
  Vector get_normal() const;
  bool is_in_triangle(Point P) const;

  friend bool operator==(const Triangle &a, const Triangle &b) {
    return (a._A == b._A && a._B == b._B && a._C == b._C);
  }

  friend std::ostream &operator<<(std::ostream &os, const Triangle &T) {
    os << "A: " << T._A << endl
       << "B: " << T._B << endl
       << "C: " << T._C << endl;
    return os;
  }
};

#endif