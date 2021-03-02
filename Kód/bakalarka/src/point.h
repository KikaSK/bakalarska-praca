#ifndef POINT_H
#define POINT_H

#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>

using namespace GiNaC;
using std::endl;

class Vector;

class Point {
private:
  numeric _x, _y, _z;

public:
  Point(numeric x, numeric y, numeric z);
  Point(const Point &v);
  Point(Point A, Vector u);
  Point() = delete;

  numeric x() const;
  numeric y() const;
  numeric z() const;

  friend std::ostream &operator<<(std::ostream &os, const Point &a) {
    os << '[' << a._x << ',' << a._y << ',' << a._z << ']';
    return os;
  }

  friend bool operator==(const Point &A, const Point &B) {
    auto diff = abs(A._x - B._x) + abs(A._y - B._y) + abs(A._z - B._z);
    return diff < 10e-8;
  }
  friend bool operator!=(const Point &A, const Point &B) {
    auto diff = abs(A._x - B._x) + abs(A._y - B._y) + abs(A._z - B._z);
    return diff >= 10e-8;
  }
};

#endif