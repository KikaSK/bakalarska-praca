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
  Point();
  Point(const Point &v);
  Point(Point A, Vector u);

  numeric x() const;
  numeric y() const;
  numeric z() const;

  friend std::ostream &operator<<(std::ostream &os, const Point &a) {
    os << '[' << a._x << ',' << a._y << ',' << a._z << ']';
    return os;
  }
};

#endif