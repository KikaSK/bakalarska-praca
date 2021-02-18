#ifndef EDGE_H
#define EDGE_H

#include "point.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>

using namespace GiNaC;
using std::endl;

class Edge {
private:
  Point _A, _B;

public:
  Edge(Point A, Point B);
  Edge();
  Edge(Point A, Vector u);

  Point A() const;
  Point B() const;

  numeric get_length() const;
  Point get_midpoint() const;

  friend std::ostream &operator<<(std::ostream &os, const Edge &e) {
    os << "A: " << e._A << endl << "B: " << e._B << endl;
    return os;
  }
};

#endif
