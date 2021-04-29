#ifndef EDGE_H
#define EDGE_H

#include "point.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>

using namespace GiNaC;
using std::endl;

class Triangle;

class Edge {
private:
  Point _A, _B;

public:
  Edge(Point A, Point B);
  Edge(Point A, Vector u);
  Edge() = delete;

  Point A() const;
  Point B() const;

  numeric get_length() const;
  Point get_midpoint() const;
  //void set_neighbour_triangle(const Triangle &T);
  //std::pair<std::optional<Triangle>, std::optional<Triangle>> get_neighbour_triangles();

  friend bool operator==(const Edge &e1, const Edge &e2) {
    return (e1._A == e2._A && e1._B == e2._B) ||
           (e1._A == e2._B && e1._B == e2._A);
  }
  friend bool operator!=(const Edge &e1, const Edge &e2) {
    return !((e1._A == e2._A && e1._B == e2._B) ||
             (e1._A == e2._B && e1._B == e2._A));
  }

  friend std::ostream &operator<<(std::ostream &os, const Edge &e) {
    os << "A: " << e._A << endl << "B: " << e._B << endl;
    return os;
  }
  /*
  private:
    std::pair<std::optional<Triangle>, std::optional<Triangle>> neighbour_triangles;
  */
};

#endif
