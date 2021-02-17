#include "edge.h"

numeric Edge::get_length() {
  numeric length =
      sqrt(pow(_A.x() - _B.x(), numeric(2)) + pow(_A.y() - _B.y(), numeric(2)) +
           pow(_A.z() - _B.z(), numeric(2)));
  return length;
}
Edge::Edge(Point A, Point B) : _A(A), _B(B) {
  if (_A.x() == _B.x() && _A.y() == _B.y() && _A.z() == _B.z())
    throw("Trying to make edge from the same points!");
}

Edge::Edge() = default;

Edge::Edge(Point A, Vector u)
    : _A(A), _B(Point(A.x() + u.x(), A.y() + u.y(), A.z() + u.z())) {
  if (u.x() == 0 && u.y() == 0 && u.z() == 0)
    throw("Trying to make edge from zero vector!");
}

Point Edge::A() const { return _A; }
Point Edge::B() const { return _B; }

Point Edge::get_midpoint() {
  Point mid((_A.x() + _B.x()) / 2, (_A.y() + _B.y()) / 2,
            (_A.z() + _B.z()) / 2);
  return mid;
}