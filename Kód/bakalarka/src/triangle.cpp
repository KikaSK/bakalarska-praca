#include "triangle.h"

#include "assertm.h"
#include <iostream>

Triangle::Triangle(Point A, Point B, Point C)
    : _A(A), _B(B), _C(C), _e1(Edge(_A, _B)), _e2(Edge(_B, _C)),
      _e3(Edge(_C, _A)){
          // if (!(_e1.get_length() + _e2.get_length() - _e3.get_length() >
          // -0.001  &&
          //_e2.get_length() + _e3.get_length() - _e1.get_length() > -0.001 &&
          //_e1.get_length() + _e3.get_length() - _e2.get_length() > -0.001 ))

          // throw "Given points don't represent a Triangle";
      };
Triangle::Triangle() = default;

Point Triangle::get_gravity_center() const {
  return Point((_A.x() + _B.x() + _C.x()) / 3, (_A.y() + _B.y() + _C.y()) / 3,
               (_A.z() + _B.z() + _C.z()) / 3);
}
bool Triangle::is_triangle() const {
  // testing triangle inequality
  if (_e1.get_length() + _e2.get_length() > _e3.get_length() &&
      _e2.get_length() + _e3.get_length() > _e1.get_length() &&
      _e1.get_length() + _e3.get_length() > _e2.get_length())
    return true;
  return false;
}
Vector Triangle::get_normal() const {
  Vector AB(_A, _B);
  Vector AC(_A, _C);

  return (AB ^ AC).unit();
}
// using formula from webpage:
// https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
Point Triangle::get_circumcenter() const {
  Vector AB(_A, _B);
  Vector AC(_A, _C);

  return Point(_A, (AC.get_length_squared() * ((AB ^ AC) ^ AB) +
                    AB.get_length_squared() * ((AC ^ AB) ^ AC)) /
                       (2 * (AB ^ AC).get_length_squared()));
}