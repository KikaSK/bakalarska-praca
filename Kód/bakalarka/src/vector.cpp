#include "vector.h"
#include "assertm.h"

numeric is_zero_precision = 10e-8;

Vector::Vector(numeric x, numeric y, numeric z) : _x(x), _y(y), _z(z){};
Vector::Vector(const Vector &v) : _x(v.x()), _y(v.y()), _z(v.z()){};
Vector::Vector(Point A, Point B)
    : _x(B.x() - A.x()), _y(B.y() - A.y()), _z(B.z() - A.z()){};

numeric Vector::x() const { return _x; }
numeric Vector::y() const { return _y; }
numeric Vector::z() const { return _z; }

numeric Vector::get_length() const { return sqrt(get_length_squared()); }
numeric Vector::get_length_squared() const {
  return pow(_x, numeric(2)) + pow(_y, numeric(2)) + pow(_z, numeric(2));
}
Vector Vector::unit() const {
  assertm(get_length() > 10e-6, "Uniting zero vector");
  return Vector(*this / get_length());
}
Vector Vector::vector_inverse() const { return Vector(-_x, -_y, -_z); }
bool Vector::is_zero() const { return (get_length() < is_zero_precision); }
Vector Vector::get_any_perpendicular() const {
  assertm(!is_zero(), "Trying to find prependicular to zero vector!");
  if (_x == 0 && _y == 0)
    return Vector(_z, 0, -_x);
  return Vector(_y, -_x, 0);
};
