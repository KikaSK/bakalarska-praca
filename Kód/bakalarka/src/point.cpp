#include "point.h"

Point::Point(numeric x, numeric y, numeric z) : _x(x), _y(y), _z(z){};
Point::Point() = default;
Point::Point(const Point &v) : _x(v.x()), _y(v.y()), _z(v.z()){};
Point::Point(Point A, Vector u)
    : _x(A.x() + u.x()), _y(A.y() + u.y()), _z(A.z() + u.z()){};

numeric Point::x() const { return _x; }
numeric Point::y() const { return _y; }
numeric Point::z() const { return _z; }
