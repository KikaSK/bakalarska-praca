#ifndef FUNCTION_H
#define FUNCTION_H

#include "assertm.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>

using namespace GiNaC;
using std::endl;
using std::vector;

class Function {
private:
  const ex _x;
  const ex _y;
  const ex _z;
  const ex _F;
  const vector<ex> _dF;

public:
  Function(const ex &x, const ex &y, const ex &z, ex F, vector<ex> dF);
  Function() = default;

  ex get_function() const;
  vector<ex> get_gradient() const;
  ex get_x() const;
  ex get_y() const;
  ex get_z() const;

  ex grad_x() const;
  ex grad_y() const;
  ex grad_z() const;

  Vector get_gradient_at_point(Point P) const;
  Vector get_tangent_at_point(Point P) const;
  numeric eval_at_point(Point P) const;
  bool is_inside(Point P) const;
  bool is_on(Point P) const;
  bool is_outside(Point P) const;
  Vector outside_normal(Triangle T) const;
  numeric substitute(GiNaC::ex il) const;
};

#endif