#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "edge.h"
#include "triangle.h"
#include "mesh.h"
#include "function.h"

#include <vector>

using std::vector;

//N-R method for root finding, not necessarily the closest root
numeric Newton_Raphson(const realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point);

//bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter);

//Bisection is called when N-R proejcts to distant point
//finds 2 points on opposite sides of surface and returns result of bisection
//on these two points
numeric Bisection(const realsymbol my_x, const ex &f, numeric starting_point,
                  numeric e_size);

//returns projected point in the direction of normal
Point project(Point point_to_project, Vector normal, const Function &F,
              const numeric e_size);

//connects two vectors of edges
vector<Edge> connect_edges(const vector<Edge> &v1, const vector<Edge> &v2);

//connects two vectors of points
vector<Point> connect_points(const vector<Point> &v1, const vector<Point> &v2);

//angle BAP in range (-Pi, Pi) with respect to neighbour triangle
numeric angle(const Edge &working_edge, const Point P, const Triangle &N);

// true if angle is between 0 and 3*pi/4 with respect to neighbour triangle
bool good_orientation(const Edge &working_edge, const Point P,
                      const Triangle &N);

// https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
// returns ditance between point and line given by working edge
numeric line_point_dist(const Edge &working_edge, const Point P,
                        const Triangle &neighbour_triangle);

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(Edge e, const Triangle &T, numeric e_size);

#endif
