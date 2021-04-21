#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "function.h"
#include "assertm.h"
#include "mesh.h"

class Mesh;

// N-R method for root finding, not necessarily the closest root
numeric Newton_Raphson(const realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point);

// bisection of function f over interval of 2 points
numeric Bisect(const realsymbol my_x, const ex &f, const numeric point1,
               const numeric point2, int iter);

// Bisection is called when N-R proejcts to distant point
// finds 2 points on opposite sides of surface and returns result of bisection
// on these two points
numeric Bisection(const realsymbol my_x, const ex &f, numeric starting_point,
                  numeric e_size);

// returns projected point in the direction of normal
Point project(Point point_to_project, Vector normal, const Function &F,
              const numeric e_size);

// connects two vectors of edges
vector<Edge> connect_edges(const vector<Edge> &v1, const vector<Edge> &v2);

// connects two vectors of points
vector<Point> connect_points(const vector<Point> &v1, const vector<Point> &v2);

// angle BAP in range (-Pi, Pi) with respect to neighbour triangle
numeric angle(const Edge &working_edge, const Point P, const Triangle &N) ;

// true if angle is between 0 and 9*pi/10 with respect to neighbour triangle
bool good_orientation(const Edge &working_edge, const Point P,
                      const Triangle &N);

// https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d

// returns ditance between point and line segment given by working edge
numeric line_point_dist(const Edge &working_edge, const Point P,
                        const Triangle &neighbour_triangle) ;

// true if edge is active
bool is_active(const Edge &edge, const vector<Edge> &active_edges);

// true if edge is checked
bool is_checked(const Edge &edge, const vector<Edge> &checked_edges);

// true if edge is active or checked
bool is_bounding(const Edge &edge, const BoundingBox &bounding_box);

bool is_border(const Edge &edge, const vector<Edge> &active_edges,
               const vector<Edge> &checked_edges,
               const BoundingBox &bounding_box) ;

// true if point is on border of mesh
bool is_border_point(Point P, const vector<Edge> &active_edges,
                     const vector<Edge> &checked_edges,
                     const BoundingBox &bounding_box);

// throws error if it is found more than once
void delete_from_active(const Edge &edge, vector<Edge> &active_edges) ;

// throws error if it is found more than once
void delete_from_checked(const Edge &edge, vector<Edge> &checked_edges);

// throws error if it is found more than once
void delete_from_bounding(const Edge &edge, BoundingBox &bounding_box);

// throws error if it is already there
void push_edge_to_active(const Edge &edge, vector<Edge> &active_edges);

// throws error if it is already there
void push_edge_to_checked(const Edge &edge, vector<Edge> &checked_edges);

void push_edge_to_bounding(const Edge &edge, BoundingBox &bounding_box);

// checks if edges of new triangle are active or are not im mesh
bool good_edges(const Mesh &my_mesh, const vector<Edge> &active_edges,
                const vector<Edge> &checked_edges, const Edge &working_edge,
                const Point &P, const BoundingBox &bounding_box);

// finds closest border point to edge
std::optional<Point> get_closest_point(const Mesh &my_mesh,
                                       const vector<Edge> &active_edges,
                                       const vector<Edge> &checked_edges,
                                       const Edge &working_edge,
                                       const Triangle &N, const numeric &e_size,
                                       const BoundingBox &bounding_box);

// finds closest border edge to point P
std::optional<pair<Edge, numeric>> get_closest_edge(
    const vector<Edge> &active_edges, const vector<Edge> &checked_edges,
    const BoundingBox &bounding_box, const Point &P, const Triangle &N);

// Returns unit vector in the plane of triangle T, pointing outside from T from
// the midpoint of edge e, perpendicular to e
Vector find_direction(Edge e, const Triangle &T, numeric e_size);

// finds neighbour of prev/next which has the smallest angle with the working
// edge
pair<std::optional<Point>, std::optional<Point>>
find_closest_prev_next(const Mesh &my_mesh, const Edge &working_edge,
                       const vector<Point> &prev, const vector<Point> &next);

pair<Point, Point> find_prev_next(const Mesh &my_mesh, const Edge &working_edge,
                                  const vector<Edge> &active_edges,
                                  const vector<Edge> &checked_edges,
                                  const BoundingBox &bounding_box);

bool is_vertex_good_possibility(const Point candidate, const Point prev,
                                const Point next, const Edge &working_edge,
                                const Triangle &neighbour_triangle,
                                const vector<Edge> &active_edges,
                                const vector<Edge> &checked_edges,
                                const Mesh &my_mesh, const Function &F,
                                const BoundingBox &bounding_box);

Point get_projected(const Edge &working_edge, const vector<Edge> &active_edges,
                    const vector<Edge> &checked_edges, const numeric e_size,
                    const Mesh &my_mesh, const Function F,
                    const BoundingBox &bounding_box);

#endif
