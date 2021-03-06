#ifndef MESH_H
#define MESH_H

#include "assertm.h"
#include "edge.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
#include <optional>

using namespace GiNaC;
using std::endl;
using std::pair;
using std::vector;

// co chcem robit
// pridavat trojuholniky, hrany, body
// vyberat a pridavat do active_edges
// mat normaly trojuholnikov ku ktorym viem pristupovat
// mat hrany, vediet povedat, ze tuto hranu maju tieto 1-2 trojuholniky
// mat body, vediet povedat ktore hrany vychadzaju z bodu
// vediet vypocitat vektor kolmy na hranu v rovine susedneho trojuholnika

class Mesh {
private:
  vector<Point> _mesh_points;
  vector<Triangle> _mesh_triangles;
  vector<pair<Edge, Triangle>> _mesh_edges;
  // vector< pair<Edge, Triangle> > _mesh_active_edges;
public:
  Mesh(Triangle T);
  Mesh() = delete;

  void cout_triangles() const;
  void add_triangle(Edge e, Point P);
  Triangle find_triangle_with_edge(Edge e) const;
  bool check_Delaunay(Triangle T) const;
  vector<Point> get_breakers(Triangle T) const;
  void output() const;
  void cout_triangles_number() const;
  void obj_format() const;
  std::optional<vector<Point>> empty_surrounding(Point P, numeric e_size) const;
};

#endif
