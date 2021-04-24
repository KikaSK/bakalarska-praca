#include "mesh.h"

#include <fstream>

Mesh::Mesh(Triangle T) {
  _mesh_triangles.push_back(T);
  _mesh_points.push_back(T.A());
  _mesh_points.push_back(T.B());
  _mesh_points.push_back(T.C());
  _mesh_edges.push_back(pair(Edge(T.A(), T.B()), T));
  _mesh_edges.push_back(pair(Edge(T.B(), T.C()), T));
  _mesh_edges.push_back(pair(Edge(T.C(), T.A()), T));
}

//prints triangles to console
void Mesh::cout_triangles() const {
  for (size_t i = 0; i < _mesh_triangles.size() - 1; ++i) {
    std::cout << "Triangle " << i << ":" << endl;
    std::cout << _mesh_triangles[i] << endl;

    std::cout << "Edge sizes: " << _mesh_triangles[i].AB().get_length() << " "
              << _mesh_triangles[i].BC().get_length() << " "
              << _mesh_triangles[i].CA().get_length() << endl;
  }
  return;
}

// adds triangle to mesh
void Mesh::add_triangle(Edge e, Point P) {
  Triangle new_triangle(e.A(), e.B(), P);

  assertm(new_triangle.is_triangle(), "Non valid triangle adding to mesh!");

  _mesh_triangles.push_back(new_triangle);
  _mesh_points.push_back(P);
  _mesh_edges.push_back(pair(Edge(P, e.A()), new_triangle));
  _mesh_edges.push_back(pair(Edge(e.B(), P), new_triangle));
}

//returns only triangle in mesh with given border edge
Triangle Mesh::find_triangle_with_edge(const Edge &e) const {
  std::optional<int> triangle_index = std::nullopt;
  for (size_t i = 0; i < _mesh_edges.size(); ++i) {
    if (_mesh_edges[i].first == e) {
      triangle_index = i;
    }
  }
  assertm(triangle_index.has_value(), "Edge is not in the mesh!");
  return _mesh_edges[triangle_index.value()].second;
}

// checks Delaunay constraint for triangle T
bool Mesh::check_Delaunay(const Triangle &T) const {
  assertm(T.is_triangle(), "Checking Delaunay of non valid triangle!");
  Point circumcenter = T.get_circumcenter();

  numeric distA = Vector(circumcenter, T.A()).get_length();
  numeric distB = Vector(circumcenter, T.B()).get_length();
  numeric distC = Vector(circumcenter, T.C()).get_length();
  assertm(abs(distA - distB) + abs(distB - distC) + abs(distC - distA) < 10e-3,
          "Wrong circumcenter in Delaunay!");
  numeric dist = distA;

  for (Triangle Tr : _mesh_triangles) {
    Point gravity_center = Tr.get_gravity_center();
    numeric gc_dist = Vector(circumcenter, gravity_center).get_length();

    if (gc_dist < dist) {
      if (!(T.AB() == Tr.AB() || T.AB() == Tr.BC() || T.AB() == Tr.CA() ||
            T.BC() == Tr.AB() || T.BC() == Tr.BC() || T.BC() == Tr.CA() ||
            T.CA() == Tr.AB() || T.CA() == Tr.BC() || T.CA() == Tr.CA())) {
        return false;
      }
    }

    if (Tr.A() != T.A() && Tr.A() != T.B() && Tr.A() != T.C()) {
      numeric dist1 = Vector(circumcenter, Tr.A()).get_length();
      if (dist1 < dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
    if (Tr.B() != T.A() && Tr.B() != T.B() && Tr.B() != T.C()) {
      numeric dist1 = Vector(circumcenter, Tr.B()).get_length();
      if (dist1 < dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
    if (Tr.C() != T.A() && Tr.C() != T.B() && Tr.C() != T.C()) {
      numeric dist1 = Vector(circumcenter, Tr.C()).get_length();
      if (dist1 < dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
  }

  // std::cout << "Delaunay returned TRUE!" << endl;
  return true;
}

// returns vector of points inside Delaunay sphere 
vector<Point> Mesh::get_breakers(Triangle T, const vector<Edge> &active_edges,
                                 const vector<Edge> &checked_edges) const {

  assertm(T.is_triangle(), "Getting breakers of non valid triangle!");
  Point circumcenter = T.get_circumcenter();
  numeric dist = Vector(circumcenter, T.A()).get_length();

  vector<Point> breakers;

  for (Point vertex : _mesh_points) {
    numeric dist1 = Vector(circumcenter, vertex).get_length();
    if (dist1 < 1.1 * dist && vertex != T.A() && vertex != T.B() &&
        vertex != T.C()) {

      bool found = false;

      for (auto edge : active_edges) {
        if (edge.A() == vertex || edge.B() == vertex) {
          breakers.push_back(vertex);
          found = true;
        }
      }
      for (auto edge : checked_edges) {
        if (!found && (edge.A() == vertex || edge.B() == vertex)) {
          breakers.push_back(vertex);
        }
      }
    }
  }
  return breakers;
}

// prints all triangles in mesh to console
void Mesh::output() const {
  for (Triangle Tri : _mesh_triangles) {
    std::cout << Tri.A().x() << " " << Tri.A().y() << " " << Tri.A().z() << " "
              << Tri.B().x() << " " << Tri.B().y() << " " << Tri.B().z() << " "
              << Tri.C().x() << " " << Tri.C().y() << " " << Tri.C().z()
              << endl;
  }
}

// prints number of triangles in mesh to console
void Mesh::cout_triangles_number() const {
  std::cout << "Number of triangles in mesh: " << _mesh_triangles.size()
            << endl;
}

// makes .obj file from mesh triangles
void Mesh::obj_format() const {
  std::ofstream out("out.obj");
  for (size_t i = 0; i < _mesh_triangles.size(); ++i) {
    out << "v " << _mesh_triangles[i].A().x() << " "
        << _mesh_triangles[i].A().y() << " " << _mesh_triangles[i].A().z()
        << endl;
    out << "v " << _mesh_triangles[i].B().x() << " "
        << _mesh_triangles[i].B().y() << " " << _mesh_triangles[i].B().z()
        << endl;
    out << "v " << _mesh_triangles[i].C().x() << " "
        << _mesh_triangles[i].C().y() << " " << _mesh_triangles[i].C().z()
        << endl;
  }
  for (size_t i = 0; i < 3 * _mesh_triangles.size(); i += 3) {
    out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << endl;
  }
}

// returns vector of points closer than 0.4*e_size to point P
// if there are no points returns std::nullopt
std::optional<vector<Point>>
Mesh::empty_surrounding(Point P, numeric e_size,
                        const Edge &working_edge,
                        const vector<Edge> &active_edges,
                        const vector<Edge> &checked_edges) const {
  numeric min_dist = 0.4 * e_size;
  vector<pair<numeric, Point>> close_points;

  for (Point point : _mesh_points) {
    if (point != P && point != working_edge.A() && point != working_edge.B()) {
      numeric dist = Vector(point, P).get_length();
      if (dist < min_dist) {
        bool found = false;
        for (auto edge : active_edges) {
          if (edge.A() == point || edge.B() == point) {
            close_points.push_back(pair(dist, point));
            found = true;
          }
        }
        for (auto edge : checked_edges) {
          if (!found && (edge.A() == point || edge.B() == point)) {
            close_points.push_back(pair(dist, point));
          }
        }
      }
    }
  }

  if (close_points.empty())
    return std::nullopt;

  sort(close_points.begin(), close_points.end(),
       [](auto p1, auto p2) { return p1.first < p2.first; });

  vector<Point> result;

  for (auto dist_point : close_points) {
    result.push_back(dist_point.second);
  }

  return result;
}

// returns true if edge e is in the mesh
bool Mesh::is_in_mesh(const Edge e) const {
  for (auto edge : _mesh_edges) {
    if (edge.first == e)
      return true;
  }
  return false;
}

// splits mesh triangle by point to two triangles
void Mesh::divide_triangle_by_point(const Edge & edge, const Point & P, const Point & new_point){
  assertm(edge.get_midpoint() == P, "Wrong call for divide function!");

  std::optional<int> e_index = std::nullopt;
  for(int i=0; i<_mesh_edges.size(); ++i){
    if(_mesh_edges[i].first == edge){
      assertm(!e_index.has_value(), "Border edge twice in mesh edges!");
      e_index = i;
    }
  } 
  assertm(e_index.has_value(), "No value!");
  std::swap(_mesh_edges[e_index.value()], _mesh_edges.back());
  Triangle T = _mesh_edges.back().second;
  _mesh_edges.pop_back();
  std::optional<int> T_index = std::nullopt;

  for(int i=0; i<_mesh_triangles.size(); ++i){
    if(_mesh_triangles[i] == T){
      assertm(!T_index.has_value(), "Border edge twice in mesh edges!");
      T_index = i;
    }
  }
  assertm(T_index.has_value(), "No value!");
  std::swap(_mesh_triangles[T_index.value()], _mesh_triangles.back());
  assertm(T == _mesh_triangles.back(), "Wrong triangle!");
  _mesh_triangles.pop_back();

  std::optional<Point> other_point = std::nullopt;

  if(T.A() != edge.A() && T.A() != edge.B()){
    other_point = T.A();
  }
  else if(T.B() != edge.A() && T.B() != edge.B()){
    other_point = T.B();
  }
  else if(T.C() != edge.A() && T.C() != edge.B()){
    other_point = T.C();
  }
  assertm(other_point.has_value(), "Point without value!");


  add_triangle(Edge(other_point.value(), edge.A()), new_point);
  add_triangle(Edge(other_point.value(), edge.B()), new_point);

  return;
}