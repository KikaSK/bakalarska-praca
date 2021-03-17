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

  /*
  _mesh_active_edges.push_back( pair(Edge(T.A(), T.B()), T) );
  _mesh_active_edges.push_back( pair(Edge(T.B(), T.C()), T) );
  _mesh_active_edges.push_back( pair(Edge(T.C(), T.A()), T) );
  */
}

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

void Mesh::add_triangle(Edge e, Point P) {

  std::cout << "Adding new triangle with edge sizes: " << endl
            << "AB: " << e.get_length() << endl
            << "BC: " << Edge(e.B(), P).get_length() << endl
            << "CA: " << Edge(e.A(), P).get_length() << endl;
  std::cout << "Edge is: " << e << endl << "Point is: " << P << endl;
  numeric e_size = 10;//0.3;
  Triangle new_triangle(e.A(), e.B(), P);

  assertm(new_triangle.is_triangle(), "Non valid triangle adding to mesh!");

  // assertm(new_triangle.CA().get_length() < 3*e_size, "Weird size of new
  // edge!"); assertm(new_triangle.BC().get_length() < 3*e_size, "Weird size of
  // new edge!"); assertm(new_triangle.AB().get_length() < 3*e_size, "Weird size
  // of new edge!");

  _mesh_triangles.push_back(new_triangle);
  _mesh_points.push_back(P);
  _mesh_edges.push_back(pair(Edge(P, e.A()), new_triangle));
  _mesh_edges.push_back(pair(Edge(e.B(), P), new_triangle));
}

Triangle Mesh::find_triangle_with_edge(const Edge & e) const {
  std::optional<int> triangle_index = std::nullopt;
  for (size_t i = 0; i < _mesh_edges.size(); ++i) {
    if (_mesh_edges[i].first == e) {
      // assertm(!triangle_index.has_value(), "More than one triangle containing
      // edge!");
      triangle_index = i;
    }
  }
  assertm(triangle_index.has_value(), "Edge is not in the mesh!");
  return _mesh_edges[triangle_index.value()].second;
}

bool Mesh::check_Delaunay(Triangle T) const {
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
      if (dist1 < 1.1 * dist) {
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
      if (dist1 < 1.1 * dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
  }

  // std::cout << "Delaunay returned TRUE!" << endl;
  return true;
}

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

void Mesh::output() const {
  for (Triangle Tri : _mesh_triangles) {
    std::cout << Tri.A().x() << " " << Tri.A().y() << " " << Tri.A().z() << " "
              << Tri.B().x() << " " << Tri.B().y() << " " << Tri.B().z() << " "
              << Tri.C().x() << " " << Tri.C().y() << " " << Tri.C().z()
              << endl;
  }
}

void Mesh::cout_triangles_number() const {
  std::cout << "Number of triangles in mesh: " << _mesh_triangles.size()
            << endl;
}
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

std::optional<vector<Point>>
Mesh::empty_surrounding(Point P, numeric e_size,
                        const vector<Edge> &active_edges,
                        const vector<Edge> &checked_edges) const {
  numeric min_dist = 0.4 * e_size;
  vector<pair<numeric, Point>> close_points;

  for (Point point : _mesh_points) {
    if (point != P) {
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

bool Mesh::is_in_mesh(const Edge e) const {
  for (auto edge : _mesh_edges) {
    if (edge.first == e)
      return true;
  }
  return false;
}
