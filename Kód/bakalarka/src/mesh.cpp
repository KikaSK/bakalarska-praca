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
  for (int i = 0; i < _mesh_triangles.size() - 1; ++i) {
    std::cout << "Triangle " << i << ":" << endl;
    std::cout << _mesh_triangles[i] << endl;

    std::cout << "Edge sizes: " << _mesh_triangles[i].AB().get_length() << " "
              << _mesh_triangles[i].BC().get_length() << " "
              << _mesh_triangles[i].CA().get_length() << endl;
  }
  return;
}

void Mesh::add_triangle(Edge e, Point P) {
  Triangle new_triangle(e.A(), e.B(), P);
  assertm(new_triangle.is_triangle(), "Non valid triangle adding to mesh!");
  _mesh_triangles.push_back(new_triangle);
  _mesh_points.push_back(P);
  _mesh_edges.push_back(pair(Edge(P, e.A()), new_triangle));
  _mesh_edges.push_back(pair(Edge(e.B(), P), new_triangle));
}

Triangle Mesh::find_triangle_with_edge(Edge e) const {
  std::optional<int> triangle_index = std::nullopt;
  for (int i = 0; i < _mesh_edges.size(); ++i) {
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

    if (Tr.A() != T.A() && Tr.A() != T.B() && Tr.A() != T.C()) {
      numeric dist1 = Vector(circumcenter, Tr.A()).get_length();
      if (dist1 < 1.1 * dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
    if (Tr.B() != T.A() && Tr.B() != T.B() && Tr.B() != T.C()) {
      numeric dist1 = Vector(circumcenter, Tr.B()).get_length();
      if (dist1 < 1.1 * dist) {
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

vector<Point> Mesh::get_breakers(Triangle T) const {

  assertm(T.is_triangle(), "Getting breakers of non valid triangle!");
  Point circumcenter = T.get_circumcenter();
  numeric dist = Vector(circumcenter, T.A()).get_length();

  vector<Point> breakers;

  for (Point vertex : _mesh_points) {
    numeric dist1 = Vector(circumcenter, vertex).get_length();
    if (dist1 < dist)
      breakers.push_back(vertex);
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
  for (int i = 0; i < _mesh_triangles.size(); ++i) {
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
  for (int i = 0; i < 3 * _mesh_triangles.size(); i += 3) {
    out << "f " << i + 1 << " " << i + 2 << " " << i + 3 << endl;
  }
}

std::optional<vector<Point>> Mesh::empty_surrounding(Point P,
                                                     numeric e_size) const {
  numeric min_dist = 0.3 * e_size;
  vector<pair<numeric, Point>> close_points;

  for (Point point : _mesh_points) {
    if (point != P) {
      numeric dist = Vector(point, P).get_length();
      if (dist < min_dist) {
        close_points.push_back(pair(dist, point));
      }
    }
  }

  if (close_points.empty())
    return std::nullopt;
  // sort(close_points.begin(), close_points.end());

  vector<Point> result;

  for (auto dist_point : close_points) {
    result.push_back(dist_point.second);
  }

  return result;
}
