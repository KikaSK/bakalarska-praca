#include "mesh.h"
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
Mesh::Mesh() = default;

void Mesh::add_triangle(Edge e, Point P) {
  Triangle new_triangle(e.A(), e.B(), P);
  _mesh_triangles.push_back(new_triangle);
  _mesh_points.push_back(P);
  _mesh_edges.push_back(pair(Edge(P, e.A()), new_triangle));
  _mesh_edges.push_back(pair(Edge(e.B(), P), new_triangle));
}

Triangle Mesh::find_triangle_with_edge(Edge e) const {
  for (int i = 0; i < _mesh_edges.size(); ++i) {
    if (_mesh_edges[i].first == e) {
      return _mesh_edges[i].second;
    }
  }
  assertm(false, "Edge is not in the mesh!");
}
