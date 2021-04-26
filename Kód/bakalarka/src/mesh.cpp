#include "mesh.h"

#include <fstream>
#include <string>

Mesh::Mesh(Triangle T) {
  _mesh_triangles.push_back(T);
  _mesh_points.push_back(T.A());
  _mesh_points.push_back(T.B());
  _mesh_points.push_back(T.C());
  _mesh_edges.push_back(pair(Edge(T.A(), T.B()), T));
  _mesh_edges.push_back(pair(Edge(T.B(), T.C()), T));
  _mesh_edges.push_back(pair(Edge(T.C(), T.A()), T));
}

// prints triangles to console
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
  //_mesh_edges.push_back(pair(e, new_triangle));
}

// returns only triangle in mesh with given border edge
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
bool Mesh::check_Delaunay(const Triangle &T, const Edge &working_edge, const Triangle &neighbour_triangle) const {
  std::optional<Point> other_point = std::nullopt;
  if(neighbour_triangle.AB() == working_edge){
    other_point = neighbour_triangle.C();
  }
  else if(neighbour_triangle.BC() == working_edge){
    other_point = neighbour_triangle.A();
  }
  else if(neighbour_triangle.CA() == working_edge){
    other_point = neighbour_triangle.B();
  }

  assertm(other_point.has_value(), "Didn't find third point in neighbour triangle!");
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

    if (Tr.A() != T.A() && Tr.A() != T.B() && Tr.A() != T.C() && Tr.A() != other_point.value()) {
      numeric dist1 = Vector(circumcenter, Tr.A()).get_length();
      if (dist1 < dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
    if (Tr.B() != T.A() && Tr.B() != T.B() && Tr.B() != T.C() && Tr.B() != other_point.value()) {
      numeric dist1 = Vector(circumcenter, Tr.B()).get_length();
      if (dist1 < dist) {
        // std::cout << "Delaunay returned FALSE!" << endl;
        return false;
      }
    }
    if (Tr.C() != T.A() && Tr.C() != T.B() && Tr.C() != T.C() && Tr.C() != other_point.value()) {
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
void Mesh::obj_format(const std::string &name) const {
  std::ofstream out(name + ".obj");
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

// returns true if edge e is in the mesh
bool Mesh::is_in_mesh(const Edge e) const {
  for (auto edge : _mesh_edges) {
    if (edge.first == e)
      return true;
  }
  return false;
}

// splits mesh triangle by point to two triangles
void Mesh::divide_triangle_by_point(const Edge &edge, const Point &P) {
  assertm(edge.get_midpoint() == P, "Wrong call for divide function!");

  std::optional<int> e_index = std::nullopt;
  for (int i = 0; i < _mesh_edges.size(); ++i) {
    if (_mesh_edges[i].first == edge) {
      assertm(!e_index.has_value(), "Border edge twice in mesh edges!");
      e_index = i;
    }
  }
  assertm(e_index.has_value(), "No value!");
  std::swap(_mesh_edges[e_index.value()], _mesh_edges.back());
  Triangle T = _mesh_edges.back().second;
  _mesh_edges.pop_back();
  std::optional<int> T_index = std::nullopt;

  for (int i = 0; i < _mesh_triangles.size(); ++i) {
    if (_mesh_triangles[i] == T) {
      assertm(!T_index.has_value(), "Border edge twice in mesh edges!");
      T_index = i;
    }
  }
  assertm(T_index.has_value(), "Triangle not in mesh triangles!");
  std::swap(_mesh_triangles[T_index.value()], _mesh_triangles.back());
  assertm(T == _mesh_triangles.back(), "Wrong triangle!");
  _mesh_triangles.pop_back();

  std::optional<Point> other_point = std::nullopt;

  if (T.BC() == edge) {
    other_point = T.A();
  } else if (T.CA() == edge) {
    other_point = T.B();
  } else if (T.AB() == edge) {
    other_point = T.C();
  }
  assertm(other_point.has_value(), "Point without value!");

  add_triangle(Edge(other_point.value(), edge.A()), P);
  add_triangle(Edge(other_point.value(), edge.B()), P);

  return;
}
std::optional<Triangle> Mesh::find_neighbour_triangle(const Edge &e, const Triangle &T) const{

  std::optional<Triangle> triangle = std::nullopt;
  for (auto MT : _mesh_triangles) {
    if (
      (MT.AB() == e || MT.BC() == e || MT.CA() == e)
    && 
    !(MT == T)
    ) 
    {
      triangle = MT;
    }
  }
  return triangle;
}
void Mesh::adaptive(const numeric &precision, const Function &F, const numeric &e_size){
  std::cout<<"Adapting.."<<endl;
  vector<int>indexes;
  vector<Triangle>new_mesh_triangles;
  for (int i = 0; i < _mesh_triangles.size(); ++i){
    Triangle T = _mesh_triangles[i];
    Point gravity_center = T.get_gravity_center();
    Vector direction = F.get_gradient_at_point(gravity_center).unit();
    Point projected_gc = project(gravity_center, direction, F, e_size);
    numeric dist = Vector(gravity_center, projected_gc).get_length();
    if(dist < precision){
      new_mesh_triangles.push_back(T);
    }
    else{
      Edge e1 = T.AB();
      Edge e2 = T.BC();
      Edge e3 = T.CA();

      std::optional<Triangle> e1T = find_neighbour_triangle(e1, T);
      Point P1 = e1.get_midpoint();
      if(e1T.has_value()){
        Point gc_e1T = e1T.value().get_gravity_center();
        Vector dir_e1T = F.get_gradient_at_point(gc_e1T).unit();
        Point projected_gc_e1T = project(gc_e1T, dir_e1T, F, e_size);
        numeric dist_e1T = Vector(gc_e1T, projected_gc_e1T).get_length();
        if(dist_e1T >= precision)
        {
          P1 = project(e1.get_midpoint(), F.get_gradient_at_point(e1.get_midpoint()), F, e_size);
        }
      }
      /*
      else{
        P1 = project(T.AB().get_midpoint(), F.get_gradient_at_point(T.AB().get_midpoint()), F, e_size);
      }
      */
      std::optional<Triangle> e2T = find_neighbour_triangle(e2, T);
      Point P2 = e2.get_midpoint();
      if(e2T.has_value()){
        Point gc_e2T = e2T.value().get_gravity_center();
        Vector dir_e2T = F.get_gradient_at_point(gc_e2T).unit();
        Point projected_gc_e2T = project(gc_e2T, dir_e2T, F, e_size);
        numeric dist_e2T = Vector(gc_e2T, projected_gc_e2T).get_length();
        if(dist_e2T >= precision)
        {
          P2 = project(e2.get_midpoint(), F.get_gradient_at_point(e2.get_midpoint()), F, e_size);
        }
      }
      /*
      else{
        P2 = project(e2.get_midpoint(), F.get_gradient_at_point(e2.get_midpoint()), F, e_size);
      }
      */
      std::optional<Triangle> e3T = find_neighbour_triangle(e3, T);
      Point P3 = e3.get_midpoint();
      if(e3T.has_value()){
        Point gc_e3T = e3T.value().get_gravity_center();
        Vector dir_e3T = F.get_gradient_at_point(gc_e3T).unit();
        Point projected_gc_e3T = project(gc_e3T, dir_e3T, F, e_size);
        numeric dist_e3T = Vector(gc_e3T, projected_gc_e3T).get_length();
        if(dist_e3T >= precision)
        {
          P3 = project(e3.get_midpoint(), F.get_gradient_at_point(e3.get_midpoint()), F, e_size);
        }
      }
      /*
      else{
        P3 = project(e3.get_midpoint(), F.get_gradient_at_point(e3.get_midpoint()), F, e_size);
      }
      */
      
      Triangle T1(T.A(), P1, P3);
      Triangle T2(T.C(), P2, P3);
      Triangle T3(T.B(), P1, P2);
      Triangle T4(P1, P2, P3);
      new_mesh_triangles.push_back(T1);
      new_mesh_triangles.push_back(T2);
      new_mesh_triangles.push_back(T3);
      new_mesh_triangles.push_back(T4);
    }
  }
  _mesh_triangles.empty();
  _mesh_triangles = new_mesh_triangles;
  return;
}