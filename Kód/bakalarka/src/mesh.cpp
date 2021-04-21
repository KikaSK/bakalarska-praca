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
  /*
    std::cout << "Adding new triangle with edge sizes: " << endl
              << "AB: " << e.get_length() << endl
              << "BC: " << Edge(e.B(), P).get_length() << endl
              << "CA: " << Edge(e.A(), P).get_length() << endl;
    std::cout << "Edge is: " << e << endl << "Point is: " << P << endl;
  */
  Triangle new_triangle(e.A(), e.B(), P);

  assertm(new_triangle.is_triangle(), "Non valid triangle adding to mesh!");

  // numeric e_size = 10; // 0.3;

  // assertm(new_triangle.CA().get_length() < 3*e_size, "Weird size of new
  // edge!"); assertm(new_triangle.BC().get_length() < 3*e_size, "Weird size of
  // new edge!"); assertm(new_triangle.AB().get_length() < 3*e_size, "Weird size
  // of new edge!");

  _mesh_triangles.push_back(new_triangle);
  _mesh_points.push_back(P);
  _mesh_edges.push_back(pair(Edge(P, e.A()), new_triangle));
  _mesh_edges.push_back(pair(Edge(e.B(), P), new_triangle));
}

Triangle Mesh::find_triangle_with_edge(const Edge &e) const {
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
    numeric gc_dist = Vector(T.C(), gravity_center).get_length();

    if (gc_dist < 0.4*(T.AB().get_length()+T.BC().get_length()+T.CA().get_length())/3) {
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

vector<Point> Mesh::get_breakers(Triangle T, const vector<Edge> &active_edges,
                                 const vector<Edge> &checked_edges,
                                 const vector<Edge> &bounding_edges) const {

  assertm(T.is_triangle(), "Getting breakers of non valid triangle!");
  Point circumcenter = T.get_circumcenter();
  numeric dist = Vector(circumcenter, T.A()).get_length();

  vector<Point> breakers;

  for (Point vertex : _mesh_points) {
    numeric dist1 = Vector(circumcenter, vertex).get_length();
    if (dist1 < dist && vertex != T.A() && vertex != T.B() &&
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
      /*
      for (auto edge : bounding_edges) {
        if (!found && (edge.A() == vertex || edge.B() == vertex)) {
          breakers.push_back(vertex);
        }
      }
      */
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
                        const vector<Edge> &checked_edges,
                        const vector<Edge> &bounding_edges) const {
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
        for (auto edge : bounding_edges) {
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

void Mesh::divide_triangle_by_point(const Edge &edge, const Point &P,
                                    const Point &new_point) {
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
  assertm(T_index.has_value(), "No value!");
  std::swap(_mesh_triangles[T_index.value()], _mesh_triangles.back());
  assertm(T == _mesh_triangles.back(), "Wrong triangle!");
  _mesh_triangles.pop_back();

  std::optional<Point> other_point = std::nullopt;

  if (T.A() != edge.A() && T.A() != edge.B()) {
    other_point = T.A();
  } else if (T.B() != edge.A() && T.B() != edge.B()) {
    other_point = T.B();
  } else if (T.C() != edge.A() && T.C() != edge.B()) {
    other_point = T.C();
  }
  assertm(other_point.has_value(), "Point without value!");

  add_triangle(Edge(other_point.value(), edge.A()), new_point);
  add_triangle(Edge(other_point.value(), edge.B()), new_point);

  return;
}

bool Mesh::intersections(const Triangle &T_new) const{
  for (auto T : _mesh_triangles){
    if(intersect(T_new, T)) return true;
  }
  return false;
}

bool intersect(const Triangle &T1, const Triangle &T2){
  //if triangles are very far away
  if(Vector(T1.A(), T2.A()).get_length()>5*T1.AB().get_length()) return false;

  Vector n1 = Vector(T1.A(), T1.B())^Vector(T1.A(), T1.C());
  realsymbol x("x"), y("y"), z("z");
  numeric a = n1.x();
  numeric b = n1.y();
  numeric c = n1.z();
  numeric d = -(n1.x()*T1.A().x() + n1.y()*T1.A().y() + n1.z()*T1.A().z());
  ex plane1 = a*x + b*y + c*z + d;
  numeric valA, valB, valC;
  valA = ex_to<numeric> (plane1.subs(
      lst{x == T2.A().x(), y == T2.A().y(), z == T2.A().z()}).evalf());
  valB = ex_to<numeric> (plane1.subs(
      lst{x == T2.B().x(), y == T2.B().y(), z == T2.B().z()}).evalf());
  valC = ex_to<numeric> (plane1.subs(
      lst{x == T2.C().x(), y == T2.C().y(), z == T2.C().z()}).evalf());

  if(valA<-10e-8 && valB<-10e-8 && valC<-10e-8) return false;
  if(valA>10e-8 && valB>10e-8 && valC>10e-8) return false;

  numeric AB = valA*valB;
  numeric BC = valB*valC;
  numeric AC = valA*valC;

  // plane
  // ax + by + cz + d = 0

  // line
  // x = p_x + u_x*t
  // y = p_y + u_y*t
  // z = p_z + u_z*t

  // t = -(a*p_x + b*p_y + c*p_z + d)/(a*u_x + b*u_y + c*u_z)

  std::optional<Point> intersectionAB = std::nullopt, intersectionBC = std::nullopt, intersectionCA = std::nullopt;

  // ak A a B sú na opačnej strane, nájdeme priesečník s rovinou
  if(AB<-10e-8){
    Point p = T2.A();
    Vector u = Vector(T2.A(), T2.B());
    numeric t = -(a*p.x() + b*p.y() + c*p.z() + d)/(a*u.x() + b*u.y() + c*u.z());
    intersectionAB = Point(p, t*u);
  }
  // ak B a C sú na opačnej strane, nájdeme priesečník s rovinou
  if(BC<-10e-8){
    Point p = T2.B();
    Vector u = Vector(T2.B(), T2.C());
    numeric t = -(a*p.x() + b*p.y() + c*p.z() + d)/(a*u.x() + b*u.y() + c*u.z());
    intersectionBC = Point(p, t*u);
  }
  // ak A a C sú na opačnej strane, nájdeme priesečník s rovinou
  if(AC<-10e-8){
    Point p = T2.A();
    Vector u = Vector(T2.A(), T2.C());
    numeric t = -(a*p.x() + b*p.y() + c*p.z() + d)/(a*u.x() + b*u.y() + c*u.z());
    intersectionCA = Point(p, t*u);
  }
  std::optional<Point> P1 = std::nullopt, P2 = std::nullopt;
  if(intersectionAB.has_value()){
    P1 = intersectionAB.value();
    if(intersectionBC.has_value())
      P2 = intersectionBC.value();
    else if(intersectionCA.has_value())
      P2 = intersectionCA.value();
    else
      assertm(false, "Not enough intersection points!");
  }
  else if(intersectionBC.has_value() && intersectionCA.has_value()){
    P1 = intersectionBC.value();
    P2 = intersectionCA.value();
  }
  else{ 
    assertm(false, "Not enough intersection points!");
  }

  assertm(!(intersectionAB.has_value() && intersectionBC.has_value() && intersectionCA.has_value()), 
  "Too much intersection points!");
  assertm(P1.has_value() && P2.has_value(), "Not enough intersection points!");

  if(T1.is_in_triangle(P1.value()) || T1.is_in_triangle(P2.value()))
    return true;

  //mame trojuholnik T v rovine a body P1 a P2 v tej istej rovine ktore nelezia v trojuholniku
  //chceme zistit ci usecka P1, P2 pretína trojuholník T

  Vector vAB = Vector(T1.A(), T1.B());
  Vector vBC = Vector(T1.B(), T1.C());
  Vector vCA = Vector(T1.C(), T1.A());
  Vector other = Vector(P1.value(), P2.value());

  // usecky sa pretinaju ak sa pretina usecka1 s priamkou2 (prva cast) a usecka2 s priamkou1 (druha cast)
  if
  (
    //pretinaju sa priamka AB a usecka P1P2
    ((Vector(T1.A(), P1.value())^vAB) * (Vector(T1.A(), P2.value())^vAB)) < 0
  &&
    //pretinaju sa priamka P1P1 a usecka AB
    ((Vector(P1.value(), T1.A())^other) * (Vector(P1.value(), T1.B())^other)) < 0
  ) 
    return true;

  // to iste pre usecku BC
  if
  (
    //pretinaju sa priamka BC a usecka P1P2
    ((Vector(T1.B(), P1.value())^vBC) * (Vector(T1.B(), P2.value())^vBC)) < 0
  &&
    //pretinaju sa priamka P1P1 a usecka BC
    ((Vector(P1.value(), T1.B())^other) * (Vector(P1.value(), T1.C())^other)) < 0
  )
    return true;

  // to iste pre usecku CA
  if
  (
    //pretinaju sa priamka CA a usecka P1P2
    ((Vector(T1.C(), P1.value())^vCA) * (Vector(T1.C(), P2.value())^vCA)) < 0
  &&
    //pretinaju sa priamka P1P1 a usecka CA
    ((Vector(P1.value(), T1.C())^other) * (Vector(P1.value(), T1.A())^other)) < 0
  )
    return true;

  // nepretinaju sa
  return false;
}