#ifndef BASIC_ALGORITHM_H
#define BASIC_ALGORITHM_H

#include "function.h"
#include "mesh.h"
#include "triangle.h"

#include <vector>
class BoundingBox{
public:
  BoundingBox(numeric _min_x, numeric _max_x, numeric _min_y, numeric _max_y, numeric _min_z, numeric _max_z):
  _min_x(_min_x), _max_x(_max_x), _min_y(_min_y), _max_y(_max_y), _min_z(_min_z), _max_z(_max_z){};
  
  numeric min_x() const {return _min_x;}
  numeric max_x() const {return _max_x;}
  numeric min_y() const {return _min_y;}
  numeric max_y() const {return _max_y;}
  numeric min_z() const {return _min_z;}
  numeric max_z() const {return _max_z;}

  bool in_interval_x(const numeric x) const {return (x >= _min_x && x <=_max_x);}
  bool in_interval_y(const numeric y) const {return (y >= _min_y && y <=_max_y);}
  bool in_interval_z(const numeric z) const {return (z >= _min_z && z <=_max_z);}

  bool is_inside(const Point P) const {
    return in_interval_x(P.x()) && in_interval_y(P.y()) && in_interval_z(P.z());
  }
  bool is_on(const Point P) const {
    bool x_wall = ((P.x()-_min_x < 10e-10 || P.x() - _max_x < 10e-10) && in_interval_y(P.y()) && in_interval_z(P.z()));
    bool y_wall = ((P.y()-_min_y < 10e-10 || P.y() - _max_y < 10e-10) && in_interval_x(P.x()) && in_interval_z(P.z()));
    bool z_wall = ((P.z()-_min_z < 10e-10 || P.z() - _max_z < 10e-10) && in_interval_y(P.y()) && in_interval_x(P.x()));
    return (x_wall || y_wall || z_wall);
  }

private:
  numeric _min_x;
  numeric _max_x;
  numeric _min_y;
  numeric _max_y;
  numeric _min_z;
  numeric _max_z;
};

class BasicAlgorithm {
public:
  BasicAlgorithm(Function f, Triangle seed_triangle, numeric e_size,
                 realsymbol x, realsymbol y, realsymbol z, BoundingBox b)
      : F(f), active_edges(), checked_edges(), my_mesh(seed_triangle),
        e_size(e_size), x(x), y(y), z(z), bounding_box(b){
    active_edges.push_back(seed_triangle.AB());
    active_edges.push_back(seed_triangle.BC());
    active_edges.push_back(seed_triangle.CA());
  }

  Mesh calculate();

private:
  Function F;
  vector<Edge> active_edges;
  vector<Edge> checked_edges;
  Mesh my_mesh;
  numeric e_size;
  realsymbol x;
  realsymbol y;
  realsymbol z;
  BoundingBox bounding_box;
};



#include "algorithms.h"
#include "assertm.h"
#include <iostream>

using std::cout;

// returns true if prev/next triangle is added to mesh else returns false
bool fix_prev_next(Mesh &my_mesh, vector<Edge> &active_edges,
                   vector<Edge> &checked_edges, const Edge &working_edge,
                   const bool is_prev, const numeric e_size) {

  // find previous and next points
  auto [prev, next] =
      find_prev_next(my_mesh, working_edge, active_edges, checked_edges);
  assertm(Vector(working_edge.A(), prev).get_length() < 3 * e_size,
          "Wrong prev point!");
  assertm(Vector(working_edge.B(), next).get_length() < 3 * e_size,
          "Wrong prev point!");
  assertm(is_border(Edge(working_edge.A(), prev), active_edges, checked_edges),
          "Prev edge not in border!");
  assertm(is_border(Edge(working_edge.B(), next), active_edges, checked_edges),
          "Next edge not in border!");

  Point vertex = prev;
  Edge edge = working_edge;

  // fix so that vertex is the new vertex and edge.A() is the adjacent vertex
  if (!is_prev) {
    edge = Edge(working_edge.B(), working_edge.A());
    vertex = next;
  }

  // potentialy new triangle
  Triangle maybe_new_T(edge.A(), edge.B(), vertex);

  // already existing neighbour triangle of working edge
  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(edge);

  if(maybe_new_T.AB().get_length()>2*e_size || maybe_new_T.CA().get_length()>2*e_size || maybe_new_T.BC().get_length()>2*e_size)
  {
    return false;
  }

  // checks if the potential triangle has good orientation and angle near A is
  // less than 90 degrees and checks Delaunay
  if (maybe_new_T.is_triangle() &&
      good_orientation(edge, vertex, neighbour_triangle) &&
      my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, edge, vertex)) {

    // cout << "Found prev triangle!" << endl;

    /*assertm(Vector(working_edge.A(), vertex).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev or next point!");
    assertm(Vector(working_edge.B(), vertex).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev or next point!");
            */

    my_mesh.add_triangle(edge, vertex);

    assertm(is_border(Edge(edge.A(), vertex), active_edges, checked_edges),
            "Neighbour edge not in border!");
    /*assertm(!(is_border(Edge(edge.B(), vertex), active_edges, checked_edges)) ||
                (prev == next),
            "New edges already in border!");*/

    Edge new_edge1(edge.A(), vertex);
    Edge new_edge2(edge.B(), vertex);
    update_border(new_edge1, new_edge2, active_edges, checked_edges);
    
    return true;
  }
  return false;
}


// returns true if overlap triangle is added to mesh else returns false
bool fix_overlap(Mesh &my_mesh, const Edge &working_edge,
                 vector<Edge> &active_edges, vector<Edge> &checked_edges,
                 Point overlap_point, const Function &F) {
  // assertm(false, "In check overlap!");

  auto [prev, next] =
      find_prev_next(my_mesh, working_edge, active_edges, checked_edges);

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  assertm(neighbour_triangle.is_triangle(), "Neighbour triangle not valid!");

  // checks if overlap point is not neighbour or working edge point also if
  // overlap is on the border and if overlap triangle has good orientation of
  // normal
  if (is_vertex_good_possibility(overlap_point, prev, next, working_edge,
                                 neighbour_triangle, active_edges,
                                 checked_edges, my_mesh, F))

  {
    Triangle maybe_new_T =
        Triangle(working_edge.A(), working_edge.B(), overlap_point);

    // if Delaunay constraint is satisfied add the triangle to
    // triangulation and end
    if (maybe_new_T.is_triangle() &&
        good_orientation(working_edge, overlap_point, neighbour_triangle) &&
        my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, working_edge, overlap_point)) {

      // assertm(false, "Found overlap triangle, hooray!");
      // cout << "Found overlap triangle" << endl;

      assertm(Vector(working_edge.A(), overlap_point).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of overlap point!");
      assertm(Vector(working_edge.B(), overlap_point).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of overlap point!");
      my_mesh.add_triangle(working_edge, overlap_point);

      Edge new_edge1 = Edge(working_edge.A(), overlap_point);
      Edge new_edge2 = Edge(working_edge.B(), overlap_point);

      assertm(new_edge1 != new_edge2, "Same edges!");
      update_border(new_edge1, new_edge2, active_edges, checked_edges);

      return true;
    }
  }
  return false;
}

bool fix_proj(Mesh &my_mesh, vector<Edge> &active_edges,
              vector<Edge> &checked_edges, const Edge &working_edge,
              numeric e_size, const Function &F, const BoundingBox & bounding_box) {

  Point projected = get_projected(working_edge, active_edges, checked_edges, e_size, my_mesh, F);

  Triangle maybe_new_T(working_edge.A(), working_edge.B(), projected);
  assertm(maybe_new_T.is_triangle(), "Proj triangle not valid!");

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);

  // checks if there are some points very close to projected point
  if (auto surrounding_points = my_mesh.empty_surrounding(projected, e_size, active_edges, checked_edges)
          ;surrounding_points.has_value()) {

    // points closer to projected point than 0.3*e_size sorted from closest
    vector<Point> close_points =
        surrounding_points
            .value();

    for (auto close_point : close_points) {

      Triangle maybe_new_T(working_edge.A(), working_edge.B(), close_point);

      if (maybe_new_T.is_triangle() &&
          good_orientation(working_edge, close_point, neighbour_triangle) &&
          my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, working_edge, close_point)) {

        auto [prev, next] =
            find_prev_next(my_mesh, working_edge, active_edges, checked_edges);

        // if close point is prev we want to try fix prev
        if (close_point == prev) {
          if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                            true, e_size))
            {
            //cout<<"Fix prev!"<<endl;
            return true;
            }
        }
        // if close point is next we want to try fix next
        else if (close_point == next) {
          if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                            false, e_size))
            {
            //cout<<"Fix next!"<<endl;
            return true;
            }
        }

        // if close point is overlap we want to try fix overlap
        else {
          if (fix_overlap(my_mesh, working_edge, active_edges, checked_edges,
                          close_point, F))
            {
            //cout<<"Fix overlap!"<<endl;
            return true;
            }
        }
      }
    }
    // if there are close points but nothing worked we want to try to
    // construct original triangle
    // cout << "There are close points but nothing worked!" << endl;
  }
  auto close_edge = get_closest_edge(active_edges, checked_edges, projected, neighbour_triangle).value();
  if(close_edge.second < e_size/3){
    Edge closest_edge = close_edge.first;
    Point P1 = closest_edge.A();
    Point P2 = closest_edge.B();
    Vector n_A = F.get_gradient_at_point(closest_edge.A()).unit();
    Vector n_B = F.get_gradient_at_point(closest_edge.B()).unit();
    Vector normal = (n_A + n_B) / 2;
    Point P3 = project(closest_edge.get_midpoint(), normal, F, e_size);
    
    if(P1!= working_edge.A() && P1!=working_edge.B()){
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P1);
      if(maybe_new_T.is_triangle() &&
      good_orientation(working_edge, P1, neighbour_triangle) &&
      my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, working_edge, P1)) {
        my_mesh.add_triangle(working_edge, P1);
        Edge new_edge1(working_edge.A(), P1);
        Edge new_edge2(working_edge.B(), P1);
        update_border(new_edge1, new_edge2, active_edges, checked_edges);
        return true;
      }
    }

    if(P2!= working_edge.A() && P2!=working_edge.B()){
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P2);
      if(maybe_new_T.is_triangle() &&
      good_orientation(working_edge, P2, neighbour_triangle) &&
      my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, working_edge, P2)) {
        my_mesh.add_triangle(working_edge, P2);
        Edge new_edge1(working_edge.A(), P2);
        Edge new_edge2(working_edge.B(), P2);
        update_border(new_edge1, new_edge2, active_edges, checked_edges);
        return true;
      }
    }
    if(P3!= working_edge.A() && P3!=working_edge.B()){
    
      Triangle maybe_new_T(working_edge.A(), working_edge.B(), P3);
      if(maybe_new_T.is_triangle() &&
      good_orientation(working_edge, P3, neighbour_triangle) &&
      my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, working_edge, P3)) {
        my_mesh.add_triangle(working_edge, P3);
        Edge new_edge1(working_edge.A(), P3);
        Edge new_edge2(working_edge.B(), P3);
        my_mesh.divide_triangle_by_point(closest_edge, closest_edge.get_midpoint(), P3);
        delete_from_active(closest_edge, active_edges);
        delete_from_checked(closest_edge, checked_edges);
        push_edge_to_active(Edge(closest_edge.A(), P3), active_edges);
        push_edge_to_active(Edge(closest_edge.B(), P3), active_edges);
        
        update_border(new_edge1, new_edge2, active_edges, checked_edges);
        return true;
      }
    }
  }

  if (my_mesh.check_Delaunay(maybe_new_T) && good_edges(my_mesh, active_edges, checked_edges, working_edge, projected) && bounding_box.is_inside(projected)) {

    // cout << "Found new triangle" << endl;

    assertm(Vector(working_edge.A(), projected).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of projected point!");
    assertm(Vector(working_edge.B(), projected).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of projected point!");

    my_mesh.add_triangle(working_edge, projected);

    Edge new_edge1(working_edge.A(), projected);
    Edge new_edge2(working_edge.B(), projected);
    assertm(new_edge1 != new_edge2, "Same edges!");
    update_border(new_edge1, new_edge2, active_edges, checked_edges);

    return true;
  }
  // if nothing worked we move on to another step
  // push_edge_to_checked(working_edge, checked_edges);
  return false;
}

// one step of the algorithm
bool step(Mesh &my_mesh, vector<Edge> &active_edges,
          vector<Edge> &checked_edges, const Edge &working_edge,
          const numeric e_size, const Function &F, const BoundingBox & bounding_box) {

  my_mesh.obj_format();

  assertm(!is_border(working_edge, active_edges, checked_edges),
          "Workind edge found in border edges!");

  if (fix_proj(my_mesh, active_edges, checked_edges, working_edge, e_size, F, bounding_box)) {
    return true;
  } else if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                           true, e_size))
    return true;
  else if (fix_prev_next(my_mesh, active_edges, checked_edges, working_edge,
                         false, e_size))
    return true;
  else {
    Point projected = get_projected(working_edge, active_edges, checked_edges, e_size, my_mesh, F);
    Triangle proj_T(working_edge.A(), working_edge.B(), projected);
    Triangle neighbour_T = my_mesh.find_triangle_with_edge(working_edge);
    vector<Point> breakers =
        my_mesh.get_breakers(proj_T, active_edges, checked_edges);
    std::sort(breakers.begin(), breakers.end(),
              [&working_edge, &neighbour_T](auto i, auto j) {
                return line_point_dist(working_edge, i, neighbour_T) <
                       line_point_dist(working_edge, j, neighbour_T);
              });

    for (auto point : breakers) {
      if (good_orientation(working_edge, point, neighbour_T)) {
        if (is_border_point(point, active_edges, checked_edges)) {
          // assertm(Vector(point, working_edge.A()).get_length() < 3 * e_size,
          //        "Too big distance of break point!");
          if (fix_overlap(my_mesh, working_edge, active_edges, checked_edges,
                          point, F))
            return true;
        }
      }
    }

    // cout << "No triangle found" << endl;
    // if no triangle was created
    assertm(!is_border(working_edge, active_edges, checked_edges),
            "Something very wrong!");
    push_edge_to_checked(working_edge, checked_edges);
    return false;
    // something is wtong with find_prev_next
    /*
        auto [prev, next] =
            find_prev_next(my_mesh, working_edge, active_edges, checked_edges);

        Triangle neighbour_T = my_mesh.find_triangle_with_edge(working_edge);

        if (prev == next) {
          Triangle new_T(working_edge.A(), working_edge.B(), prev);
          if (new_T.is_triangle() && !(neighbour_T == new_T) &&
              good_orientation(working_edge, prev, neighbour_T)) {
            my_mesh.add_triangle(working_edge, prev);
            delete_from_active(Edge(working_edge.A(), prev), active_edges);
            delete_from_active(Edge(working_edge.B(), prev), active_edges);
            delete_from_checked(Edge(working_edge.A(), prev), checked_edges);
            delete_from_checked(Edge(working_edge.B(), prev), checked_edges);
            return;
          }
        } else {
          std::optional<Point> min_dist_point = closest_point(working_edge,
       breakers); if(!min_dist_point.has_value()){ Triangle
       new_prev_T(working_edge.A(), working_edge.B(), prev); Triangle
       new_next_T(working_edge.A(), working_edge.B(), next);
            if(new_prev_T.is_triangle() && !(neighbour_T == new_prev_T) &&
              good_orientation(working_edge, prev, neighbour_T)){
                delete_from_active(Edge(working_edge.A(), prev), active_edges);
                delete_from_checked(Edge(working_edge.A(), prev),
       checked_edges); push_edge_to_active(Edge(working_edge.B(), prev),
       active_edges);
              }
            else if(new_next_T.is_triangle() && !(neighbour_T == new_next_T) &&
              good_orientation(working_edge, next, neighbour_T)){
                delete_from_active(Edge(working_edge.B(), next), active_edges);
                delete_from_checked(Edge(working_edge.B(), next),
       checked_edges); push_edge_to_active(Edge(working_edge.A(), next),
       active_edges);
            }
            else{
              push_edge_to_checked(working_edge, checked_edges);
            }
            return;
          }
          assertm(min_dist_point.value() != working_edge.A() &&
                      min_dist_point.value() != working_edge.B(),
                  "Closest_point function is wrong!");

          Triangle closest_point_T(working_edge.A(), working_edge.B(),
                                   min_dist_point.value());
          if (closest_point_T.is_triangle() &&
              good_orientation(working_edge, min_dist_point.value(),
       neighbour_T)) { my_mesh.add_triangle(working_edge,
       min_dist_point.value());

            Edge new_edge1(working_edge.A(), min_dist_point.value());
            Edge new_edge2(working_edge.B(), min_dist_point.value());

            assertm(new_edge1 != new_edge2, "Same edges!");
            assertm(prev != next, "Prev should not be next!");

            if(is_border(new_edge1, active_edges, checked_edges)){
              delete_from_active(new_edge1, active_edges);
              delete_from_checked(new_edge1, checked_edges);
            }
            else{
              push_edge_to_active(new_edge1, active_edges);
            }

            if(is_border(new_edge2, active_edges, checked_edges)){
              delete_from_active(new_edge2, active_edges);
              delete_from_checked(new_edge2, checked_edges);
            }
            else{
              push_edge_to_active(new_edge2, active_edges);
            }
            return;
          }
          else {*/
    //}
    //}
  }
  assertm(false, "Should not get there!");
  push_edge_to_checked(working_edge, checked_edges);
  return false;
}

//just help in determining all border edges in .obj viewer
void add_marks(Mesh &my_mesh, const vector<Edge> &active_edges,
               const vector<Edge> &checked_edges, const numeric & e_size) {
  auto border = connect_edges(active_edges, checked_edges);
  for (auto edge : border) {
    auto dir =
        e_size/5 * find_direction(edge, my_mesh.find_triangle_with_edge(edge),
                                    e_size);
    Vector edge_dir = e_size/5 * Vector(edge.A(), edge.B()).unit();
    Edge new_e(Point(edge.get_midpoint(), edge_dir.vector_inverse() / 2),
               Point(edge.get_midpoint(), edge_dir / 2));
    auto new_p = Point(edge.get_midpoint(), dir);
    my_mesh.add_triangle(new_e, new_p);
  }
}

int fix_holes(Mesh &my_mesh, const Function &F, const Edge &working_edge,
              vector<Edge> &active_edges, vector<Edge> &checked_edges,
              const numeric e_size) {

  int number_of_new_edges = 0;

  cout << "In fix_holes!" << endl;
  checked_edges.push_back(working_edge); add_marks(my_mesh, active_edges, checked_edges, e_size); my_mesh.obj_format(); return 0;

  assertm(!is_border(working_edge, active_edges, checked_edges),
          "Working edge found in border!");
  Point projected = get_projected(working_edge, active_edges, checked_edges, e_size, my_mesh, F);

  Triangle maybe_new_T(working_edge.A(), working_edge.B(), projected);
  assertm(maybe_new_T.is_triangle(), "Proj triangle not valid!");

  Triangle neighbour_triangle = my_mesh.find_triangle_with_edge(working_edge);
  auto breakers =
      my_mesh.get_breakers(maybe_new_T, active_edges, checked_edges);
  auto close_points =
      my_mesh.empty_surrounding(projected, e_size, active_edges, checked_edges);
  if (close_points.has_value()) {
    if (breakers.empty()) {
      breakers = close_points.value();
    }
    breakers = connect_points(breakers, close_points.value());
    /*breakers.insert(breakers.end(), close_points.value().begin(),
                    close_points.value().end());*/
  }
  //if (breakers.empty()) {
    //if (fix_proj(my_mesh, active_edges, checked_edges, working_edge, e_size, F))
      //return 2;
    //else {
      /*std::sort(breakers.begin(), breakers.end(), [&working_edge,
    &neighbour_triangle](auto i, auto j){ return line_point_dist(working_edge,
    i, neighbour_triangle) < line_point_dist(working_edge, j,
    neighbour_triangle);
    });*/
      std::optional<Point> maybe_new_point = // breakers[0];
          get_closest_point(my_mesh, active_edges, checked_edges, working_edge,
                            neighbour_triangle, e_size);

      if (maybe_new_point.has_value()) {
        Point closest_point = maybe_new_point.value();
        Triangle new_triangle(working_edge.A(), working_edge.B(),
                              closest_point);
        my_mesh.add_triangle(working_edge, closest_point);
        Edge new_edge1(working_edge.A(), closest_point);
        Edge new_edge2(working_edge.B(), closest_point);
        number_of_new_edges =
            update_border(new_edge1, new_edge2, active_edges, checked_edges);
        return number_of_new_edges;
      }
    push_edge_to_checked(working_edge, checked_edges);
    return false;

  std::optional<pair<Point, numeric>> closest_point = std::nullopt;

  for (auto point : breakers) {
    if (!closest_point.has_value()) {
      closest_point =
          pair(point, line_point_dist(working_edge, point, neighbour_triangle));
    } else if (Triangle(working_edge.A(), working_edge.B(), point)
                   .is_triangle() &&
               line_point_dist(working_edge, point, neighbour_triangle) <
                   closest_point.value().second &&
               good_orientation(working_edge, point, neighbour_triangle) && good_edges(my_mesh, active_edges, checked_edges, working_edge, point)) {
      closest_point =
          pair(point, line_point_dist(working_edge, point, neighbour_triangle));
    }
  }
  if (closest_point.has_value()) {
    Triangle new_T(working_edge.A(), working_edge.B(),
                   closest_point.value().first);
    my_mesh.add_triangle(working_edge, closest_point.value().first);
    my_mesh.obj_format();
    cout << "New triangle1!" << endl;
    Edge new_edge1(working_edge.A(), closest_point.value().first);
    Edge new_edge2(working_edge.B(), closest_point.value().first);
    number_of_new_edges =
        update_border(new_edge1, new_edge2, active_edges, checked_edges);
    return number_of_new_edges;
  }

  auto [prev, next] =
      find_prev_next(my_mesh, working_edge, active_edges, checked_edges);
  Triangle T = Triangle(working_edge.A(), working_edge.B(), prev);
  Triangle NeighbourT = my_mesh.find_triangle_with_edge(working_edge);
  if (T.is_triangle() && good_orientation(working_edge, prev, NeighbourT) && good_edges(my_mesh, active_edges, checked_edges, working_edge, prev)) {
    assertm(Vector(working_edge.A(), prev).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev point!");
    assertm(Vector(working_edge.B(), prev).get_length() <
                3 * working_edge.get_length(),
            "Weird distance of prev point!");

    my_mesh.add_triangle(working_edge, prev);
    my_mesh.obj_format();
    cout << "New triangle2!" << endl;
    Edge new_edge1(working_edge.A(), prev);
    Edge new_edge2(working_edge.B(), prev);
    number_of_new_edges =
        update_border(new_edge1, new_edge2, active_edges, checked_edges);
    return number_of_new_edges;
  } else if (!is_border(Edge(working_edge.B(), prev), active_edges,
                        checked_edges)) {
    T = Triangle(working_edge.A(), working_edge.B(), next);
    if (T.is_triangle() && good_orientation(working_edge, next, NeighbourT) && good_edges(my_mesh, active_edges, checked_edges, working_edge, next)) {

      assertm(Vector(working_edge.A(), next).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of next point!");
      assertm(Vector(working_edge.B(), next).get_length() <
                  3 * working_edge.get_length(),
              "Weird distance of next point!");

      my_mesh.add_triangle(working_edge, next);
      my_mesh.obj_format();
      cout << "New triangle3!" << endl;
      Edge new_edge1(Edge(working_edge.B(), next));
      Edge new_edge2(Edge(working_edge.A(), next));
      number_of_new_edges =
          update_border(new_edge1, new_edge2, active_edges, checked_edges);
      return number_of_new_edges;
    }
  } else {
    push_edge_to_checked(working_edge, checked_edges);
    return 0;
  }
  assertm(false, "Should not be here!");
  return 1000;
}

void ending(Mesh &my_mesh, vector<Edge> &active_edges,
            vector<Edge> &checked_edges, const Function &F, numeric e_size, const BoundingBox & bounding_box);

void starting(Mesh &my_mesh, vector<Edge> &active_edges,
              vector<Edge> &checked_edges, const Function &F, numeric e_size, const BoundingBox & bounding_box) {
  int round = 0;
  while (!active_edges.empty()) {

      
    round++;
    std::optional<Edge> working_edge = std::nullopt;

    // std::random_shuffle(active_edges.begin(), active_edges.end());

    working_edge = active_edges.back();

    assertm(working_edge.has_value(), "No working edge!");

    delete_from_active(working_edge.value(), active_edges);

    // active_edges.pop_back();

    // cout << "Current working edge: " << endl << working_edge.value() <<
    // endl;

    if (step(my_mesh, active_edges, checked_edges, working_edge.value(), e_size,
             F, bounding_box)) {
      assertm(!is_border(working_edge.value(), active_edges, checked_edges),
              "Working edge found in border!");
    } else {
      assertm(is_checked(working_edge.value(), checked_edges),
              "Checked edge not checked");
    }
    // my_mesh.cout_triangles();
    if (round % 10 == 0) {
      my_mesh.cout_triangles_number();
      
      cout << "Number of edges in active_edges: " << active_edges.size()
           << endl;
      cout << endl;
   
    }
    /*if (round == 700) {
      add_marks(my_mesh, active_edges, checked_edges);
      my_mesh.obj_format();
      return;
    }*/
    my_mesh.obj_format();
  }
  if (checked_edges.empty()) {
    return;
  }
  assertm(my_mesh.is_in_mesh(checked_edges.back()),
          "Checked edge not im mesh!");
  active_edges.clear();
  active_edges = checked_edges;
  checked_edges.clear();

  //cout << "Fixing holes!" << endl;
  //cout << "Number of active edges: " << active_edges.size() << endl;

  ending(my_mesh, active_edges, checked_edges, F, e_size, bounding_box);

  return;
}
void ending(Mesh &my_mesh, vector<Edge> &active_edges,
            vector<Edge> &checked_edges, const Function &F, numeric e_size, const BoundingBox & bounding_box) {
  int round = 0;
  //my_mesh.obj_format();
  //return; 
  while (!active_edges.empty()) {

    round++;
    std::optional<Edge> working_edge = std::nullopt;

    // std::random_shuffle(active_edges.begin(), active_edges.end());

    working_edge = active_edges.back();
    active_edges.pop_back();

    assertm(working_edge.has_value(), "No working edge!");
    assertm(my_mesh.is_in_mesh(working_edge.value()),
            "Working edge not in mesh!");

    Triangle neighbour_triangle =
        my_mesh.find_triangle_with_edge(working_edge.value());

    int new_edges = fix_holes(my_mesh, F, working_edge.value(), active_edges,
                              checked_edges, e_size);
    assertm(new_edges == 0 || new_edges == 1 || new_edges == 2,
            "Wrong number of new edges!");

    if (new_edges > 0) {
      vector<Edge> new_active_edges;
      push_edge_to_active(active_edges.back(), new_active_edges);
      active_edges.pop_back();
      if (new_edges == 2) {
        push_edge_to_active(active_edges.back(), new_active_edges);
        active_edges.pop_back();
      }
      vector<Edge> new_checked_edges =
          connect_edges(checked_edges, active_edges);

      // new_checked_edges.insert(new_checked_edges.end(), active_edges.begin(),
      // active_edges.end());

      starting(my_mesh, new_active_edges, new_checked_edges, F, e_size, bounding_box);
      return;
      active_edges = connect_edges(new_active_edges, new_checked_edges);
      checked_edges.clear();
    }
    // my_mesh.cout_triangles();
    if (round % 10 == 0) {
      my_mesh.cout_triangles_number();
      cout << "Number of edges in active_edges: " << active_edges.size()
           << endl;
      cout << endl;
    }
    
    my_mesh.obj_format();
  }
  return;
}

// a
// b
// c

Mesh BasicAlgorithm::calculate() {

  //my_mesh.obj_format();
  //return my_mesh;
  starting(my_mesh, active_edges, checked_edges, F, e_size, bounding_box);
  vector<Edge> empty_vector;
  //ending(my_mesh, checked_edges, empty_vector, F, e_size);

  my_mesh.obj_format();
  // my_mesh.output();
  return my_mesh;
}

#endif