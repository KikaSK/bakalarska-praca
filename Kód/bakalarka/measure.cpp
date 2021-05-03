
#include <cassert>
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "assertm.h"
#include "point.h"
#include "vector.h"
#include "edge.h"
#include "triangle.h"
#include "function.h"
#include "algorithms.h"
#include "bounding_box.h"

using namespace std;
using namespace GiNaC;

int is_in_points(const Point P, const vector<pair<Point, vector<int> > >points){
    for(int i = 0; i < points.size(); ++i){
        if(points[i].first == P)
            return i;
    }
    return -1;
}

int is_in_edges(const Edge & E, vector<pair<Edge, vector<Triangle> > > edges){
    for (int i = 0; i < edges.size(); ++i){
        if(edges[i].first == E){
            return i;
        }
    }
    return -1;
}
bool is_in_pointers(const int p, const vector<int>pointers){
    for(auto i : pointers){
        if(i==p){
            return true;
        }
    }
    return false;
}


vector<numeric> average_gc_distance(const vector<Triangle> &mesh_triangles, 
                                    const Function &F, const numeric e_size) {
  int n = mesh_triangles.size();
  vector<numeric>result;
  numeric sum = 0;
  numeric min = 1000*e_size;
  numeric max = -1;
  for(auto T: mesh_triangles){
    Point gc = T.get_gravity_center();
    Vector direction = F.get_gradient_at_point(gc).unit();
    Point gc_proj = project(gc, direction, F, e_size);
    numeric dist = Vector(gc, gc_proj).get_length();
    sum += dist;
    if(dist<min) min = dist;
    if(dist>max) max = dist;
  }
  result.push_back(sum/n);
  result.push_back(max);
  result.push_back(min);
  return result;
}

numeric average_side_length(const vector<Triangle> &mesh_triangles, const vector<Edge> &bounding_edges) {
  int n = mesh_triangles.size()*3 + bounding_edges.size();
  numeric sum = 0;
  for(auto T: mesh_triangles){
    sum += T.AB().get_length();
    sum += T.BC().get_length();
    sum += T.CA().get_length();
  }
  for(auto e : bounding_edges){
    sum += e.get_length();
  }
  return sum/n;
}
numeric average_side_ratio(const vector<Triangle> &mesh_triangles) {
  int n = mesh_triangles.size();
  numeric sum = 0;
  for(auto T: mesh_triangles){
    numeric a = T.BC().get_length();
    numeric b = T.CA().get_length();
    numeric c = T.AB().get_length();

    numeric biggest = std::max(a, std::max(b, c));
    numeric smallest = std::min(a, std::min(b, c));
    numeric ratio = biggest/smallest;
    sum = sum+ratio;
  }
  return sum/n;
}

pair<numeric, numeric> min_max_normal_angle(const vector<Triangle> &mesh_triangles, 
const vector<pair<Edge, vector<Triangle>>> &mesh_edges, const Function &F, const numeric e_size){
    numeric min_angle = 4;
    numeric max_angle =-1;
    for(auto E : mesh_edges)
    {
        if(E.second.size() == 1){
            continue;
        }
        Edge e = E.first;
        vector<Triangle> T = E.second;
        Vector n1 = F.outside_normal(T[0], e_size);
        Vector n2 = F.outside_normal(T[1], e_size);
        numeric angle = acos(n1*n2);
        if(angle>max_angle) max_angle = angle;
        if(angle<min_angle) min_angle = angle;
    }
    return pair(min_angle, max_angle);
}

pair<numeric, numeric> mean_std(const vector< pair<Point, vector<int> > > & points){
    numeric sum_avg = 0;
    numeric sum_std = 0;
    for (auto p : points){
        Point P = p.first;
        numeric sum = 0;
        for(auto index : p.second){
            numeric dist = Vector(P, points[index].first).get_length();
            sum += dist;
        }
        numeric avg = sum/p.second.size();
        sum = 0;
        for (auto index : p.second){
            numeric clen = (Vector(P, points[index].first).get_length() - avg);
            clen *= clen; 
            sum += clen;
        }
        numeric std = sqrt(sum/p.second.size());
        sum_avg += avg;
        sum_std += std;
    }
    numeric avg_avg = sum_avg/points.size();
    numeric std_avg = sum_std/points.size();
    return pair(avg_avg, std_avg);
}




void measure (const vector<pair<Point, vector<int> > > &mesh_points, 
            const vector<Triangle> &mesh_triangles, 
            const vector<pair<Edge, vector<Triangle>> > &mesh_edges,
            const vector<Edge> &bounding_edges, const Function &F, 
            const numeric e_size, const string &name){
    
    std::ofstream out("measure_data/" + name + ".out");
    numeric avg_side_length = average_side_length(mesh_triangles, bounding_edges);
    numeric avg_max_side_ratio = average_side_ratio(mesh_triangles);
    //numeric avg_tri_area = average_triangle_area();

    auto avg_max_min_gc_dist = average_gc_distance(mesh_triangles, F, e_size);
    numeric avg_gc_dist = avg_max_min_gc_dist[0];
    numeric max_gc_dist = avg_max_min_gc_dist[1];
    numeric min_gc_dist = avg_max_min_gc_dist[2];
    //numeric avg_gc_dist = average_gc_distance(mesh_triangles, F, e_size);
    numeric ideal_triangle_area = e_size*e_size*sqrt(numeric(3))/4;
    auto [avg_mean_neighbour_points_dist, avg_std_neighbour_points_dist] = mean_std(mesh_points);
    auto [min_normal_angle, max_normal_angle] = min_max_normal_angle(mesh_triangles, mesh_edges, F, e_size);

    out << "DATA:" << endl;
    out << "Edge size: " << e_size << endl;
    //out << "Area of triangle with edge size " << e_size << " is: " << ideal_triangle_area << endl;
    out << "Number of triangles in mesh: " << mesh_triangles.size() << endl;
    out << "Number of points in mesh: " << mesh_points.size() << endl << endl;
    
    out << "AVERAGES: " << endl;
    out << "Average side length: " << avg_side_length << endl;
    out << "Average maximum triangle side ratio: " << avg_max_side_ratio << endl;
    out << "Average of mean of nieghbour points distance: " << avg_mean_neighbour_points_dist << endl;
    out << "Average of standard deviation of nieghbour points distance: " << avg_std_neighbour_points_dist << endl << endl;
//   out << "Average area of triangles: " << avg_tri_area << endl << endl;

    out << "MIN and MAX: " << endl;
    out << "Hausdorff distance - Max gc distance: " << max_gc_dist << endl;
    out << "Min neighbour triangles normals angle: " << min_normal_angle << endl;
    out << "Max neighbour triangles normals angle: " << max_normal_angle << endl << endl;

    out << "RATIOS:" << endl;
//   out << "Avergae trinagle area to ideal triangle area: " << avg_tri_area/ideal_triangle_area << endl;
    out << "Average gravity center distance to edge size: " << avg_gc_dist/e_size << endl; 
    out << "Average edge size to wished edge size: " << avg_side_length/e_size << endl;
    out << "Max gravity center distance to edge size: " << max_gc_dist/e_size << endl;

}


void run_input(const int i, const string folder, const string index){

    string index_str = index;
    
    std::ifstream input_file ("./inputs" + folder +"/input" + to_string(i), std::ifstream::in);
    assertm(input_file.is_open(), "Failed opening the file!");
    vector<string>parsed_input;
    for(string str; getline(input_file, str);){
        string s;
        bool writing = false;
        for(char c : str){
        if(writing && c!='"'){
            s.push_back(c);
        }
        if(c=='"') writing = !writing;
        }
        parsed_input.push_back(s);
    }

    realsymbol x("x"), y("y"), z("z");

    symtab table;
    table["x"] = x;
    table["y"] = y;
    table["z"] = z;
    
    parser reader(table);

    string name = index_str + "_" + parsed_input[0] + "_" + parsed_input[1];

    ex input_F = reader(parsed_input[2]);
    
    numeric min_x = ex_to<numeric>(stod(parsed_input[3]));
    numeric max_x = ex_to<numeric>(stod(parsed_input[4]));
    numeric min_y = ex_to<numeric>(stod(parsed_input[5]));
    numeric max_y = ex_to<numeric>(stod(parsed_input[6]));
    numeric min_z = ex_to<numeric>(stod(parsed_input[7]));
    numeric max_z = ex_to<numeric>(stod(parsed_input[8]));
    BoundingBox my_bounding_box(min_x, max_x, min_y, max_y, min_z, max_z);

    numeric e_size = ex_to<numeric>(stod(parsed_input[9]));

    vector<ex> input_dF;
    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));

    Function F(x, y, z, input_F, input_dF);

    string file_name = "/" + index_str + "_" + name + ".obj";


    std::ifstream obj_file ("./outputs" + folder + file_name, std::ifstream::in);
    assertm(obj_file.is_open(), "Failed opening the file!");
    vector<pair<Point, vector<int> > >mesh_points;
    vector<Triangle>mesh_triangles;
    vector<pair<Edge, vector<Triangle>> >mesh_edges;
    vector<Edge>bounding_edges;
    string type;
    int round = 0;
    vector<Point>last_points;
    for(string type; obj_file >> type;){
        if(type == "v")
        {   
            double a, b, c;
            input_file >> a >> b >> c;
            numeric n_a, n_b, n_c;
            n_a = ex_to<numeric>(a);
            n_b = ex_to<numeric>(b);
            n_c = ex_to<numeric>(c);

            Point P(n_a, n_b, n_c);
            if(round == 2)
            {
                Triangle T(last_points[0], last_points[1], P);
                last_points.clear();
                mesh_triangles.push_back(T);
                int edge_index0 = is_in_edges(T.AB(), mesh_edges);
                int edge_index1 = is_in_edges(T.AB(), mesh_edges);
                int edge_index2 = is_in_edges(T.AB(), mesh_edges);
                if(edge_index0 == -1){
                    vector<Triangle>Tvec;
                    Tvec.push_back(T);
                    mesh_edges.push_back(pair(T.AB(), Tvec));
                }
                else{
                    mesh_edges[edge_index0].second.push_back(T);
                }
                if(edge_index1 == -1){
                    vector<Triangle>Tvec;
                    Tvec.push_back(T);
                    mesh_edges.push_back(pair(T.BC(), Tvec));
                }
                else{
                    mesh_edges[edge_index1].second.push_back(T);
                }
                if(edge_index2 == -1){
                    vector<Triangle>Tvec;
                    Tvec.push_back(T);
                    mesh_edges.push_back(pair(T.CA(), Tvec));
                }
                else{
                    mesh_edges[edge_index2].second.push_back(T);
                }
                if(my_bounding_box.new_bounding_edge(T.AB())){
                    bounding_edges.push_back(T.AB());
                }
                if(my_bounding_box.new_bounding_edge(T.BC())){
                    bounding_edges.push_back(T.BC());
                }
                if(my_bounding_box.new_bounding_edge(T.CA())){
                    bounding_edges.push_back(T.CA());
                }
            }

            int ind = is_in_points(P, mesh_points);
            if(ind == -1){
                pair<Point, vector<int> > par = pair(P, vector<int>());
                mesh_points.push_back(par);
            }
        round = (round+1)%3;
        }
    }
    for(Triangle T : mesh_triangles){
        int indA = is_in_points(T.A(), mesh_points);
        int indB = is_in_points(T.B(), mesh_points);
        int indC = is_in_points(T.C(), mesh_points);
        
        if(!is_in_pointers(indB, mesh_points[indA].second))
            mesh_points[indA].second.push_back(indB);
        if(!is_in_pointers(indC, mesh_points[indA].second))
            mesh_points[indA].second.push_back(indC);
        
        if(!is_in_pointers(indA, mesh_points[indB].second))
            mesh_points[indB].second.push_back(indA);
        if(!is_in_pointers(indC, mesh_points[indB].second))
            mesh_points[indB].second.push_back(indC);

        if(!is_in_pointers(indA, mesh_points[indC].second))
            mesh_points[indC].second.push_back(indA);
        if(!is_in_pointers(indB, mesh_points[indC].second))
            mesh_points[indC].second.push_back(indB);
    }

    for( auto E : mesh_edges){

    }

    measure(mesh_points, mesh_triangles, mesh_edges, bounding_edges, F, e_size, name);
    
}

void run_all(const string folder, const string name){
    int beg = 0;
    int end = 1;

    for(int i = beg; i<=end; ++i){
        run_input(i, folder, name);
    }
}


int main(){
    run_all("/sphere", "final_measures");
}