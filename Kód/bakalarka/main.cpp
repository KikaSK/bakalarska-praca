#include <iostream>
#include <cmath>
#include <ostream>
#include <vector>
#include <ginac/ginac.h>
#include <cassert>
#include <queue>

#include "point.h"
#include "vector.h"
#include "edge.h"
#include "triangle.h"
#include "function.h"
#include "assertm.h"

using namespace std;
using namespace GiNaC;


numeric Newton_Rhapson(realsymbol my_x, ex f, ex df, Point P, int my_case)
{
    numeric precision = 1e-6;

    numeric k;
    if(my_case == 0) k = P.x();
    else if(my_case == 1) k = P.y();
    else k = P.z(); 

    int iterations = 0;
    while (abs(f.subs(my_x == k).evalf())>precision && iterations < 10000000)
    {
        iterations++;
        k -= ex_to<numeric>(f.subs(my_x == k).evalf()/df.subs(my_x == k).evalf());
    }
    cout<<"Number of iterations: " << iterations <<endl;

    return k;
}

Point project(Point P, Vector normal, Function F)
{
    realsymbol my_x("my_x");

    ex parametric_x, parametric_y, parametric_z;
    int my_case = -1;
    if (normal.x() != 0)
    {
        my_case = 0;
        parametric_x = my_x;
        parametric_y = P.y() - normal.y()*P.x()/normal.x() + (normal.y()/normal.x())*my_x;
        parametric_z = P.z() - normal.z()*P.x()/normal.x() + (normal.z()/normal.x())*my_x;
    }
    else if(normal.y() != 0)
    {
        my_case = 1;
        parametric_x = P.x();
        parametric_y = my_x;
        parametric_z = P.z() - normal.z()*P.y()/normal.y() + (normal.z()/normal.y())*my_x;    
    }
    else if(normal.z() != 0)
    {
        my_case = 2;
        parametric_x = P.x();
        parametric_y = P.y();
        parametric_z = my_x;   
    }
    else
    {
        assertm(false, "Normal is a zero vector!");
    }


    ex f = F.get_function().subs( lst{F.get_x() == parametric_x, F.get_y() == parametric_y, F.get_z() == parametric_z} );
    ex df = f.diff(my_x, 1);


    numeric k = Newton_Rhapson(my_x, f, df, P, my_case);
    numeric xx = ex_to<numeric>(parametric_x.subs(my_x == k).evalf());
    numeric yy = ex_to<numeric>(parametric_y.subs(my_x == k).evalf());
    numeric zz = ex_to<numeric>(parametric_z.subs(my_x == k).evalf());

    return Point(xx, yy, zz);

}

Edge get_seed_edge(Point seed, Function F, numeric e_size)
{
    // point to project
    Point P(seed, e_size*(F.get_tangent_at_point(seed).unit()));
    cout<< "Point to proejct: " << P <<endl;
    // direction of projection
    Vector normal = F.get_gradient_at_point(seed).unit();
    cout<<"normal:" << normal<<endl;
    // checking sign of direction
    Point to_check_normal(P, Vector((e_size/2)*normal));
    if (F.is_inside(to_check_normal)) normal = numeric(-1)*normal;
    
    // getting function f(x) from F(x, y, z) and condition that it liest on the line given by P and normal

    Point Q = project(P, normal, F);

    return Edge(seed, Q);
}

Point get_seed_triangle(Edge e, numeric e_size, Function F)
{
    Point center = e.get_midpoint();
    cout<< "My edge: " << e << endl;
    cout<< "Edge midpoint: " << center << endl;
    Vector n_A = F.get_gradient_at_point(e.A()).unit();
    Vector n_B = F.get_gradient_at_point(e.B()).unit();

    Vector normal_center((n_A + n_B)/2);
    Vector edge_vector(e.A(), e.B());
    Vector center_tangent = normal_center^edge_vector;

    assertm(abs(normal_center*center_tangent) < 1e-6, "Not perpendicular!!");
    
    numeric height = sqrt(numeric(3))/2*e_size;
    Point P(center, height*center_tangent.unit());

    Vector dist1(e.A(), P), dist2(e.B(), P);
    cout<< endl << "Distances are: " << dist1.get_length() << " and " << dist2.get_length() << endl;

    cout<< "To project: " << P << endl;

    Point Q = project(P, normal_center, F);
    
    return Q;
}

int main()
{
    realsymbol x("x"), y("y"), z("z");
    Digits = 15;

    ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) -1;
    vector<ex> input_dF;
    input_dF.push_back(2*x);
    input_dF.push_back(2*y);
    input_dF.push_back(2*z);

    queue<Edge> active_edges;

    Function F(x, y, z, input_F, input_dF);
    
    numeric e_size = 0.05;
    Point seed(sqrt(numeric(119)/144), numeric(1)/4, numeric(1)/3);

    //cout<< "Seed point: " << seed <<endl;
    

    //cout<< "Gradient at seed point: " << F.get_gradient_at_point(seed)<<endl;
    //cout<< "Is seed outside: " << F.is_outside(seed) << endl;
    //cout<< F.get_tangent_at_point(seed)<<endl;
    Edge seed_edge = get_seed_edge(seed, F, e_size);
    Point projected = seed_edge.B();
    //cout<< "Projected point: " << projected << endl;
    ex my_func = F.get_function();
    //functional value of projected
    //cout<< "Value at projected point: " << ex_to<numeric>(my_func.subs(lst{x == projected.x(), y == projected.y(), z == projected.z()}).evalf()) <<endl;
    //cout<<endl;
    Vector dist(seed, projected); 
    //cout<< "Distance from seed: " << dist.get_length();
    cout<<endl;
    active_edges.push(seed_edge);

    Point Q = get_seed_triangle(seed_edge, e_size, F);
    cout<< "New point: " << Q << endl;
    Vector dist1(seed, Q), dist2(projected, Q);
    cout<< "Distances from original points: " << dist1.get_length() << " and " << dist2.get_length() << endl;
    cout<< "Value at new point: " << ex_to<numeric>(my_func.subs(lst{x == Q.x(), y == Q.y(), z == Q.z()}).evalf()) << endl;
    Edge e1(seed, Q), e2(projected, Q);
    active_edges.push(e1);
    active_edges.push(e2);
}
