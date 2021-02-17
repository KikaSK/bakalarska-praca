#include <iostream>
#include <cmath>
#include <ostream>
#include <vector>
#include <ginac/ginac.h>
#include <cassert>

#include "point.h"
#include "vector.h"
#include "edge.h"
#include "triangle.h"

using namespace std;
using namespace GiNaC;


int main()
{
    realsymbol x("x"), y("y"), z("z");
    ex F = pow(x, 2) + pow(y, 2) + pow(z, 2) -1;
    
    //Point P(2, 2, 7);
    Point P(1, 7, 6);
    Point Q(1, 7, 6);
    Point R(3, 4, 5);

    Triangle T(P, Q, R);
    cout<< T.get_circumcenter()<<endl;

    cout<<T<<endl;

    /*Vector A(P,R);
    Vector B(1, 2, 3);

    Edge e(P, R);

    cout << e <<endl;
    cout << Vector(Point(0, 0, 0), e.get_midpoint()) <<endl;
    cout<< B.unit() <<endl;

    cout << A << " " << B << endl;

    cout << P.x() << " " << P.y() << " " << P.z() << endl;
    
    cout << F << endl;*/
}
