#include "gtest/gtest.h"
#include "point.h"
#include "vector.h"
#include "edge.h"
#include "triangle.h"
#include "function.h"
#include "mesh.h"

//#include "basic_algorithm.h"

namespace
{
 //TEST examples

 //TEST(ADD, Trivial) { EXPECT_EQ(5 + 3, 20); }
 TEST(DIVISION, Zero) {EXPECT_EQ(5/1, 5); }

 //TEST Point

 TEST(POINT, Equality_True) {
    Point p1(1, 2, 3); Point q1(1, 2, 3); EXPECT_TRUE(p1 == q1);
    Point p2(0, 0, 0); Point q2(0, 0, 0); EXPECT_TRUE(p2 == q2);
    Point p3(0, 0, 0); Point q3(10e-9, 10e-9, 10e-9); EXPECT_TRUE(p3 == q3);
    Point p4(-1.984769347, 2.98764978, -6.56316); Point q4(-1.984769349, 2.98764978, -6.56316); EXPECT_TRUE(p4 == q4);
    Point p5(5, -4, -385683746); Point q5(5.000000002, -4.000000001, -385683746.000000006); EXPECT_TRUE(p5 == q5);
    Point p6(5, -4, -385683746); Point q6(5.00000000, -4.00000000, -385683746.000000006); EXPECT_TRUE(p6 == q6);
    }
 
 
 TEST(POINT, Equality_False) {
    Point p1(0, 0, 0); Point q1(10e-8, 10e-8, 10e-8); EXPECT_FALSE(p1 == q1);
    Point p2(0, 0, 0); Point q2(9*10e-8, 10e-8, 0); EXPECT_FALSE(p2 == q2);
    Point p3(0, 0, 0); Point q3(numeric(10e-8), 0, 0); EXPECT_FALSE(p3 == q3);
    Point p4(0, 0, 0); Point q4(numeric(0.0000001), 0, 0); EXPECT_FALSE(p4 == q4);
    Point p5(1, 2, 3); Point q5(-1, 2, 3); EXPECT_FALSE(p5 == q5);
    Point p6(numeric(-3485.5683206), numeric(8765.985630496), numeric(-9615001.8760505)); Point q6(numeric(-3485.56832060), numeric(8765.9856304960), numeric(-9615001.8760506)); EXPECT_FALSE(p6 == q6);
    Point p7(1, 2, 3); Point q7(-1, 2, 3); EXPECT_FALSE(p7 == q7);
    }

 TEST(POINT, Inequality_False) {
    Point p1(1, 2, 3); Point q1(1, 2, 3); EXPECT_FALSE(p1 != q1);
    Point p2(0, 0, 0); Point q2(0, 0, 0); EXPECT_FALSE(p2 != q2);
    Point p3(0, 0, 0); Point q3(10e-9, 10e-9, 10e-9); EXPECT_FALSE(p3 != q3);
    Point p4(-1.984769347, 2.98764978, -6.56316); Point q4(-1.984769349, 2.98764978, -6.56316); EXPECT_FALSE(p4 != q4);
    Point p5(5, -4, -385683746); Point q5(5.000000002, -4.000000001, -385683746.000000006); EXPECT_FALSE(p5 != q5);
    Point p6(5, -4, -385683746); Point q6(5.00000000, -4.00000000, -385683746.000000006); EXPECT_FALSE(p6 != q6);
    }
 
 
 TEST(POINT, Inequality_True) {
    Point p1(1, 2, 3); Point q1(1.0001, 2, 3); EXPECT_TRUE(p1 != q1);
    Point p2(0, 0, 0); Point q2(10e-8, 10e-8, 10e-8); EXPECT_TRUE(p2 != q2);
    Point p3(0, 0, 0); Point q3(9*10e-8, 10e-8, 0); EXPECT_TRUE(p3!= q3);
    Point p4(0, 0, 0); Point q4(10e-8, 0, 0); EXPECT_TRUE(p4!= q4);
    Point p5(0, 0, 0); Point q5(0.0000001, 0, 0); EXPECT_TRUE(p5!= q5);
    Point p6(1, 2, 3); Point q6(-1, 2, 3); EXPECT_TRUE(p6!= q6);
    Point p7(-3485.5683206, 8765.985630496, -9615001.8760505); Point q7(-3485.56832060, 8765.9856304960, -9615001.8760506); EXPECT_TRUE(p7!= q7);
    Point p8(1, 2, 3); Point q8(-1, 2, 3); EXPECT_TRUE(p8!= q8);
 }

TEST(POINT, Get) {
    Point p1(1, 2, 3); EXPECT_TRUE(p1.x() == 1); EXPECT_TRUE(p1.y() == 2); EXPECT_TRUE(p1.z() == 3); 
    Point p2(1.098954325076408372837608, 2, 3); EXPECT_NEAR(p2.x().to_double(), 1.09895432507640837846406, 1e-15);
    Point p4(1.098954325076408372837608, 2, 3); EXPECT_TRUE(p4.x() == 1.09895432507640837846406);
    Point p3(-1.38536296, 2.28375802, 3.3285740386264); 
        EXPECT_TRUE(p3.x() == -1.38536296); 
        EXPECT_TRUE(p3.y() == 2.28375802); 
        EXPECT_TRUE(p3.z() == 3.3285740386264);
        EXPECT_FALSE(p3.x() == -1.385362961); 
        EXPECT_FALSE(p3.y() == 2.28375802376); 
        EXPECT_FALSE(p3.z() == 3.328574038626);
}

//TEST Vector

TEST(VECTOR, Length)
{
    Vector v1(1, 2, 3); EXPECT_TRUE(v1.get_length() == sqrt(numeric(14)));
    Vector v2(Point(1, 2, 3), Point(0, 0, 0)); EXPECT_TRUE(v2.get_length() == sqrt(numeric(14)));

}

TEST(VECTOR, Unit)
{
    Vector v1(1, 2, 3); EXPECT_TRUE(v1.unit().get_length() - 1 < 10e-8);
}

TEST(VECTOR, Perpendicular)
{
    Vector v1(1, 2, 3); EXPECT_TRUE(v1.get_any_perpendicular()*v1 < 10e-8);
    Vector v2(10e-6, 10e-6, 10e-6); EXPECT_TRUE(v2.get_any_perpendicular()*v2 < 10e-8);
    Vector v3(-1.984769347, 2.98764978, -6.56316); EXPECT_TRUE(v3.get_any_perpendicular()*v3 < 10e-8);
    Vector v4(5, -4, -385683746); EXPECT_TRUE(v4.get_any_perpendicular()*v4 < 10e-8);
    Vector v5(5, -4, -385683746); EXPECT_TRUE(v5.get_any_perpendicular()*v5 < 10e-8);
}

TEST(VECTOR, Cross_Product)
{
    Vector v1(numeric(1), numeric(2), numeric(3)); Vector q1(numeric(8), numeric(9), numeric(3)); 
        EXPECT_TRUE((v1^q1) == Vector(numeric(-21), numeric(21), numeric(-7)));
    Vector v2(numeric(10e-3), numeric(10e-3), numeric(10e-3)); 
        EXPECT_TRUE((v2^q1) == Vector(-numeric(3)/50, numeric(1)/20, numeric(1)/100));
    Vector v4(numeric(5), numeric(-4), numeric(-385683746)); 
        EXPECT_TRUE((v4^q1) == Vector(numeric(3471153702), numeric(-3085469983), numeric(77)));
}

//TEST Edge

TEST(EDGE, Length_Midpoint)
{
    Edge e1(Point(1, 2, 3), Point(0, 0, 0)); 
        EXPECT_TRUE(e1.get_length() == sqrt(numeric(14)));
        EXPECT_TRUE(e1.get_midpoint() == Point(numeric(1)/2, numeric(1), numeric(3)/2));
}

TEST(EDGE, Eq_Ineq)
{
    Edge e1(Point(1, 2, 3), Point(0, 0, 0));
        EXPECT_TRUE(e1.get_length() == sqrt(numeric(14)));
        EXPECT_TRUE(e1.get_midpoint() == Point(numeric(1)/2, numeric(1), numeric(3)/2));
         
    Edge e2(Point(1, 2, 3), Point(0, 0, 0));
        EXPECT_TRUE(e1 == e2 && !(e1 != e2));
    Edge e3(Point(1, 2, 3), Point(0, 10e-7, 0));
        EXPECT_TRUE(!(e1 == e3) && e1 != e3);
    Edge e4(Point(3, 4, 0), Point(6, 7, 3));
        EXPECT_TRUE(e1 != e4 && !(e1 == e4));
}

//TEST Triangle

//Dorobit lepsie testy na is_in_triangle
TEST(TRIANGLE, Functions){
    Point P1(0, 0, 0);

    Point P21(1, 0, 0);
    Point P31(0, 1, 0);
    Triangle T1(P1, P21, P31);
        EXPECT_TRUE(T1.is_triangle());
        EXPECT_TRUE(T1.get_gravity_center() == Point(numeric(1)/3, numeric(1)/3, 0));
        EXPECT_TRUE(T1.get_circumcenter() == Point(numeric(1)/2, numeric(1)/2, 0));
    Point P41(0.25, 0.25, 0);
        EXPECT_TRUE(T1.is_in_triangle(P41));
    Point P51(0.25, 0.25, 1e-8);
        EXPECT_TRUE(T1.is_in_triangle(P51));
    Point P61(0.25, 0.25, 10e-8);
        EXPECT_FALSE(T1.is_in_triangle(P61));
    
    EXPECT_TRUE(T1.get_normal() == Vector(0, 0, 1) || T1.get_normal() == Vector(0, 0, -1));

    Point P22(2, 3, 4);
    Point P32(4, 2, 3);
    Triangle T2(P1, P22, P32);
        EXPECT_TRUE(T2.is_triangle());
        EXPECT_TRUE(T2.get_circumcenter() == Point(numeric(87)/55, numeric(29)/22, numeric(203)/110));
        EXPECT_TRUE(T2.get_gravity_center() == Point(numeric(2), numeric(5)/3, numeric(7)/3));
    Point P42 = T2.get_gravity_center();
        EXPECT_TRUE(T2.is_in_triangle(P42));
    Point P52 = T2.get_circumcenter();
        EXPECT_TRUE(T2.is_in_triangle(P52));
    Point P62 = Point(P52, Vector(10e-8, 10e-9, 10e-10));
        EXPECT_TRUE(T2.is_in_triangle(P62));
    Point P72 = Point(P52, Vector(10e-6, 10e-9, 10e-10));
        EXPECT_FALSE(T2.is_in_triangle(P72));
    
    EXPECT_TRUE(T2.get_normal().get_length() - 1 < 10e-8);
    EXPECT_TRUE(T2.get_normal()*Vector(P1, P22) < 10e-8);
    EXPECT_TRUE(T2.get_normal()*Vector(P1, P32) < 10e-8);
}

//TEST Function
TEST(FUNCTION, Evrything){
    realsymbol x("x"), y("y"), z("z");

    ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
    std::vector<ex> input_dF;

    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));

    Function F(x, y, z, input_F, input_dF);

    Point P11(numeric(0), numeric(0), numeric(1));
    Point P12(numeric(0), numeric(1), numeric(0));
    Point P13(numeric(-1), numeric(0), numeric(0));
    Point P14(numeric(1)/2, numeric(1)/2, numeric(1)/sqrt(numeric(2)));
    Point P15(-numeric(1)/2, 0, sqrt(numeric(3))/2);
    Point P16(numeric(0.644), numeric(0.1), -sqrt(numeric(0.575264)));
    Point P17(numeric(-0.432), numeric(1)/2, -sqrt(numeric(0.563376)));
    Point P18(numeric(-0.432), sqrt(numeric(0.773376)), -numeric(0.2));

    // get

    EXPECT_TRUE(F.get_x() == x);
    EXPECT_TRUE(F.grad_x() == 2*x);
    // ggradient (2x, 2y, 2z) 

    // gradient
    EXPECT_TRUE(F.get_gradient_at_point(P11) == 2*Vector(Point(0, 0, 0), P11));
    EXPECT_TRUE(F.get_gradient_at_point(P12) == 2*Vector(Point(0, 0, 0), P12));
    EXPECT_TRUE(F.get_gradient_at_point(P13) == 2*Vector(Point(0, 0, 0), P13));
    EXPECT_TRUE(F.get_gradient_at_point(P14) == 2*Vector(Point(0, 0, 0), P14));
    EXPECT_TRUE(F.get_gradient_at_point(P15) == 2*Vector(Point(0, 0, 0), P15));
    EXPECT_TRUE(F.get_gradient_at_point(P16) == 2*Vector(Point(0, 0, 0), P16));
    EXPECT_TRUE(F.get_gradient_at_point(P17) == 2*Vector(Point(0, 0, 0), P17));
    EXPECT_TRUE(F.get_gradient_at_point(P18) == 2*Vector(Point(0, 0, 0), P18));

    // gradient and tangent are perpendicular
    EXPECT_TRUE(F.get_tangent_at_point(P11)*F.get_gradient_at_point(P11) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P12)*F.get_gradient_at_point(P12) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P13)*F.get_gradient_at_point(P13) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P14)*F.get_gradient_at_point(P14) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P15)*F.get_gradient_at_point(P15) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P16)*F.get_gradient_at_point(P16) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P17)*F.get_gradient_at_point(P17) < 10e-8);
    EXPECT_TRUE(F.get_tangent_at_point(P18)*F.get_gradient_at_point(P18) < 10e-8);

    // is_inside, is_outside, is_on, eval_at_point
    EXPECT_TRUE(F.is_on(P11));
    EXPECT_TRUE(F.is_on(P12));
    EXPECT_TRUE(F.is_on(P13));
    EXPECT_TRUE(F.is_on(P14));
    EXPECT_TRUE(F.is_on(P15));
    EXPECT_TRUE(F.is_on(P16));
    EXPECT_TRUE(F.is_on(P17));
    EXPECT_TRUE(F.is_on(P18));

    EXPECT_TRUE(F.eval_at_point(P11) - (P11.x()*P11.x() + P11.y()*P11.y() + P11.z()*P11.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P12) - (P12.x()*P12.x() + P12.y()*P12.y() + P12.z()*P12.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P13) - (P13.x()*P13.x() + P13.y()*P13.y() + P13.z()*P13.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P14) - (P14.x()*P14.x() + P14.y()*P14.y() + P14.z()*P14.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P15) - (P15.x()*P15.x() + P15.y()*P15.y() + P15.z()*P15.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P16) - (P16.x()*P16.x() + P16.y()*P16.y() + P16.z()*P16.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P17) - (P17.x()*P17.x() + P17.y()*P17.y() + P17.z()*P17.z() -numeric(1) ) < 10e-8);
    EXPECT_TRUE(F.eval_at_point(P18) - (P18.x()*P18.x() + P18.y()*P18.y() + P18.z()*P18.z() -numeric(1) ) < 10e-8);

    Point P21(numeric(1), numeric(1), numeric(1));
    Point P22(numeric(-1), numeric(3), numeric(0.2));
    Point P23(numeric(0.564), numeric(1.05), numeric(0));
    Point P24(numeric(0.5), numeric(0.5), numeric(0.9));

    EXPECT_TRUE(F.is_outside(P21));
    EXPECT_TRUE(F.is_outside(P22));
    EXPECT_TRUE(F.is_outside(P23));
    EXPECT_TRUE(F.is_outside(P24));
    
    EXPECT_TRUE(F.eval_at_point(P21) - (P21.x()*P21.x() + P21.y()*P21.y() + P21.z()*P21.z() -numeric(1) ) <10e-8);
    EXPECT_TRUE(F.eval_at_point(P22) - (P22.x()*P22.x() + P22.y()*P22.y() + P22.z()*P22.z() -numeric(1) ) <10e-8);
    EXPECT_TRUE(F.eval_at_point(P23) - (P23.x()*P23.x() + P23.y()*P23.y() + P23.z()*P23.z() -numeric(1) ) <10e-8);
    EXPECT_TRUE(F.eval_at_point(P24) - (P24.x()*P24.x() + P24.y()*P24.y() + P24.z()*P24.z() -numeric(1) ) <10e-8);

    Point P31(numeric(0.1), numeric(-0.3), numeric(-0.23212312413234));
    Point P32(numeric(0.5), numeric(0.5), numeric(0.5));
    Point P33(numeric(-0.3), numeric(-0.321), numeric(-0.51));
    Point P34(numeric(1)/4, numeric(-2)/5, numeric(-1)/7);
    
    EXPECT_TRUE(F.is_inside(P31));
    EXPECT_TRUE(F.is_inside(P32));
    EXPECT_TRUE(F.is_inside(P33));
    EXPECT_TRUE(F.is_inside(P34));

    EXPECT_TRUE(F.eval_at_point(P31) - (P31.x()*P31.x() + P31.y()*P31.y() + P31.z()*P31.z() -numeric(1) ) <10e-8);
    EXPECT_TRUE(F.eval_at_point(P32) - (P32.x()*P32.x() + P32.y()*P32.y() + P32.z()*P32.z() -numeric(1) ) <10e-8);
    EXPECT_TRUE(F.eval_at_point(P33) - (P33.x()*P33.x() + P33.y()*P33.y() + P33.z()*P33.z() -numeric(1) ) <10e-8);
    EXPECT_TRUE(F.eval_at_point(P34) - (P34.x()*P34.x() + P34.y()*P34.y() + P34.z()*P34.z() -numeric(1) ) <10e-8);


    //TEST outside normal, imported some triangles from triangulation

    Point T11(1, 0, 0), T12(0.9797959183673469365, -0.2000000000000000111, 0), T13(0.9798984706880384728, -0.098989846549567533174, -0.17320508075688773897);
    Point T21(0.9797959183673469365, -0.2000000000000000111, 0), T22(0.9798984706880384728, -0.098989846549567533174, -0.17320508075688773897), T23(0.93995085126291774415, -0.29507759059056836527, -0.17152741299174244788);
    Point T31(0.9798984706880384728, -0.098989846549567533174, -0.17320508075688773897), T32(1, 0, 0), T33(0.9799755048755642622, 0.10112564628025310584, -0.17152741522629678692);

    Triangle T1(T11, T12, T13);
    Vector normal1 = F.outside_normal(T1);
    numeric delta = 0.05;
    EXPECT_TRUE(F.is_outside(Point(T11, delta*normal1)));
    EXPECT_TRUE(F.is_outside(Point(T12, delta*normal1)));
    EXPECT_TRUE(F.is_outside(Point(T13, delta*normal1)));
    EXPECT_TRUE(F.is_inside(Point(T11, -delta*normal1)));
    EXPECT_TRUE(F.is_inside(Point(T12, -delta*normal1)));
    EXPECT_TRUE(F.is_inside(Point(T13, -delta*normal1)));

    Triangle T2(T21, T22, T23);
    Vector normal2 = F.outside_normal(T2);
    EXPECT_TRUE(F.is_outside(Point(T21, delta*normal2)));
    EXPECT_TRUE(F.is_outside(Point(T22, delta*normal2)));
    EXPECT_TRUE(F.is_outside(Point(T23, delta*normal2)));
    EXPECT_TRUE(F.is_inside(Point(T21, -delta*normal2)));
    EXPECT_TRUE(F.is_inside(Point(T22, -delta*normal2)));
    EXPECT_TRUE(F.is_inside(Point(T23, -delta*normal2)));

    Triangle T3(T31, T32, T33);
    Vector normal3 = F.outside_normal(T3);
    EXPECT_TRUE(F.is_outside(Point(T31, delta*normal3)));
    EXPECT_TRUE(F.is_outside(Point(T32, delta*normal3)));
    EXPECT_TRUE(F.is_outside(Point(T33, delta*normal3)));
    EXPECT_TRUE(F.is_inside(Point(T31, -delta*normal3)));
    EXPECT_TRUE(F.is_inside(Point(T32, -delta*normal3)));
    EXPECT_TRUE(F.is_inside(Point(T33, -delta*normal3)));
    
    EXPECT_TRUE(normal1.get_length() - 1 < 10e-8);
    EXPECT_TRUE(normal2.get_length() - 1 < 10e-8);
    EXPECT_TRUE(normal3.get_length() - 1 < 10e-8);

    EXPECT_TRUE(normal1*Vector(T11, T12) < 10e-8);
    EXPECT_TRUE(normal1*Vector(T13, T12) < 10e-8);
    EXPECT_TRUE(normal1*Vector(T11, T13) < 10e-8);

    EXPECT_TRUE(normal2*Vector(T21, T22) < 10e-8);
    EXPECT_TRUE(normal2*Vector(T23, T22) < 10e-8);
    EXPECT_TRUE(normal2*Vector(T21, T23) < 10e-8);

    EXPECT_TRUE(normal3*Vector(T31, T32) < 10e-8);
    EXPECT_TRUE(normal3*Vector(T33, T32) < 10e-8);
    EXPECT_TRUE(normal3*Vector(T31, T33) < 10e-8);

/* 
another triangles to test:

0.9797959183673469365 -0.2000000000000000111 0
0.93995085126291774415 -0.29507759059056836527 -0.17152741299174244788
0.9203935026188246701 -0.3909825131671295809 0.0029181045385228633043

1 0 0
0.9797959183673469365 -0.2000000000000000111 0
0.97994933756728075245 -0.098994984108335860966 0.17291422376552230617

0.9797959183673469365 -0.2000000000000000111 0
0.9203935026188246701 -0.3909825131671295809 0.0029181045385228633043
0.9410062808848314827 -0.29002292853999145017 0.17433852328699591384

0.9797959183673469365 -0.2000000000000000111 0
0.97994933756728075245 -0.098994984108335860966 0.17291422376552230617
0.9410062808848314827 -0.29002292853999145017 0.17433852328699591384
*/
}

//TEST Mesh
TEST(MESH, Evrything){
    Point T11(1, 0, 0), T12(0.9797959183673469365, -0.2000000000000000111, 0), T13(0.9798984706880384728, -0.098989846549567533174, -0.17320508075688773897);
    Point T21(0.9797959183673469365, -0.2000000000000000111, 0), T22(0.9798984706880384728, -0.098989846549567533174, -0.17320508075688773897), T23(0.93995085126291774415, -0.29507759059056836527, -0.17152741299174244788);
    Point T31(0.9798984706880384728, -0.098989846549567533174, -0.17320508075688773897), T32(1, 0, 0), T33(0.9799755048755642622, 0.10112564628025310584, -0.17152741522629678692);
    Point T41(0.9797959183673469365, -0.2000000000000000111, 0), T42(0.93995085126291774415, -0.29507759059056836527, -0.17152741299174244788), T43(0.9203935026188246701, -0.3909825131671295809, 0.0029181045385228633043);
    Point T51(1, 0, 0), T52(0.9797959183673469365, -0.2000000000000000111, 0), T53(0.97994933756728075245, -0.098994984108335860966, 0.17291422376552230617);
    Point T61(0.9797959183673469365, -0.2000000000000000111, 0), T62(0.9203935026188246701, -0.3909825131671295809, 0.0029181045385228633043), T63(0.9410062808848314827, -0.29002292853999145017, 0.17433852328699591384);
    Point T71(0.9797959183673469365, -0.2000000000000000111, 0), T72(0.97994933756728075245, -0.098994984108335860966, 0.17291422376552230617), T73(0.9410062808848314827, -0.29002292853999145017, 0.17433852328699591384);

    Triangle T1(T11, T12, T13);
    Triangle T2(T21, T22, T32);
    Triangle T3(T31, T32, T33);
    Triangle T4(T41, T42, T43);
    Triangle T5(T51, T52, T52);
    Triangle T6(T61, T62, T63);
    Triangle T7(T71, T72, T73);

    Mesh my_mesh(T1);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T11, T12)) == T1);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T12, T13)) == T1);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T11, T13)) == T1);

    my_mesh.add_triangle(Edge(T21, T22), T23);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T11, T12)) == T1);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T11, T13)) == T1);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T22, T23)) == T2);
    EXPECT_TRUE(my_mesh.find_triangle_with_edge(Edge(T21, T23)) == T2);
    
    my_mesh.add_triangle(Edge(T31, T32), T33);
    my_mesh.add_triangle(Edge(T41, T42), T43);
    my_mesh.add_triangle(Edge(T51, T52), T53);
    my_mesh.add_triangle(Edge(T61, T62), T63);
    my_mesh.add_triangle(Edge(T71, T72), T73);

}

}

// namespace
