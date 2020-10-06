// Gmsh project created on Fri Sep 28 14:27:49 2018
SetFactory("OpenCASCADE");
//+
// L longueur plaque en mètre, W largeur en mètre, R=rayon hole en mètre
L=0.125; W=0.018; R=0.003; rc=0.00186;
//+
X1=0.0;X2=W;X3=W/2;X4=R+2*rc;
//+
Y1=0.0;Y2=L;Y3=L/2;Y4=R+2*rc;Y5=3*R+2*rc;
//+
Point(1) = {X1, Y1, 0, 1.0};
//+
Point(2) = {X1, Y2, 0, 1.0};
//+
Point(3) = {X1, Y3, 0, 1.0};
//+
Point(4) = {X1, Y4, 0, 1.0};
//+
Point(5) = {X1, Y5, 0, 1.0};
//+
Point(6) = {X3, Y1, 0, 1.0};
//+
Point(7) = {X3, Y2, 0, 1.0};
//+
Point(8) = {X3, Y3, 0, 1.0};
//+
Point(9) = {X3, Y5, 0, 1.0};
//+
Point(10) = {X2, Y1, 0, 1.0};
//+
Point(11) = {X2, Y2, 0, 1.0};
//+
Point(12) = {X2, Y3, 0, 1.0};
//+
Point(13) = {X2, Y4, 0, 1.0};
//+
Point(14) = {X2, Y5, 0, 1.0};
//+
Point(15) = {X4, Y1, 0, 1.0};
//+
Point(16) = {X1, R, 0, 1.0};
//+
Point(17) = {R, Y1, 0, 1.0};
//+
Point(18) = {X2, R, 0, 1.0};
//+
Point(19) = {X2-R, Y1, 0, 1.0};
//+
Point(20) = {2*X3-Y4, Y1, 0, 1.0};
//+
//+
Line(1) = {17, 15};
//+
Line(2) = {15, 6};
//+
Line(3) = {6, 20};
//+
Line(4) = {20, 19};
//+
Line(5) = {19, 10};
//+
Line(6) = {10, 18};
//+
Line(7) = {18, 13};
//+
Line(8) = {13, 14};
//+
Line(9) = {14, 12};
//+
Line(10) = {12, 11};
//+
Line(11) = {11, 7};
//+
Line(12) = {7, 2};
//+
Line(13) = {2, 3};
//+
Line(14) = {3, 5};
//+
Line(15) = {5, 4};
//+
Line(16) = {4, 16};
//+
Line(17) = {5, 9};
//+
Line(18) = {9, 14};
//+
Line(19) = {3, 8};
//+
Line(20) = {8, 12};
//+
Line(21) = {9, 8};
//+
Line(22) = {8, 7};
//+
Circle(23) = {17, 1, 16};
//+
Circle(24) = {15, 1, 4};
//+
Line Loop(1) = {1, 24, 16, -23};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {2, 3, 4, 5, 6, 7, 8, -18, -17, 15, -24};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {18, 9, -20, -21};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {17, 21, -19, 14};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {20, 10, 11, -22};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {19, 22, 12, 13};
//+
Plane Surface(6) = {6};
//+
Physical Surface(1) = {1};
//+
Physical Surface(2) = {2};
//+
Physical Surface(3) = {3};
//+
Physical Surface(4) = {4};
//+
Physical Surface(5) = {5};
//+
Physical Surface(6) = {6};
//+

//+
Transfinite Line {1, 23} = 300 Using Progression 1;
//+
Transfinite Line {16} = 150 Using Progression 1;
//+
Transfinite Line {24} = 100 Using Progression 1;
//+
Transfinite Line {2, 3} = 20 Using Progression 1;
//+
Transfinite Line {4, 5, 6, 7} = 40 Using Progression 1;
//+
Transfinite Line {15, 8, 17, 18} = 50 Using Progression 1;
//+
Transfinite Line {14, 21, 9, 13, 22, 10} = 100 Using Progression 1;
//+
Transfinite Line {19, 20, 12, 11} = 50 Using Progression 1;
