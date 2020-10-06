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