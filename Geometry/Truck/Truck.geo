// Gmsh project created on Fri Sep  7 09:27:28 2018
//+


L=1;
Lpml = 1.5;

sizemesh = 0.8;

SetFactory("OpenCASCADE");

Box(1) = {0, 0, 0.5, 1.6, 2.3, 2.5};

Cylinder(2) = {1.1, 0, 0.5, 0, 0.3, 0, 0.5, 2*Pi};

Cylinder(20) = {1.1, 2, 0.5, 0, 0.3, 0, 0.5, 2*Pi};

Box(4) = {2, 0, 1, 9, 2.3, 2.8};
//+
Box(5) = {5.75, 0, 0.5, 2.5, 2.3, 0.5};
//+
Cylinder(6) = {6.3, 0, 0.5, 0, 0.3, 0, 0.5, 2*Pi};
//+
Cylinder(7) = {7.7, 0, 0.5, 0, 0.3, 0, 0.5, 2*Pi};
//+
Cylinder(8) = {6.3, 2, 0.5, 0, 0.3, 0, 0.5, 2*Pi};
//+
Cylinder(9) = {7.7, 2, 0.5, 0, 0.3, 0, 0.5, 2*Pi};
//+
Box(10) = {1.6, 1.7, 0.8, 4.15, 0.4, 0.2};
//+
Box(11) = {1.6, 0.2, 0.8, 4.15, 0.4, 0.2};
//+
Box(13) = {-L, -L, 0, 11+2*L, 2+4+L, 3.8+L};
//+
Surface Loop(15) = {1, 3, 5, 4, 2, 6};
//+
Surface Loop(16) = {6, 1, 3, 5, 4, 2};
//+
Surface Loop(17) = {21, 23, 19, 24, 22, 20};

BooleanDifference{ Volume{13}; Delete; }{ Volume{4}; Volume{9}; Volume{7}; Volume{6}; Volume{8}; Volume{11}; Volume{10}; Volume{5}; Volume{2}; Volume{20}; Volume{1}; Delete; }


Characteristic Length {3, 1, 5, 7, 2, 4, 45, 46, 56, 55, 58, 38, 37, 43, 44, 30, 29, 42, 41, 57, 49, 48, 47, 50, 51, 54, 53, 52, 75, 61, 64, 70, 76, 59, 60, 71, 39, 40, 62, 63, 66, 65, 68, 69, 72, 73, 74, 95, 94, 79, 80, 93, 92, 96, 35, 36, 90, 91, 34, 33, 88, 89, 87, 98, 100, 99, 78, 81, 85, 86, 84, 82, 83, 97, 13, 15, 19, 17, 26, 18, 20, 28, 25, 22, 24, 23, 21, 27, 31, 32, 67} = sizemesh;


//+
Point(101) = {-1-Lpml, 6, 0, sizemesh};
//+
Point(102) = {-1-Lpml, 6, 4.8+Lpml, sizemesh};
//+
Point(103) = {-1, 6, 4.8+Lpml, sizemesh};
//+
Point(104) = {-1-Lpml, -1-Lpml, 0, sizemesh};
//+
Point(105) = {-1-Lpml, -1-Lpml, 4.8, sizemesh};
//+
Point(106) = {-1-Lpml, -1-Lpml, 4.8+Lpml, sizemesh};
//+
Point(107) = {-1-Lpml, -1, 0, sizemesh};
//+
Point(108) = {-1-Lpml, -1, 4.8, sizemesh};
//+
Point(109) = {-1-Lpml, -1, 4.8+Lpml, sizemesh};
//+
Point(110) = {-1, -1-Lpml, 0, sizemesh};
//+
Point(111) = {-1, -1-Lpml, 4.8, sizemesh};
//+
Point(112) = {-1, -1-Lpml, 4.8+Lpml, sizemesh};
//+
Point(113) = {-1, -1, 4.8+Lpml, sizemesh};
//+
Point(114) = {12, -1, 4.8+Lpml, sizemesh};
//+
Point(115) = {12, -1-Lpml, 0, sizemesh};
//+
Point(116) = {12, -1-Lpml, 4.8, sizemesh};
//+
Point(117) = {12, -1-Lpml, 4.8+Lpml, sizemesh};
//+
Point(118) = {12, 6, 4.8+Lpml, sizemesh};
//+
Point(119) = {-1-Lpml, 6, 4.8, sizemesh};



//+
Line(173) = {102, 109};
//+
Line(174) = {109, 106};
//+
Line(175) = {106, 112};
//+
Line(176) = {112, 113};
//+
Line(177) = {113, 109};
//+
Line(178) = {109, 108};
//+
Line(179) = {108, 105};
//+
Line(180) = {105, 111};
//+
Line(181) = {111, 112};
//+
Line(182) = {108, 107};
//+
Line(183) = {107, 104};
//+
Line(184) = {104, 105};
//+
Line(185) = {106, 105};
//+
Line(186) = {111, 110};
//+
Line(187) = {104, 110};
//+
Line(188) = {22, 110};
//+
Line(189) = {107, 22};
//+
Line(190) = {24, 101};
//+
Line(192) = {102, 103};
//+
Line(193) = {103, 23};
//+
Line(194) = {102, 119};
//+
Line(195) = {119, 101};
//+
Line(196) = {101, 107};
//+
Line(197) = {119, 108};
//+
Line(198) = {103, 113};
//+
Line(199) = {113, 114};
//+
Line(200) = {114, 117};
//+
Line(201) = {117, 116};
//+
Line(202) = {116, 115};
//+
Line(203) = {115, 110};
//+
Line(204) = {111, 116};
//+
Line(205) = {117, 112};
//+
Line(206) = {114, 118};
//+
Line(207) = {118, 103};
//+
Line(208) = {119, 23};
//+
Line(209) = {21, 113};
//+
Line(210) = {26, 114};
//+
Line(211) = {26, 116};
//+
Line(212) = {25, 115};
//+
Line(213) = {108, 21};
//+
Line(214) = {21, 111};
//+
Line(215) = {118, 27};
//+
Line Loop(91) = {207, 198, 199, 206};
//+
Plane Surface(83) = {91};
//+
Line Loop(92) = {38, -39, -37, 32};
//+
Plane Surface(84) = {92};
//+
Line Loop(94) = {199, -210, -37, 209};
//+
Plane Surface(85) = {94};
//+
Line Loop(95) = {206, 215, -39, 210};
//+
Plane Surface(86) = {95};
//+
Line Loop(96) = {207, 193, 38, -215};
//+
Plane Surface(87) = {96};
//+
Line Loop(97) = {198, -209, 32, -193};
//+
Plane Surface(88) = {97};
//+
Line Loop(98) = {192, 193, -208, -194};
//+
Plane Surface(89) = {98};
//+
Line Loop(99) = {194, 197, -178, -173};
//+
Plane Surface(90) = {99};
//+
Line Loop(100) = {178, 213, 209, 177};
//+
Plane Surface(91) = {100};
//+
Line Loop(101) = {209, -198, 193, -32};
//+
Plane Surface(92) = {101};
//+
Line Loop(102) = {192, 198, 177, -173};
//+
Plane Surface(93) = {102};
//+
Line Loop(103) = {197, 213, 32, -208};
//+
Plane Surface(94) = {103};
//+
Line Loop(104) = {177, 174, 175, 176};
//+
Plane Surface(95) = {104};
//+
Line Loop(105) = {179, 180, -214, -213};
//+
Plane Surface(96) = {105};
//+
Line Loop(106) = {178, 179, -185, -174};
//+
Plane Surface(97) = {106};
//+
Line Loop(107) = {185, 180, 181, -175};
//+
Plane Surface(98) = {107};
//+
Line Loop(108) = {214, 181, 176, -209};
//+
Plane Surface(99) = {108};
//+
Line Loop(109) = {213, 209, 177, 178};
//+
Plane Surface(105) = {109};
//+
Line Loop(110) = {199, 200, 205, 176};
//+
Plane Surface(100) = {110};
//+
Line Loop(111) = {204, -211, -37, 214};
//+
Plane Surface(101) = {111};
//+
Line Loop(112) = {204, -201, 205, -181};
//+
Plane Surface(102) = {112};
//+
Line Loop(113) = {199, -210, -37, 209};
//+
Plane Surface(103) = {113};
//+
Line Loop(114) = {210, 200, 201, -211};
//+
Plane Surface(104) = {114};
//+
Line Loop(115) = {184, 180, 186, -187};
//+
Plane Surface(106) = {115};
//+
Line Loop(116) = {182, 183, 184, -179};
//+
Plane Surface(107) = {116};
//+
Line Loop(117) = {182, 189, 31, -213};
//+
Plane Surface(108) = {117};
//+
Line Loop(118) = {31, 214, 186, -188};
//+
Plane Surface(109) = {118};
//+
Line Loop(119) = {183, 187, -188, -189};
//+
Plane Surface(110) = {119};
//+
Line Loop(120) = {196, 189, 34, 190};
//+
Plane Surface(111) = {120};
//+
Line Loop(121) = {33, -32, -31, 34};
//+
Plane Surface(112) = {121};
//+
Line Loop(123) = {197, 182, -196, -195};
//+
Plane Surface(113) = {123};
//+
Line Loop(124) = {190, -195, 208, -33};
//+
Plane Surface(114) = {124};
//+
Line Loop(125) = {203, -188, 35, 212};
//+
Plane Surface(115) = {125};
//+
Line Loop(126) = {203, -186, 204, 202};
//+
Plane Surface(116) = {126};
//+
Line Loop(127) = {37, -36, -35, 31};
//+
Plane Surface(117) = {127};
//+
Line Loop(129) = {36, 211, 202, -212};
//+
Plane Surface(118) = {129};
//+
Surface Loop(18) = {83, 87, 88, 85, 21, 86};
//+
Volume(14) = {18};
//+
Surface Loop(19) = {90, 89, 93, 94, 88, 91};
//+
Volume(15) = {19};
//+
Surface Loop(20) = {98, 97, 95, 91, 96, 99};
//+
Volume(16) = {20};
//+
Surface Loop(21) = {102, 104, 100, 99, 85, 101};
//+
Volume(17) = {21};
//+
Surface Loop(22) = {115, 116, 118, 109, 20, 101};
//+
Volume(18) = {22};
//+
Surface Loop(23) = {110, 107, 106, 96, 108, 109};
//+
Volume(19) = {23};
//+
Surface Loop(24) = {114, 111, 113, 94, 19, 108};
//+
Volume(20) = {24};

//+
Physical Volume("acoustic",1) = {13};
//+
Physical Volume("PML",2) = {14, 17, 15, 20, 19, 16, 18};
//+
Physical Surface("surf",3) = {14};
