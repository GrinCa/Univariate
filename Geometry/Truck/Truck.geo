// Gmsh project created on Fri Sep  7 09:27:28 2018
//+

SetFactory("OpenCASCADE");


L=1;
Lpml = 1.5;
coeff_reduction = 4;

sizemesh = 1.5;


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
Surface Loop(15) = {1, 3, 5, 4, 2, 6};
//+
Surface Loop(16) = {6, 1, 3, 5, 4, 2};
//+
Surface Loop(17) = {21, 23, 19, 24, 22, 20};




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
Point(120) = {7, 6, 1.5, sizemesh};
//+
Point(121) = {7, 6, 2.5, sizemesh};
//+
Point(122) = {8, 6, 1.5, sizemesh};
//+
Point(123) = {8, 6, 2.5, sizemesh};
//+
Point(124) = {-1, 6, 0, sizemesh};
//+
Point(125) = {-1, 6, 4.8, sizemesh};
//+
Point(126) = {12, 6, 0, sizemesh};
//+
Point(127) = {12, 6, 4.8, sizemesh};
//+
Point(128) = {-1, -1, 0, sizemesh};
//+
Point(129) = {-1, -1, 4.8, sizemesh};
//+
Point(130) = {12, -1, 0, sizemesh};
//+
Point(131) = {12, -1, 4.8, sizemesh};



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
Line(192) = {102, 103};
//+
Line(194) = {102, 119};
//+
Line(195) = {119, 101};
//+
Line(196) = {101, 107};
//+
Line(197) = {119, 108};
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
Line(216) = {123, 122};
//+
Line(217) = {122, 120};
//+
Line(218) = {120, 121};
//+
Line(219) = {121, 123};
//+
Line(220) = {116, 131};
//+
Line(221) = {131, 114};
//+
Line(222) = {131, 127};
//+
Line(223) = {127, 118};
//+
Line(224) = {115, 130};
//+
Line(225) = {130, 126};
//+
Line(226) = {126, 124};
//+
Line(227) = {124, 101};
//+
Line(228) = {124, 125};
//+
Line(229) = {125, 103};
//+
Line(230) = {125, 119};
//+
Line(231) = {125, 129};
//+
Line(232) = {129, 113};
//+
Line(233) = {129, 108};
//+
Line(234) = {129, 111};
//+
Line(235) = {103, 113};
//+
Line(236) = {113, 114};
//+
Line(237) = {127, 125};
//+
Line(238) = {129, 131};
//+
Line(239) = {130, 131};
//+
Line(240) = {127, 126};
//+
Line(241) = {110, 128};
//+
Line(242) = {128, 107};
//+
Line(243) = {128, 124};
//+
Line(244) = {130, 128};
//+
Line(245) = {129, 128};


//+
Line Loop(49) = {219, 216, 217, 218};
//+
Plane Surface(49) = {49};
//+
Line Loop(50) = {237, -228, -226, -240};
//+
Line Loop(51) = {219, 216, 217, 218};
//+
Plane Surface(50) = {50, 51};
//+
Line Loop(52) = {222, 240, -225, 239};
//+
Plane Surface(51) = {52};

//+
Line Loop(53) = {238, -239, 244, -245};
//+
Plane Surface(52) = {53};
//+
Line Loop(54) = {243, 228, 231, 245};
//+
Plane Surface(53) = {54};
//+
Line Loop(55) = {237, 231, 238, 222};
//+
Plane Surface(54) = {55};
//+
Line Loop(56) = {225, 226, -243, -244};
//+
Plane Surface(55) = {56};
//+
Surface Loop(18) = {52, 54, 50, 53, 55, 51, 49};
//+
Surface Loop(19) = {51, 54, 53, 55, 52, 50, 49};
//+
Volume(21) = {19};

BooleanDifference{ Volume{21} ; Delete; }{ Volume{4}; Volume{9}; Volume{7}; Volume{6}; Volume{8}; Volume{11}; Volume{10}; Volume{5}; Volume{2}; Volume{20}; Volume{1}; Delete; }



//+
Line Loop(131) = {230, -194, 192, -229};
//+
Plane Surface(122) = {131};
//+
Line Loop(132) = {197, -178, -173, 194};
//+
Plane Surface(123) = {132};
//+
Line Loop(133) = {178, -233, 232, 177};
//+
Plane Surface(124) = {133};
//+
Line Loop(134) = {232, -235, -229, 246};
//+
Plane Surface(125) = {134};
//+
Line Loop(135) = {230, 197, -233, -246};
//+
Plane Surface(126) = {135};
//+
Line Loop(136) = {177, -173, 192, 235};
//+
Plane Surface(127) = {136};
//+
Line Loop(137) = {174, 185, -179, -178};
//+
Plane Surface(128) = {137};
//+
Line Loop(138) = {180, 181, -175, 185};
//+
Plane Surface(129) = {138};
//+
Line Loop(139) = {232, -176, -181, -234};
//+
Plane Surface(130) = {139};
//+
Line Loop(140) = {174, 175, 176, 177};
//+
Plane Surface(131) = {140};
//+
Line Loop(141) = {180, -234, 233, 179};
//+
Plane Surface(132) = {141};
//+
Line Loop(142) = {236, 200, 205, 176};
//+
Plane Surface(133) = {142};
//+
Line Loop(143) = {181, -205, 201, -204};
//+
Plane Surface(134) = {143};
//+
Line Loop(144) = {234, 204, 220, -247};
//+
Plane Surface(135) = {144};
//+
Line Loop(145) = {247, 221, -236, -232};
//+
Plane Surface(136) = {145};
//+
Line Loop(146) = {221, 200, 201, 220};
//+
Plane Surface(137) = {146};
//+
Line Loop(147) = {206, -223, -243, 221};
//+
Plane Surface(138) = {147};
//+
Line Loop(148) = {237, 229, -207, -223};
//+
Plane Surface(139) = {148};
//+
Line Loop(149) = {235, 236, 206, 207};
//+
Plane Surface(140) = {149};
//+
Line Loop(150) = {245, -247, 256, -248};
//+
Plane Surface(141) = {150};

//+
Line Loop(152) = {248, -241, -203, 224};
//+
Plane Surface(142) = {152};
//+
Line Loop(153) = {245, -220, 202, 224};
//+
Plane Surface(143) = {153};
//+
Line Loop(155) = {241, -256, 234, 186};
//+
Plane Surface(144) = {155};
//+
Line Loop(156) = {203, -186, 204, 202};
//+
Plane Surface(145) = {156};
//+
Line Loop(157) = {242, 183, 187, 241};
//+
Plane Surface(146) = {157};
//+
Line Loop(158) = {184, 180, 186, -187};
//+
Plane Surface(147) = {158};
//+
Line Loop(159) = {179, -184, -183, -182};
//+
Plane Surface(148) = {159};
//+
Line Loop(160) = {242, -182, -233, 256};
//+
Plane Surface(149) = {160};
//+
Line Loop(161) = {249, 227, 196, -242};
//+
Plane Surface(150) = {161};
//+
Line Loop(162) = {227, -195, -230, -228};
//+
Plane Surface(151) = {162};
//+
Line Loop(163) = {196, -182, -197, 195};
//+
Plane Surface(152) = {163};
//+
Surface Loop(19) = {150, 151, 152, 62, 126, 149};
//+
Volume(22) = {19};
//+
Surface Loop(20) = {146, 148, 147, 149, 144, 132};
//+
Volume(23) = {20};
//+
Surface Loop(21) = {142, 145, 143, 61, 135, 144};
//+
Volume(24) = {21};
//+
Surface Loop(22) = {134, 133, 137, 136, 130, 135};
//+
Volume(25) = {22};
//+
Surface Loop(23) = {132, 124, 128, 131, 129, 130};
//+
Volume(26) = {23};
//+
Surface Loop(24) = {126, 123, 127, 122, 125, 124};
//+
Volume(27) = {24};
//+
Surface Loop(25) = {58, 139, 140, 138, 125, 136};
//+
Volume(28) = {25};
//+
Characteristic Length {113, 112, 106, 109, 108, 129, 105, 111, 131, 116, 117, 114, 118, 127, 103, 125, 119, 102, 101, 124, 128, 107, 104, 110, 130, 115, 126, 192, 191, 189, 190, 180, 193, 188, 187, 182, 181, 184, 186, 185, 183, 201, 199, 15, 13, 17, 19, 20, 18, 167, 168, 164, 165, 179, 178, 177, 176, 150, 149, 145, 144, 147, 146, 148, 161, 163, 162, 160, 159, 157, 158, 156, 169, 154, 155, 153, 152, 151, 172, 171, 175, 174, 173, 166, 5, 7, 3, 1} = sizemesh;
//+
Characteristic Length {137, 136, 140, 141, 132, 133, 142, 143, 138, 139, 135, 134, 2, 4, 196,200,203,197,202,194,198,195} = sizemesh/coeff_reduction;


//+
Physical Volume("PML",2) = {28, 25, 27, 26, 23, 24, 22};
//+
Physical Volume("acoustic",1) = {21};
//+
Physical Surface("surface",3) = {63};
