function obj = coef_Jacob4_qt_pnp_func_new(in1)
coef_J4 = in1(:,4);
coef_J7 = in1(:,7);
coef_J9 = in1(:,9);
coef_J10 = in1(:,10);
coef_J30 = in1(:,30);
coef_J31 = in1(:,31);
coef_J32 = in1(:,32);
coef_J50 = in1(:,50);
coef_J24 = in1(:,24);
coef_J51 = in1(:,51);
coef_J60 = in1(:,60);
coef_J16 = in1(:,16);
coef_J25 = in1(:,25);
coef_J43 = in1(:,43);
coef_J61 = in1(:,61);
coef_J35 = in1(:,35);
coef_J44 = in1(:,44);
coef_J53 = in1(:,53);
coef_J62 = in1(:,62);
coef_J18 = in1(:,18);
coef_J54 = in1(:,54);
coef_J63 = in1(:,63);
coef_J19 = in1(:,19);
coef_J37 = in1(:,37);
coef_J64 = in1(:,64);
coef_J29 = in1(:,29);
coef_J38 = in1(:,38);
coef_J65 = in1(:,65);
coef_J48 = in1(:,48);
coef_J49 = in1(:,49);
coef_J58 = in1(:,58);
coef_J59 = in1(:,59);
obj = [-coef_J4+coef_J16,coef_J18,coef_J19.*2.0,-coef_J4+coef_J24,coef_J25.*2.0,-coef_J4+coef_J29.*3.0,coef_J30,coef_J31,coef_J32,coef_J4,-coef_J7+coef_J35,-coef_J9+coef_J37,coef_J10.*-2.0+coef_J38.*2.0,-coef_J7+coef_J43,coef_J44.*2.0,-coef_J7+coef_J48.*3.0,coef_J49,coef_J50,coef_J51,coef_J7,-coef_J9+coef_J53,coef_J10.*-2.0+coef_J54.*2.0,-coef_J9+coef_J58.*3.0,coef_J59,coef_J60,coef_J61,coef_J9,coef_J10.*-2.0+coef_J62.*4.0,coef_J63.*2.0,coef_J64.*2.0,coef_J65.*2.0,coef_J10.*2.0];