function obj = D_func_hand_eye_new(in1,in2,in3,in4)
coef_f0_q_sym1 = in1(:,1);
coef_f0_q_sym2 = in1(:,2);
coef_f0_q_sym3 = in1(:,3);
coef_f0_q_sym4 = in1(:,4);
coef_f0_q_sym5 = in1(:,5);
coef_f0_q_sym6 = in1(:,6);
coef_f0_q_sym7 = in1(:,7);
coef_f0_q_sym8 = in1(:,8);
coef_f0_q_sym9 = in1(:,9);
coef_f0_q_sym10 = in1(:,10);
coef_f1_q_sym1 = in2(:,1);
coef_f0_q_sym12 = in1(:,12);
coef_f1_q_sym2 = in2(:,2);
coef_f0_q_sym13 = in1(:,13);
coef_f1_q_sym3 = in2(:,3);
coef_f0_q_sym14 = in1(:,14);
coef_f1_q_sym4 = in2(:,4);
coef_f0_q_sym15 = in1(:,15);
coef_f1_q_sym5 = in2(:,5);
coef_f0_q_sym16 = in1(:,16);
coef_f1_q_sym6 = in2(:,6);
coef_f0_q_sym17 = in1(:,17);
coef_f1_q_sym7 = in2(:,7);
coef_f1_q_sym8 = in2(:,8);
coef_f0_q_sym19 = in1(:,19);
coef_f1_q_sym9 = in2(:,9);
coef_f0_q_sym20 = in1(:,20);
coef_f0_q_sym21 = in1(:,21);
coef_f2_q_sym1 = in3(:,1);
coef_f2_q_sym2 = in3(:,2);
coef_f0_q_sym23 = in1(:,23);
coef_f2_q_sym3 = in3(:,3);
coef_f2_q_sym4 = in3(:,4);
coef_f2_q_sym5 = in3(:,5);
coef_f2_q_sym6 = in3(:,6);
coef_f2_q_sym7 = in3(:,7);
coef_f2_q_sym8 = in3(:,8);
coef_f2_q_sym9 = in3(:,9);
coef_f3_q_sym1 = in4(:,1);
coef_f3_q_sym2 = in4(:,2);
coef_f3_q_sym3 = in4(:,3);
coef_f3_q_sym4 = in4(:,4);
coef_f3_q_sym5 = in4(:,5);
coef_f3_q_sym6 = in4(:,6);
coef_f3_q_sym7 = in4(:,7);
coef_f3_q_sym8 = in4(:,8);
coef_f3_q_sym9 = in4(:,9);
coef_f1_q_sym10 = in2(:,10);
coef_f1_q_sym12 = in2(:,12);
coef_f1_q_sym13 = in2(:,13);
coef_f1_q_sym14 = in2(:,14);
coef_f1_q_sym15 = in2(:,15);
coef_f1_q_sym16 = in2(:,16);
coef_f1_q_sym17 = in2(:,17);
coef_f1_q_sym19 = in2(:,19);
coef_f1_q_sym20 = in2(:,20);
coef_f1_q_sym21 = in2(:,21);
coef_f1_q_sym23 = in2(:,23);
coef_f2_q_sym10 = in3(:,10);
coef_f2_q_sym12 = in3(:,12);
coef_f2_q_sym13 = in3(:,13);
coef_f2_q_sym14 = in3(:,14);
coef_f2_q_sym15 = in3(:,15);
coef_f2_q_sym16 = in3(:,16);
coef_f2_q_sym17 = in3(:,17);
coef_f2_q_sym19 = in3(:,19);
coef_f2_q_sym20 = in3(:,20);
coef_f2_q_sym21 = in3(:,21);
coef_f2_q_sym23 = in3(:,23);
coef_f3_q_sym10 = in4(:,10);
coef_f3_q_sym12 = in4(:,12);
coef_f3_q_sym13 = in4(:,13);
coef_f3_q_sym14 = in4(:,14);
coef_f3_q_sym15 = in4(:,15);
coef_f3_q_sym16 = in4(:,16);
coef_f3_q_sym17 = in4(:,17);
coef_f3_q_sym19 = in4(:,19);
coef_f3_q_sym20 = in4(:,20);
coef_f3_q_sym21 = in4(:,21);
coef_f3_q_sym23 = in4(:,23);
t2 = coef_f1_q_sym1.*2.0;
t3 = coef_f2_q_sym1.*2.0;
t4 = coef_f3_q_sym1.*2.0;
t5 = -coef_f0_q_sym2;
t6 = -coef_f0_q_sym3;
t7 = -coef_f0_q_sym4;
t8 = -coef_f1_q_sym1;
t10 = -coef_f2_q_sym1;
t12 = -coef_f3_q_sym1;
t9 = -t2;
t11 = -t3;
t13 = -t4;
obj = reshape([coef_f0_q_sym1-coef_f1_q_sym2,-coef_f2_q_sym2,-coef_f3_q_sym2,-coef_f1_q_sym3,coef_f0_q_sym1-coef_f2_q_sym3,-coef_f3_q_sym3,-coef_f1_q_sym4,-coef_f2_q_sym4,coef_f0_q_sym1-coef_f3_q_sym4,coef_f0_q_sym5-coef_f1_q_sym12,-coef_f2_q_sym12,-coef_f3_q_sym12,coef_f0_q_sym6-coef_f1_q_sym13,coef_f0_q_sym5-coef_f2_q_sym13,-coef_f3_q_sym13,coef_f0_q_sym7-coef_f1_q_sym14,-coef_f2_q_sym14,coef_f0_q_sym5-coef_f3_q_sym14,coef_f0_q_sym8-coef_f1_q_sym15,coef_f0_q_sym6-coef_f2_q_sym15,-coef_f3_q_sym15,coef_f0_q_sym9-coef_f1_q_sym16,coef_f0_q_sym7-coef_f2_q_sym16,coef_f0_q_sym6-coef_f3_q_sym16,coef_f0_q_sym10-coef_f1_q_sym17,-coef_f2_q_sym17,coef_f0_q_sym7-coef_f3_q_sym17,-coef_f1_q_sym19,coef_f0_q_sym8-coef_f2_q_sym19,-coef_f3_q_sym19,-coef_f1_q_sym20,coef_f0_q_sym9-coef_f2_q_sym20,coef_f0_q_sym8-coef_f3_q_sym20,-coef_f1_q_sym21,coef_f0_q_sym10-coef_f2_q_sym21,coef_f0_q_sym9-coef_f3_q_sym21,-coef_f1_q_sym23,-coef_f2_q_sym23,coef_f0_q_sym10-coef_f3_q_sym23,coef_f0_q_sym12+coef_f1_q_sym5+t5+t8,coef_f2_q_sym5+t10,coef_f3_q_sym5+t12,coef_f0_q_sym13+coef_f1_q_sym6+t6,coef_f0_q_sym12+coef_f2_q_sym6+t5,coef_f3_q_sym6,coef_f0_q_sym14+coef_f1_q_sym7+t7,coef_f2_q_sym7,coef_f0_q_sym12+coef_f3_q_sym7+t5,coef_f0_q_sym15+coef_f1_q_sym5+coef_f1_q_sym8+t5+t9,coef_f0_q_sym13+coef_f2_q_sym5+coef_f2_q_sym8+t6+t11,coef_f3_q_sym5+coef_f3_q_sym8+t13,coef_f0_q_sym16+coef_f1_q_sym9,coef_f0_q_sym14+coef_f2_q_sym9+t7,coef_f0_q_sym13+coef_f3_q_sym9+t6,coef_f1_q_sym5+coef_f0_q_sym17+coef_f1_q_sym10+t5+t9,coef_f2_q_sym5+coef_f2_q_sym10+t11,coef_f0_q_sym14+coef_f3_q_sym5+coef_f3_q_sym10+t7+t13,coef_f1_q_sym6+coef_f0_q_sym19+t6,coef_f0_q_sym15+coef_f2_q_sym6+t5,coef_f3_q_sym6,coef_f1_q_sym7+coef_f0_q_sym20+t7,coef_f0_q_sym16+coef_f2_q_sym7,coef_f0_q_sym15+coef_f3_q_sym7+t5,coef_f1_q_sym6+coef_f0_q_sym21+t6,coef_f0_q_sym17+coef_f2_q_sym6+t5,coef_f0_q_sym16+coef_f3_q_sym6,coef_f1_q_sym7+coef_f0_q_sym23+t7,coef_f2_q_sym7,coef_f0_q_sym17+coef_f3_q_sym7+t5,coef_f1_q_sym8+t8,coef_f0_q_sym19+coef_f2_q_sym8+t6+t10,coef_f3_q_sym8+t12,coef_f1_q_sym9,coef_f0_q_sym20+coef_f2_q_sym9+t7,coef_f0_q_sym19+coef_f3_q_sym9+t6,coef_f1_q_sym8+coef_f1_q_sym10+t9,coef_f0_q_sym21+coef_f2_q_sym8+coef_f2_q_sym10+t6+t11,coef_f0_q_sym20+coef_f3_q_sym8+coef_f3_q_sym10+t7+t13,coef_f1_q_sym9,coef_f0_q_sym23+coef_f2_q_sym9+t7,coef_f0_q_sym21+coef_f3_q_sym9+t6,coef_f1_q_sym10+t8,coef_f2_q_sym10+t10,coef_f0_q_sym23+coef_f3_q_sym10+t7+t12],[3,28]);
