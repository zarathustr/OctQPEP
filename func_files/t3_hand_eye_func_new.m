function obj = t3_hand_eye_func_new(in1,in2,in3)
coefs_tq1_1 = in2(1);
coefs_tq1_2 = in2(4);
coefs_tq1_3 = in2(7);
coefs_tq1_4 = in2(10);
coefs_tq1_5 = in2(13);
coefs_tq1_6 = in2(16);
coefs_tq1_7 = in2(19);
coefs_tq1_8 = in2(22);
coefs_tq1_9 = in2(25);
coefs_tq2_1 = in2(2);
coefs_tq2_2 = in2(5);
coefs_tq2_3 = in2(8);
coefs_tq2_4 = in2(11);
coefs_tq2_5 = in2(14);
coefs_tq2_6 = in2(17);
coefs_tq2_7 = in2(20);
coefs_tq2_8 = in2(23);
coefs_tq2_9 = in2(26);
coefs_tq3_1 = in2(3);
coefs_tq3_2 = in2(6);
coefs_tq3_3 = in2(9);
coefs_tq3_4 = in2(12);
coefs_tq3_5 = in2(15);
coefs_tq3_6 = in2(18);
coefs_tq3_7 = in2(21);
coefs_tq3_8 = in2(24);
coefs_tq3_9 = in2(27);
coefs_tq1_10 = in2(28);
coefs_tq1_11 = in2(31);
coefs_tq2_10 = in2(29);
coefs_tq2_11 = in2(32);
coefs_tq3_10 = in2(30);
coefs_tq3_11 = in2(33);
pinvG3_1 = in1(3);
pinvG3_2 = in1(6);
pinvG3_3 = in1(9);
q0 = in3(1,:);
q1 = in3(2,:);
q2 = in3(3,:);
q3 = in3(4,:);
obj = q0.^2.*(coefs_tq1_1.*pinvG3_1+coefs_tq2_1.*pinvG3_2+coefs_tq3_1.*pinvG3_3)+q1.^2.*(coefs_tq1_5.*pinvG3_1+coefs_tq2_5.*pinvG3_2+coefs_tq3_5.*pinvG3_3)+q2.^2.*(coefs_tq1_8.*pinvG3_1+coefs_tq2_8.*pinvG3_2+coefs_tq3_8.*pinvG3_3)+q3.^2.*(coefs_tq1_10.*pinvG3_1+coefs_tq2_10.*pinvG3_2+coefs_tq3_10.*pinvG3_3)+coefs_tq1_11.*pinvG3_1+coefs_tq2_11.*pinvG3_2+coefs_tq3_11.*pinvG3_3+q0.*q1.*(coefs_tq1_2.*pinvG3_1+coefs_tq2_2.*pinvG3_2+coefs_tq3_2.*pinvG3_3)+q0.*q2.*(coefs_tq1_3.*pinvG3_1+coefs_tq2_3.*pinvG3_2+coefs_tq3_3.*pinvG3_3)+q0.*q3.*(coefs_tq1_4.*pinvG3_1+coefs_tq2_4.*pinvG3_2+coefs_tq3_4.*pinvG3_3)+q1.*q2.*(coefs_tq1_6.*pinvG3_1+coefs_tq2_6.*pinvG3_2+coefs_tq3_6.*pinvG3_3)+q1.*q3.*(coefs_tq1_7.*pinvG3_1+coefs_tq2_7.*pinvG3_2+coefs_tq3_7.*pinvG3_3)+q2.*q3.*(coefs_tq1_9.*pinvG3_1+coefs_tq2_9.*pinvG3_2+coefs_tq3_9.*pinvG3_3);