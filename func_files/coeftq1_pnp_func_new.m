function obj = coeftq1_pnp_func_new(in1)
coef_J11 = in1(:,11);
coef_J20 = in1(:,20);
coef_J30 = in1(:,30);
coef_J26 = in1(:,26);
coef_J45 = in1(:,45);
coef_J63 = in1(:,63);
coef_J55 = in1(:,55);
coef_J39 = in1(:,39);
coef_J49 = in1(:,49);
coef_J59 = in1(:,59);
obj = [coef_J11,coef_J20,coef_J26,coef_J30,coef_J39,coef_J45,coef_J49,coef_J55,coef_J59,coef_J63];
