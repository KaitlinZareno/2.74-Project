function Corr_Joint_Sp = Corr_leg(in1,in2)
%CORR_LEG
%    CORR_JOINT_SP = CORR_LEG(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    09-Dec-2021 15:44:45

c5 = in2(20,:);
dth1 = in1(11,:);
dth2 = in1(12,:);
dth3 = in1(13,:);
dth21 = in1(14,:);
dth22 = in1(15,:);
dth23 = in1(16,:);
l_AC = in2(23,:);
l_A_m3 = in2(18,:);
l_B_m2 = in2(17,:);
l_C_m4 = in2(19,:);
l_DE = in2(24,:);
l_OA = in2(21,:);
l_OB = in2(22,:);
l_O_m1 = in2(16,:);
m1 = in2(2,:);
m2 = in2(3,:);
m3 = in2(4,:);
m4 = in2(5,:);
m5 = in2(6,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
th21 = in1(6,:);
th22 = in1(7,:);
th23 = in1(8,:);
t2 = cos(th1);
t3 = cos(th21);
t4 = sin(th1);
t5 = sin(th2);
t6 = sin(th3);
t7 = sin(th21);
t8 = sin(th22);
t9 = sin(th23);
t10 = th1+th2;
t11 = th2+th3;
t12 = th21+th22;
t13 = th22+th23;
t14 = dth1.^2;
t15 = dth2.^2;
t16 = dth3.^2;
t17 = dth21.^2;
t18 = dth22.^2;
t19 = dth23.^2;
t20 = cos(t10);
t21 = cos(t12);
t22 = sin(t10);
t23 = sin(t11);
t24 = sin(t12);
t25 = sin(t13);
t26 = t10+th3;
t27 = t12+th23;
t32 = c5.*l_AC.*m5.*t6.*t16;
t33 = c5.*l_AC.*m5.*t9.*t19;
t34 = c5.*dth1.*dth3.*l_AC.*m5.*t6.*2.0;
t35 = c5.*dth2.*dth3.*l_AC.*m5.*t6.*2.0;
t36 = c5.*dth21.*dth23.*l_AC.*m5.*t9.*2.0;
t37 = c5.*dth22.*dth23.*l_AC.*m5.*t9.*2.0;
t28 = cos(t26);
t29 = cos(t27);
t30 = sin(t26);
t31 = sin(t27);
t38 = -t34;
t39 = -t35;
t40 = -t36;
t41 = -t37;
t42 = -t32;
t43 = -t33;
Corr_Joint_Sp = [-c5.*m5.*t14.*t30-c5.*m5.*t15.*t30-c5.*m5.*t16.*t30-c5.*m5.*t17.*t31-c5.*m5.*t18.*t31-c5.*m5.*t19.*t31-l_AC.*m4.*t14.*t22-l_AC.*m4.*t15.*t22-l_AC.*m5.*t14.*t22-l_AC.*m5.*t15.*t22-l_AC.*m4.*t17.*t24-l_AC.*m4.*t18.*t24-l_AC.*m5.*t17.*t24-l_AC.*m5.*t18.*t24-l_A_m3.*m3.*t14.*t22-l_A_m3.*m3.*t15.*t22-l_A_m3.*m3.*t17.*t24-l_A_m3.*m3.*t18.*t24-l_B_m2.*m2.*t14.*t22-l_B_m2.*m2.*t15.*t22-l_B_m2.*m2.*t17.*t24-l_B_m2.*m2.*t18.*t24-l_C_m4.*m4.*t4.*t14-l_C_m4.*m4.*t7.*t17-l_DE.*m5.*t4.*t14-l_DE.*m5.*t7.*t17-l_OA.*m3.*t4.*t14-l_OB.*m2.*t4.*t14-l_OA.*m4.*t4.*t14-l_OB.*m5.*t4.*t14-l_OA.*m3.*t7.*t17-l_OB.*m2.*t7.*t17-l_OA.*m4.*t7.*t17-l_OB.*m5.*t7.*t17-l_O_m1.*m1.*t4.*t14-l_O_m1.*m1.*t7.*t17-c5.*dth1.*dth2.*m5.*t30.*2.0-c5.*dth1.*dth3.*m5.*t30.*2.0-c5.*dth2.*dth3.*m5.*t30.*2.0-c5.*dth21.*dth22.*m5.*t31.*2.0-c5.*dth21.*dth23.*m5.*t31.*2.0-c5.*dth22.*dth23.*m5.*t31.*2.0-dth1.*dth2.*l_AC.*m4.*t22.*2.0-dth1.*dth2.*l_AC.*m5.*t22.*2.0-dth21.*dth22.*l_AC.*m4.*t24.*2.0-dth21.*dth22.*l_AC.*m5.*t24.*2.0-dth1.*dth2.*l_A_m3.*m3.*t22.*2.0-dth21.*dth22.*l_A_m3.*m3.*t24.*2.0-dth1.*dth2.*l_B_m2.*m2.*t22.*2.0-dth21.*dth22.*l_B_m2.*m2.*t24.*2.0;c5.*m5.*t14.*t28+c5.*m5.*t15.*t28+c5.*m5.*t16.*t28+c5.*m5.*t17.*t29+c5.*m5.*t18.*t29+c5.*m5.*t19.*t29+l_AC.*m4.*t14.*t20+l_AC.*m4.*t15.*t20+l_AC.*m5.*t14.*t20+l_AC.*m5.*t15.*t20+l_AC.*m4.*t17.*t21+l_AC.*m4.*t18.*t21+l_AC.*m5.*t17.*t21+l_AC.*m5.*t18.*t21+l_A_m3.*m3.*t14.*t20+l_A_m3.*m3.*t15.*t20+l_A_m3.*m3.*t17.*t21+l_A_m3.*m3.*t18.*t21+l_B_m2.*m2.*t14.*t20+l_B_m2.*m2.*t15.*t20+l_B_m2.*m2.*t17.*t21+l_B_m2.*m2.*t18.*t21+l_C_m4.*m4.*t2.*t14+l_C_m4.*m4.*t3.*t17+l_DE.*m5.*t2.*t14+l_DE.*m5.*t3.*t17+l_OA.*m3.*t2.*t14+l_OB.*m2.*t2.*t14+l_OA.*m4.*t2.*t14+l_OB.*m5.*t2.*t14+l_OA.*m3.*t3.*t17+l_OB.*m2.*t3.*t17+l_OA.*m4.*t3.*t17+l_OB.*m5.*t3.*t17+l_O_m1.*m1.*t2.*t14+l_O_m1.*m1.*t3.*t17+c5.*dth1.*dth2.*m5.*t28.*2.0+c5.*dth1.*dth3.*m5.*t28.*2.0+c5.*dth2.*dth3.*m5.*t28.*2.0+c5.*dth21.*dth22.*m5.*t29.*2.0+c5.*dth21.*dth23.*m5.*t29.*2.0+c5.*dth22.*dth23.*m5.*t29.*2.0+dth1.*dth2.*l_AC.*m4.*t20.*2.0+dth1.*dth2.*l_AC.*m5.*t20.*2.0+dth21.*dth22.*l_AC.*m4.*t21.*2.0+dth21.*dth22.*l_AC.*m5.*t21.*2.0+dth1.*dth2.*l_A_m3.*m3.*t20.*2.0+dth21.*dth22.*l_A_m3.*m3.*t21.*2.0+dth1.*dth2.*l_B_m2.*m2.*t20.*2.0+dth21.*dth22.*l_B_m2.*m2.*t21.*2.0;t38+t39+t42-c5.*l_DE.*m5.*t15.*t23-c5.*l_DE.*m5.*t16.*t23-c5.*l_OB.*m5.*t15.*t23-c5.*l_OB.*m5.*t16.*t23-l_AC.*l_C_m4.*m4.*t5.*t15-l_AC.*l_DE.*m5.*t5.*t15-l_AC.*l_OA.*m4.*t5.*t15-l_AC.*l_OB.*m5.*t5.*t15-l_A_m3.*l_OA.*m3.*t5.*t15-l_B_m2.*l_OB.*m2.*t5.*t15-c5.*dth1.*dth2.*l_DE.*m5.*t23.*2.0-c5.*dth1.*dth3.*l_DE.*m5.*t23.*2.0-c5.*dth2.*dth3.*l_DE.*m5.*t23.*2.0-c5.*dth1.*dth2.*l_OB.*m5.*t23.*2.0-c5.*dth1.*dth3.*l_OB.*m5.*t23.*2.0-c5.*dth2.*dth3.*l_OB.*m5.*t23.*2.0-dth1.*dth2.*l_AC.*l_C_m4.*m4.*t5.*2.0-dth1.*dth2.*l_AC.*l_DE.*m5.*t5.*2.0-dth1.*dth2.*l_AC.*l_OA.*m4.*t5.*2.0-dth1.*dth2.*l_AC.*l_OB.*m5.*t5.*2.0-dth1.*dth2.*l_A_m3.*l_OA.*m3.*t5.*2.0-dth1.*dth2.*l_B_m2.*l_OB.*m2.*t5.*2.0;t38+t39+t42+c5.*l_DE.*m5.*t14.*t23+c5.*l_OB.*m5.*t14.*t23+l_AC.*l_C_m4.*m4.*t5.*t14+l_AC.*l_DE.*m5.*t5.*t14+l_AC.*l_OA.*m4.*t5.*t14+l_AC.*l_OB.*m5.*t5.*t14+l_A_m3.*l_OA.*m3.*t5.*t14+l_B_m2.*l_OB.*m2.*t5.*t14;c5.*m5.*(l_AC.*t6.*t14+l_AC.*t6.*t15+l_DE.*t14.*t23+l_OB.*t14.*t23+dth1.*dth2.*l_AC.*t6.*2.0);t40+t41+t43-c5.*l_DE.*m5.*t18.*t25-c5.*l_DE.*m5.*t19.*t25-c5.*l_OB.*m5.*t18.*t25-c5.*l_OB.*m5.*t19.*t25-l_AC.*l_C_m4.*m4.*t8.*t18-l_AC.*l_DE.*m5.*t8.*t18-l_AC.*l_OA.*m4.*t8.*t18-l_AC.*l_OB.*m5.*t8.*t18-l_A_m3.*l_OA.*m3.*t8.*t18-l_B_m2.*l_OB.*m2.*t8.*t18-c5.*dth21.*dth22.*l_DE.*m5.*t25.*2.0-c5.*dth21.*dth23.*l_DE.*m5.*t25.*2.0-c5.*dth22.*dth23.*l_DE.*m5.*t25.*2.0-c5.*dth21.*dth22.*l_OB.*m5.*t25.*2.0-c5.*dth21.*dth23.*l_OB.*m5.*t25.*2.0-c5.*dth22.*dth23.*l_OB.*m5.*t25.*2.0-dth21.*dth22.*l_AC.*l_C_m4.*m4.*t8.*2.0-dth21.*dth22.*l_AC.*l_DE.*m5.*t8.*2.0-dth21.*dth22.*l_AC.*l_OA.*m4.*t8.*2.0-dth21.*dth22.*l_AC.*l_OB.*m5.*t8.*2.0-dth21.*dth22.*l_A_m3.*l_OA.*m3.*t8.*2.0-dth21.*dth22.*l_B_m2.*l_OB.*m2.*t8.*2.0;t40+t41+t43+c5.*l_DE.*m5.*t17.*t25+c5.*l_OB.*m5.*t17.*t25+l_AC.*l_C_m4.*m4.*t8.*t17+l_AC.*l_DE.*m5.*t8.*t17+l_AC.*l_OA.*m4.*t8.*t17+l_AC.*l_OB.*m5.*t8.*t17+l_A_m3.*l_OA.*m3.*t8.*t17+l_B_m2.*l_OB.*m2.*t8.*t17;c5.*m5.*(l_AC.*t9.*t17+l_AC.*t9.*t18+l_DE.*t17.*t25+l_OB.*t17.*t25+dth21.*dth22.*l_AC.*t9.*2.0)];
