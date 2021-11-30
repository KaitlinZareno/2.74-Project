function Corr_Joint_Sp = Corr_leg(in1,in2)
%CORR_LEG
%    CORR_JOINT_SP = CORR_LEG(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Nov-2021 15:22:38

dth1 = in1(4,:);
dth2 = in1(5,:);
dths = in1(6,:);
l_AC = in2(26,:);
l_CE = in2(28,:);
l_DE = in2(27,:);
l_OB = in2(25,:);
l_m2 = in2(19,:);
l_m3 = in2(20,:);
l_m4 = in2(21,:);
l_ms = in2(22,:);
m0 = in2(1,:);
m2 = in2(3,:);
m3 = in2(4,:);
m4 = in2(5,:);
ms = in2(7,:);
th1 = in1(1,:);
th2 = in1(2,:);
ths = in1(3,:);
t2 = sin(th2);
t3 = dth1.^2;
t4 = dth2.^2;
t5 = dths.^2;
t6 = -ths;
t7 = t6+th1;
t8 = sin(t7);
t9 = t7+th2;
t10 = sin(t9);
t11 = l_AC.*l_ms.*ms.*t5.*t10;
t12 = -t11;
Corr_Joint_Sp = [t12+l_AC.*l_DE.*m0.*t2.*t4+l_AC.*l_DE.*m4.*t2.*t4+l_AC.*l_OB.*m0.*t2.*t4+l_AC.*l_m4.*m4.*t2.*t4+l_CE.*l_m2.*m2.*t2.*t4+l_DE.*l_m3.*m3.*t2.*t4+l_AC.*l_DE.*ms.*t2.*t4+l_AC.*l_OB.*ms.*t2.*t4+l_DE.*l_ms.*ms.*t5.*t8+l_OB.*l_ms.*ms.*t5.*t8+dth1.*dth2.*l_AC.*l_DE.*m0.*t2.*2.0+dth1.*dth2.*l_AC.*l_DE.*m4.*t2.*2.0+dth1.*dth2.*l_AC.*l_OB.*m0.*t2.*2.0+dth1.*dth2.*l_AC.*l_m4.*m4.*t2.*2.0+dth1.*dth2.*l_CE.*l_m2.*m2.*t2.*2.0+dth1.*dth2.*l_DE.*l_m3.*m3.*t2.*2.0+dth1.*dth2.*l_AC.*l_DE.*ms.*t2.*2.0+dth1.*dth2.*l_AC.*l_OB.*ms.*t2.*2.0;t12-l_AC.*l_DE.*m0.*t2.*t3-l_AC.*l_DE.*m4.*t2.*t3-l_AC.*l_OB.*m0.*t2.*t3-l_AC.*l_m4.*m4.*t2.*t3-l_CE.*l_m2.*m2.*t2.*t3-l_DE.*l_m3.*m3.*t2.*t3-l_AC.*l_DE.*ms.*t2.*t3-l_AC.*l_OB.*ms.*t2.*t3;l_ms.*ms.*(l_AC.*t3.*t10+l_AC.*t4.*t10-l_DE.*t3.*t8-l_OB.*t3.*t8+dth1.*dth2.*l_AC.*t10.*2.0)];
