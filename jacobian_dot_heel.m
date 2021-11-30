function dJheel = jacobian_dot_heel(in1,in2)
%JACOBIAN_DOT_HEEL
%    DJHEEL = JACOBIAN_DOT_HEEL(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Nov-2021 15:09:52

dth1 = in1(9,:);
dth2 = in1(10,:);
dth3 = in1(11,:);
l_AC = in2(23,:);
l_DE = in2(24,:);
l_IG = in2(25,:);
l_OB = in2(22,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = cos(t4);
t6 = sin(t4);
t7 = t4+th3;
t8 = cos(t7);
t9 = sin(t7);
t10 = l_AC.*t5;
t11 = l_AC.*t6;
t12 = l_IG.*t8;
t13 = l_IG.*t9;
t14 = dth3.*t12;
t15 = dth3.*t13;
t16 = -t12;
t17 = -t13;
t18 = -t14;
t19 = t10+t16;
t20 = t11+t17;
t21 = dth2.*t19;
t22 = dth2.*t20;
t23 = -t22;
dJheel = reshape([0.0,0.0,0.0,0.0,t15+t23-dth1.*(t20+l_DE.*t3+l_OB.*t3),t18+t21+dth1.*(t19+l_DE.*t2+l_OB.*t2),t15+t23-dth1.*t20,t18+t21+dth1.*t19,t15+dth1.*t13+dth2.*t13,t18+dth1.*t16+dth2.*t16,0.0,0.0],[2,6]);
