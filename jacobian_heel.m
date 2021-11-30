function Jheel = jacobian_heel(in1,in2)
%JACOBIAN_HEEL
%    JHEEL = JACOBIAN_HEEL(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Nov-2021 15:09:51

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
t14 = -t12;
t15 = -t13;
Jheel = reshape([1.0,0.0,0.0,1.0,t10+t14+l_DE.*t2+l_OB.*t2,t11+t15+l_DE.*t3+l_OB.*t3,t10+t14,t11+t15,t14,t15,0.0,0.0],[2,6]);
