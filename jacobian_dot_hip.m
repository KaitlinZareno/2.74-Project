function dJhip = jacobian_dot_hip(in1,in2)
%JACOBIAN_DOT_HIP
%    DJHIP = JACOBIAN_DOT_HIP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    29-Nov-2021 15:42:50

dth1 = in1(4,:);
dth2 = in1(5,:);
l_AC = in2(26,:);
l_DE = in2(27,:);
l_OB = in2(25,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = cos(t4);
t6 = sin(t4);
t7 = dth2.*l_AC.*t5;
t8 = dth2.*l_AC.*t6;
t9 = -t8;
dJhip = reshape([t7-dth1.*(-l_AC.*t5+l_DE.*t2+l_OB.*t2),t9+dth1.*(-l_AC.*t6+l_DE.*t3+l_OB.*t3),0.0,t7+dth1.*l_AC.*t5,t9-dth1.*l_AC.*t6,0.0,0.0,0.0,0.0],[3,3]);
