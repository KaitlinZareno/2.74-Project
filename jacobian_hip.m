function Jhip = jacobian_hip(in1,in2)
%JACOBIAN_HIP
%    JHIP = JACOBIAN_HIP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    29-Nov-2021 14:25:00

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
t7 = l_AC.*t5;
t8 = l_AC.*t6;
Jhip = reshape([t8-l_DE.*t3-l_OB.*t3,t7-l_DE.*t2-l_OB.*t2,0.0,t8,t7,0.0,0.0,0.0,0.0],[3,3]);
