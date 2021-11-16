function Jankle = jacobian_ankle(in1,in2)
%JACOBIAN_ANKLE
%    JANKLE = JACOBIAN_ANKLE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Nov-2021 19:10:08

l_AC = in2(21,:);
l_DE = in2(22,:);
l_OB = in2(20,:);
th1 = in1(3,:);
th2 = in1(4,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = cos(t4);
t6 = sin(t4);
t7 = l_AC.*t5;
t8 = l_AC.*t6;
Jankle = reshape([1.0,0.0,0.0,0.0,1.0,0.0,t7+l_DE.*t2+l_OB.*t2,t8+l_DE.*t3+l_OB.*t3,0.0,t7,t8,0.0,0.0,0.0,0.0],[3,5]);
