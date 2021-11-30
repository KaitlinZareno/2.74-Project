function rI = position_heel(in1,in2)
%POSITION_HEEL
%    RI = POSITION_HEEL(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Nov-2021 15:09:51

l_AC = in2(23,:);
l_DE = in2(24,:);
l_IG = in2(25,:);
l_OB = in2(22,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
x = in1(1,:);
y = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = t4+th3;
rI = [x+l_DE.*t3+l_OB.*t3+l_AC.*sin(t4)-l_IG.*sin(t5);y-l_DE.*t2-l_OB.*t2-l_AC.*cos(t4)+l_IG.*cos(t5);0.0];
