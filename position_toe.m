function rH = position_toe(in1,in2)
%POSITION_TOE
%    RH = POSITION_TOE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Nov-2021 19:10:06

l_AC = in2(21,:);
l_DE = in2(22,:);
l_GH = in2(24,:);
l_OB = in2(20,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
x = in1(1,:);
y = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = t4+th3;
rH = [x+l_DE.*t3+l_OB.*t3+l_AC.*sin(t4)+l_GH.*sin(t5);y-l_DE.*t2-l_OB.*t2-l_AC.*cos(t4)-l_GH.*cos(t5);0.0];
