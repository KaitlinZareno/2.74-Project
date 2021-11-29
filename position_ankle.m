function rE = position_ankle(in1,in2)
%POSITION_ANKLE
%    RE = POSITION_ANKLE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    29-Nov-2021 14:34:50

l_AC = in2(23,:);
l_DE = in2(24,:);
l_OB = in2(22,:);
th1 = in1(3,:);
th2 = in1(4,:);
x = in1(1,:);
y = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
rE = [x+l_DE.*t3+l_OB.*t3+l_AC.*sin(t4);y-l_DE.*t2-l_OB.*t2-l_AC.*cos(t4);0.0];
