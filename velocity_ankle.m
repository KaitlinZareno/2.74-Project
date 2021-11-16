function drE = velocity_ankle(in1,in2)
%VELOCITY_ANKLE
%    DRE = VELOCITY_ANKLE(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Nov-2021 19:10:08

dth1 = in1(8,:);
dth2 = in1(9,:);
dx = in1(6,:);
dy = in1(7,:);
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
drE = [dx+dth1.*(l_AC.*t5+l_DE.*t2+l_OB.*t2)+dth2.*l_AC.*t5;dy+dth1.*(l_AC.*t6+l_DE.*t3+l_OB.*t3)+dth2.*l_AC.*t6;0.0];
