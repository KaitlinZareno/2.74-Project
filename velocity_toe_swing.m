function drHs = velocity_toe_swing(in1,in2)
%VELOCITY_TOE_SWING
%    DRHS = VELOCITY_TOE_SWING(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    07-Dec-2021 16:33:56

dth1 = in1(4,:);
dth2 = in1(5,:);
dths = in1(6,:);
l_AC = in2(26,:);
l_DE = in2(27,:);
l_OB = in2(25,:);
ls = in2(23,:);
th1 = in1(1,:);
th2 = in1(2,:);
ths = in1(3,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = cos(t4);
t6 = sin(t4);
drHs = [-dth1.*(-l_AC.*t6+l_DE.*t3+l_OB.*t3)-dths.*ls.*sin(ths)+dth2.*l_AC.*t6;-dth1.*(-l_AC.*t5+l_DE.*t2+l_OB.*t2)-dths.*ls.*cos(ths)+dth2.*l_AC.*t5;0.0];
