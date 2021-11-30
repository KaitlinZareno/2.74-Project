function r0 = position_hip(in1,in2)
%POSITION_HIP
%    R0 = POSITION_HIP(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Nov-2021 15:22:39

l_AC = in2(26,:);
l_DE = in2(27,:);
l_OB = in2(25,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
r0 = [l_DE.*t2+l_OB.*t2-l_AC.*cos(t4);-l_DE.*t3-l_OB.*t3+l_AC.*sin(t4);0.0];
