function r2_heela = position_heel_attachment2(in1,in2)
%POSITION_HEEL_ATTACHMENT2
%    R2_HEELA = POSITION_HEEL_ATTACHMENT2(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    09-Dec-2021 10:40:32

l_AC = in2(23,:);
l_OA = in2(21,:);
l_heela = in2(27,:);
th1 = in1(3,:);
th2 = in1(4,:);
x = in1(1,:);
y = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
r2_heela = [x+l_OA.*t3+l_heela.*t3+l_AC.*sin(t4);y-l_OA.*t2-l_heela.*t2-l_AC.*cos(t4);0.0];
