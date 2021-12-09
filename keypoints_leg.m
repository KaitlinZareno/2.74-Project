function keypoints = keypoints_leg(in1,in2)
%KEYPOINTS_LEG
%    KEYPOINTS = KEYPOINTS_LEG(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    09-Dec-2021 10:40:33

l_AC = in2(23,:);
l_DE = in2(24,:);
l_GH = in2(26,:);
l_IG = in2(25,:);
l_OA = in2(21,:);
l_OB = in2(22,:);
l_heela = in2(27,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
x = in1(1,:);
y = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = th1+th3;
t6 = l_DE.*t2;
t7 = l_OA.*t2;
t8 = l_OB.*t2;
t9 = l_heela.*t2;
t10 = cos(t4);
t11 = cos(t5);
t12 = l_DE.*t3;
t13 = l_OA.*t3;
t14 = l_OB.*t3;
t15 = l_heela.*t3;
t16 = sin(t4);
t17 = sin(t5);
t18 = l_AC.*t10;
t19 = l_GH.*t11;
t20 = l_IG.*t11;
t21 = l_AC.*t16;
t22 = l_GH.*t17;
t23 = l_IG.*t17;
t24 = t13+x;
t25 = t14+x;
t26 = -t6;
t27 = -t7;
t28 = -t8;
t29 = -t9;
t30 = -t18;
t31 = -t19;
t32 = t27+y;
t33 = t28+y;
t34 = -t23;
t35 = t21+t24;
t36 = t21+t25;
t37 = t12+t36;
t38 = t15+t35;
t39 = t30+t32;
t40 = t30+t33;
t41 = t22+t37;
t42 = t26+t40;
t43 = t29+t39;
t44 = t34+t37;
t45 = t20+t42;
t46 = t31+t42;
keypoints = reshape([x,y,t24,t32,t25,t33,t35,t39,t36,t40,t37,t42,t44,t45,t41,t46,t38,t43,t24,t32,t25,t33,t35,t39,t36,t40,t37,t42,t44,t45,t41,t46,t38,t43],[2,17]);
