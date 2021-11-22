function A = A_leg(in1,in2)
%A_LEG
%    A = A_LEG(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    21-Nov-2021 18:24:32

I1 = in2(8,:);
I2 = in2(9,:);
I3 = in2(10,:);
I4 = in2(11,:);
I5 = in2(12,:);
Is = in2(13,:);
c5 = in2(20,:);
l_AC = in2(23,:);
l_A_m3 = in2(18,:);
l_B_m2 = in2(17,:);
l_C_m4 = in2(19,:);
l_C_s = in2(29,:);
l_DE = in2(24,:);
l_OA = in2(21,:);
l_OB = in2(22,:);
l_O_m1 = in2(16,:);
m0 = in2(1,:);
m1 = in2(2,:);
m2 = in2(3,:);
m3 = in2(4,:);
m4 = in2(5,:);
m5 = in2(6,:);
ms = in2(7,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
ths = in1(6,:);
t2 = cos(th1);
t3 = cos(ths);
t4 = sin(th1);
t5 = sin(ths);
t6 = th1+th2;
t7 = c5.^2;
t8 = l_AC.^2;
t9 = l_A_m3.^2;
t10 = l_B_m2.^2;
t11 = l_C_s.^2;
t12 = l_O_m1.^2;
t61 = m0+m1+m2+m3+m4+m5+ms;
t13 = l_C_m4.*t2;
t14 = l_DE.*t2;
t15 = l_OA.*t2;
t16 = l_OB.*t2;
t17 = cos(t6);
t18 = l_C_m4.*t4;
t19 = l_DE.*t4;
t20 = l_OA.*t4;
t21 = l_OB.*t4;
t22 = sin(t6);
t23 = t6+th3;
t24 = l_O_m1.*m1.*t2;
t25 = l_C_s.*ms.*t3;
t27 = l_O_m1.*m1.*t4;
t28 = l_C_s.*ms.*t5;
t26 = cos(t23);
t29 = sin(t23);
t30 = t13.*2.0;
t31 = t14.*2.0;
t32 = t15.*2.0;
t33 = t16.*2.0;
t34 = t18.*2.0;
t35 = t19.*2.0;
t36 = t20.*2.0;
t37 = t21.*2.0;
t38 = t17.^2;
t39 = t22.^2;
t40 = l_AC.*t17;
t41 = l_A_m3.*t17;
t42 = l_B_m2.*t17;
t43 = l_AC.*t22;
t44 = l_A_m3.*t22;
t45 = l_B_m2.*t22;
t46 = t40.*2.0;
t47 = t41.*2.0;
t48 = t42.*2.0;
t49 = t43.*2.0;
t50 = t44.*2.0;
t51 = t45.*2.0;
t52 = c5.*t26;
t53 = m4.*t40;
t54 = m3.*t41;
t55 = m2.*t42;
t56 = c5.*t29;
t57 = m4.*t43;
t58 = m3.*t44;
t59 = m2.*t45;
t65 = t15+t41;
t66 = t16+t42;
t67 = t20+t44;
t68 = t21+t45;
t75 = t13+t15+t40;
t76 = t18+t20+t43;
t60 = m5.*t52;
t62 = m5.*t56;
t63 = t52.*2.0;
t64 = t56.*2.0;
t69 = t32+t47;
t70 = t33+t48;
t71 = t36+t50;
t72 = t37+t51;
t73 = t40+t52;
t74 = t43+t56;
t79 = t30+t32+t46;
t80 = t50.*t67;
t81 = t51.*t68;
t84 = t34+t36+t49;
t87 = t47.*t65;
t88 = t48.*t66;
t93 = t46.*t75;
t96 = t49.*t76;
t77 = t46+t63;
t78 = t49+t64;
t82 = (m3.*t69)./2.0;
t83 = (m2.*t70)./2.0;
t85 = (m3.*t71)./2.0;
t86 = (m2.*t72)./2.0;
t91 = t14+t16+t73;
t92 = t19+t21+t74;
t94 = t63.*t73;
t95 = t64.*t74;
t97 = (m4.*t84)./2.0;
t98 = (m4.*t79)./2.0;
t105 = t80+t87;
t106 = t81+t88;
t114 = t93+t96;
t89 = (m5.*t77)./2.0;
t90 = (m5.*t78)./2.0;
t99 = t31+t33+t77;
t100 = t35+t37+t78;
t101 = t64.*t92;
t102 = t63.*t91;
t107 = t73.*t91.*2.0;
t108 = t74.*t92.*2.0;
t109 = (m3.*t105)./2.0;
t110 = (m2.*t106)./2.0;
t113 = t94+t95;
t116 = (m4.*t114)./2.0;
t103 = (m5.*t100)./2.0;
t104 = (m5.*t99)./2.0;
t111 = t53+t54+t55+t89;
t112 = t57+t58+t59+t90;
t115 = (m5.*t113)./2.0;
t117 = t101+t102;
t119 = t107+t108;
t118 = (m5.*t117)./2.0;
t120 = (m5.*t119)./2.0;
t121 = t27+t85+t86+t97+t103;
t122 = t24+t82+t83+t98+t104;
t123 = t109+t110+t116+t120;
A = reshape([t61,0.0,t122,t111,t60,t25,0.0,t61,t121,t112,t62,t28,t122,t121,I1+I4+(m1.*(t2.^2.*t12.*2.0+t4.^2.*t12.*2.0))./2.0+(m3.*(t65.^2.*2.0+t67.^2.*2.0))./2.0+(m2.*(t66.^2.*2.0+t68.^2.*2.0))./2.0+(m4.*(t75.^2.*2.0+t76.^2.*2.0))./2.0+(m5.*(t91.^2.*2.0+t92.^2.*2.0))./2.0,t123,t118,0.0,t111,t112,t123,I2+I3+(m5.*(t73.^2.*2.0+t74.^2.*2.0))./2.0+(m4.*(t8.*t38.*2.0+t8.*t39.*2.0))./2.0+(m3.*(t9.*t38.*2.0+t9.*t39.*2.0))./2.0+(m2.*(t10.*t38.*2.0+t10.*t39.*2.0))./2.0,t115,0.0,t60,t62,t118,t115,I5+(m5.*(t7.*t26.^2.*2.0+t7.*t29.^2.*2.0))./2.0,0.0,t25,t28,0.0,0.0,0.0,Is+(ms.*(t3.^2.*t11.*2.0+t5.^2.*t11.*2.0))./2.0],[6,6]);
