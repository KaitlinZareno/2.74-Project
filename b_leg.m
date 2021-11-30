function b = b_leg(in1,in2,in3)
%B_LEG
%    B = B_LEG(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Nov-2021 15:22:36

dth1 = in1(4,:);
dth2 = in1(5,:);
dths = in1(6,:);
g = in3(16,:);
l_AC = in3(26,:);
l_CE = in3(28,:);
l_DE = in3(27,:);
l_OB = in3(25,:);
l_m1 = in3(18,:);
l_m2 = in3(19,:);
l_m3 = in3(20,:);
l_m4 = in3(21,:);
l_ms = in3(22,:);
m0 = in3(1,:);
m1 = in3(2,:);
m2 = in3(3,:);
m3 = in3(4,:);
m4 = in3(5,:);
ms = in3(7,:);
tau1 = in2(1,:);
tau2 = in2(2,:);
taus = in2(3,:);
th1 = in1(1,:);
th2 = in1(2,:);
ths = in1(3,:);
t2 = cos(th1);
t3 = cos(ths);
t4 = sin(th1);
t5 = sin(ths);
t6 = th1+th2;
t7 = l_CE.*t2;
t8 = l_DE.*t2;
t9 = l_OB.*t2;
t10 = l_m4.*t2;
t11 = cos(t6);
t12 = l_CE.*t4;
t13 = l_DE.*t4;
t14 = l_OB.*t4;
t15 = l_m4.*t4;
t16 = sin(t6);
t17 = dths.*l_ms.*t3;
t18 = dths.*l_ms.*t5;
t19 = l_AC.*t11;
t20 = l_m2.*t11;
t21 = l_m3.*t11;
t22 = l_AC.*t16;
t23 = l_m2.*t16;
t24 = l_m3.*t16;
t25 = -t7;
t26 = -t8;
t27 = -t12;
t28 = -t13;
t29 = dth1.*t19;
t30 = dth2.*t19;
t31 = dth1.*t20;
t32 = dth1.*t21;
t33 = dth2.*t20;
t34 = dth2.*t21;
t35 = dth1.*t22;
t36 = dth2.*t22;
t37 = dth1.*t23;
t38 = dth1.*t24;
t39 = dth2.*t23;
t40 = dth2.*t24;
t41 = -t19;
t42 = -t22;
t45 = t20+t25;
t46 = t21+t26;
t47 = t23+t27;
t48 = t24+t28;
t55 = -dth1.*(t7-t20);
t56 = -dth1.*(t8-t21);
t57 = -dth1.*(t12-t23);
t58 = -dth1.*(t13-t24);
t43 = -t30;
t44 = -t36;
t49 = t29+t30;
t50 = t31+t33;
t51 = t32+t34;
t52 = t35+t36;
t53 = t37+t39;
t54 = t38+t40;
t59 = t8+t9+t41;
t60 = t8+t10+t41;
t61 = t13+t14+t42;
t62 = t13+t15+t42;
t69 = t33+t55;
t70 = t34+t56;
t71 = t39+t57;
t72 = t40+t58;
t63 = dth1.*t61;
t64 = dth1.*t62;
t65 = dth1.*t59;
t66 = dth1.*t60;
t67 = t22.*t49.*2.0;
t68 = t19.*t52.*2.0;
t78 = t23.*t69.*2.0;
t79 = t24.*t70.*2.0;
t80 = t20.*t71.*2.0;
t81 = t21.*t72.*2.0;
t84 = t52.*t59.*2.0;
t85 = t49.*t61.*2.0;
t73 = -t68;
t74 = t43+t65;
t75 = t43+t66;
t76 = t44+t63;
t77 = t44+t64;
t86 = -t85;
t87 = t22.*(t30-t65).*-2.0;
t88 = t22.*(t30-t66).*-2.0;
t89 = t19.*(t36-t63).*-2.0;
t90 = t19.*(t36-t64).*-2.0;
t91 = t19.*(t36-t63).*2.0;
t92 = t19.*(t36-t64).*2.0;
t82 = t17+t74;
t83 = t18+t76;
t93 = t22.*t82.*2.0;
t94 = t19.*t83.*2.0;
t95 = -t94;
b = [tau1-(ms.*(t83.*(t30-t65).*2.0-t82.*(t36-t63).*2.0))./2.0+dth2.*((m2.*(t78-t80-t53.*(t7-t20).*2.0+t50.*(t12-t23).*2.0))./2.0+(m3.*(t79-t81-t54.*(t8-t21).*2.0+t51.*(t13-t24).*2.0))./2.0+(m4.*(t90+t49.*t62.*2.0-t52.*t60.*2.0+t22.*(t30-t66).*2.0))./2.0-(m0.*(t84+t86+t87+t91))./2.0-(ms.*(t84+t86+t93+t95))./2.0)-(dth1.*ms.*(t59.*t83.*2.0-t61.*t82.*2.0-t61.*(t30-t65).*2.0+t59.*(t36-t63).*2.0))./2.0+(dths.*ms.*(t18.*t59.*2.0-t17.*t61.*2.0))./2.0-g.*m4.*t60-g.*ms.*t59-g.*m2.*(t7-t20)-g.*m3.*(t8-t21)-g.*l_m1.*m1.*t2;tau2-dth2.*(m2.*(t78-t80+t20.*t53.*2.0-t23.*t50.*2.0).*(-1.0./2.0)-(m3.*(t79-t81+t21.*t54.*2.0-t24.*t51.*2.0))./2.0+(m0.*(t67+t73+t87+t91))./2.0+(m4.*(t67+t73+t88+t92))./2.0+(ms.*(t67+t73+t93+t95))./2.0)-(m0.*(t52.*(t30-t65).*2.0-t49.*(t36-t63).*2.0))./2.0-(m4.*(t52.*(t30-t66).*2.0-t49.*(t36-t64).*2.0))./2.0+(m2.*(t50.*t71.*2.0-t53.*t69.*2.0))./2.0+(m3.*(t51.*t72.*2.0-t54.*t70.*2.0))./2.0-(ms.*(t49.*t83.*2.0-t52.*t82.*2.0))./2.0-(dths.*ms.*(t18.*t19.*2.0-t17.*t22.*2.0))./2.0+g.*m2.*t20+g.*m4.*t19+g.*m3.*t21+g.*ms.*t19+(dth1.*ms.*(t87+t91-t93+t94))./2.0;taus+(ms.*(t17.*t83.*2.0-t18.*t82.*2.0))./2.0+(dth1.*ms.*(l_ms.*t5.*(t30-t65).*2.0-l_ms.*t3.*(t36-t63).*2.0))./2.0+(dth2.*ms.*(l_ms.*t5.*t49.*2.0-l_ms.*t3.*t52.*2.0))./2.0-(dths.*ms.*(l_ms.*t3.*t83.*2.0-l_ms.*t5.*t82.*2.0))./2.0-g.*l_ms.*ms.*t3];
