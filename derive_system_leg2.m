name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms x y t th1 th2 th3 th21 th22 th23 dx dy dth1 dth2 dth3 dth21 dth22 dth23 ddx ddy ddth1 ddth2 ddth3 ddth21 ddth22 ddth23 real

syms m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_C_s ls l_anklerest g k real
syms l_OA l_OB l_AC l_DE l_IG l_GH l_heela real

syms tau1 tau2 tau3 tau21 tau22 tau23 Fx Fy real
syms Ir N real

% Group them
q   = [x y th1 th2 th3 th21 th22 th23]';      % generalized coordinates
dq  = [dx dy dth1 dth2 dth3 dth21 dth22 dth23]';    % first time derivatives
ddq = [ddx ddy ddth1 ddth2 ddth3 ddth21 ddth22 ddth23]';  % second time derivatives
u   = [tau1 tau2 tau3 tau21 tau22 tau23]';     % controls
F   = [Fx ; Fy];

p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB ...
       l_AC l_DE l_IG l_GH l_heela ls l_C_s g k l_anklerest]';        % parameters

% Generate Vectors and Derivatives

% [DEFINE HELPFUL UNIT VECTORS HERE] %
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

%LEG 1
th1hat = sin(th1)*ihat - cos(th1)*jhat;
th2hat = sin(th2 + th1)*ihat - cos(th2 + th1)*jhat;
th3hat = sin(th3 + th2 + th1)*ihat - cos(th3 + th2+ th1)*jhat; %th2 hat not relevant 

%LEG 2
th1hat_s = sin(th21)*ihat - cos(th21)*jhat;
th2hat_s = sin(th22 + th21)*ihat - cos(th22 + th21)*jhat;
th3hat_s = sin(th23 + th22 + th21)*ihat - cos(th23 + th22+ th21)*jhat;
% ---------------------------------- %

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

% [DEFINE KINEMATIC VECTORS HERE] %
r0 = x*ihat + y*jhat;  %ADD X VALUE TO INDICATE MOVEMTN

rA = r0 + l_OA*th1hat;  %ADD X VALUE TO INDICATE MOVEMTN
rM1 = r0 + l_O_m1*th1hat;
rB = r0 + l_OB*th1hat;

rM3 = rA+ l_A_m3 * th2hat;
rM2 = rB + l_B_m2 * th2hat;

rC = rA + l_AC * th2hat;
rD = rB + l_AC*th2hat; 
rE = rD + l_DE*th1hat; %defined rc+lce*th1hat?

rM4 = rC + l_C_m4*th1hat;

%FOOT
r_heela = rC + l_heela*th1hat;
rH = rE + l_GH*th3hat; %two havles of the same foot
rI = rE-l_IG*th3hat; %second half MAKE NEGATIVE TO FLIP DIRECTION OF VECTOR, heel position
rM5 = rE + c5*th3hat;

%SECOND LEG -- same variables denoted with a r2x
r2A = r0 + l_OA*th1hat_s;  %ADD X VALUE TO INDICATE MOVEMTN
r2M1 = r0 + l_O_m1*th1hat_s;
r2B = r0 + l_OB*th1hat_s;

r2M3 = r2A+ l_A_m3 * th2hat_s;
r2M2 = r2B + l_B_m2 * th2hat_s;

r2C = r2A + l_AC * th2hat_s;
r2D = r2B + l_AC*th2hat_s; 
r2E = r2D + l_DE*th1hat_s;

r2M4 = r2C + l_C_m4*th1hat_s;

%FOOT
r2_heela = r2C + l_heela*th1hat_s;
r2H = r2E + l_GH*th3hat_s; %two havles of the same foot
r2I = r2E-l_IG*th3hat_s; %second half MAKE NEGATIVE TO FLIP DIRECTION OF VECTOR, heel position
r2M5 = r2E + c5*th3hat_s;

% -------------------------------------------------------------------------%        
 
%CALCULATE DERIVATIVES
%velocity
vR0 = ddt(r0);
vM1 = ddt(rM1);
vM2 = ddt(rM2);
vM3 = ddt(rM3);
vM4 = ddt(rM4);
drE = ddt(rE);

%foot
vM5 = ddt(rM5);
drH = ddt(rH);
drI = ddt(rI);

%SECOND LEG
v2M1 = ddt(r2M1);
v2M2 = ddt(r2M2);
v2M3 = ddt(r2M3);
v2M4 = ddt(r2M4);
dr2E = ddt(r2E);

%foot
v2M5 = ddt(r2M5);
dr2H = ddt(r2H);
dr2I = ddt(r2I);


% ---------------------------------- %

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

%KINETIC ENERGY
omega1 = dth1;
omega2 = dth1 + dth2;
omega3 = dth1 + dth2;
omega4 = dth1;
omega5 = dth1 + dth3; 

T0 = (1/2)*m0*dot(vR0', vR0); %INCLUDE ANOTHER ROTATIONAL INTERTIA THING FOR BOOM?
T1 = (1/2)*m1*dot(vM1', vM1)+ (1/2)*(I1)*omega1^2; 
T2 = (1/2)*m2*dot(vM2', vM2)+ (1/2)*(I2)*omega2^2;
T3 = (1/2)*m3*dot(vM3', vM3)+ (1/2)*(I3)*omega3^2;
T4 = (1/2)*m4*dot(vM4', vM4)+ (1/2)*(I4)*omega4^2;

T5 = (1/2)*m5*dot(vM5', vM5)+ (1/2)*(I5)*omega5^2;

%SECOND LEG
omega21 = dth21;
omega22 = dth21 + dth22;
omega23 = dth21 + dth22;
omega24 = dth21;
omega25 = dth21 + dth23; 

T21 = (1/2)*m1*dot(v2M1', v2M1)+ (1/2)*(I1)*omega21^2; % will have same masses and inertias as leg 1
T22 = (1/2)*m2*dot(v2M2', v2M2)+ (1/2)*(I2)*omega22^2;
T23 = (1/2)*m3*dot(v2M3', v2M3)+ (1/2)*(I3)*omega23^2;
T24 = (1/2)*m4*dot(v2M4', v2M4)+ (1/2)*(I4)*omega24^2;

T25 = (1/2)*m5*dot(v2M5', v2M5)+ (1/2)*(I5)*omega25^2;


%ROTOR INERTIA ENERGY                CHECK
T1r = (1/2)*Ir*(N*dth1)^2;           %hip 
T2r = (1/2)*Ir*(dth1 + N*dth2)^2;    %knee
T3r = (1/2)*Ir*(dth1 + N*dth3)^2;    %foot

%ROTOR INERTIA FOOT 2
T21r = (1/2)*Ir*(N*dth21)^2;
T22r = (1/2)*Ir*(dth21 + N*dth22)^2;
T23r = (1/2)*Ir*(dth21 + N*dth23)^2; 

%Compute system kinetic
T1 = T0 + T1 + T2 + T3 + T4 + T5 + T1r + T2r + T3r;
T2 = T21 + T22 + T23 + T24 + T25 + T21r + T22r + T23r; %x,y movement included in T1
T = simplify(T1+T2);

% ---------------------------------- %

%POTENTIAL ENERGY
P1 = m1*g*dot(rM1,jhat);
P2 = m2*g*dot(rM2,jhat);
P3 = m3*g*dot(rM3,jhat);
P4 = m4*g*dot(rM4,jhat);

%FOOT
P5 = m5*g*dot(rM5,jhat);
% Pky = 1/2*k*(dot(rI,-ihat)^2-dot(r_heela,-ihat)^2); %break down spring force into x and y components
% Pkx = 1/2*k*(dot(rI,jhat)^2-dot(r_heela,jhat)^2); 
dist = r_heela-rI;
Pa = 1/2*k*((dist-l_anklerest)'*(dist-l_anklerest)); %break down spring force into x and y components


%LEG 2
P21 = m1*g*dot(r2M1,jhat);
P22 = m2*g*dot(r2M2,jhat);
P23 = m3*g*dot(r2M3,jhat);
P24 = m4*g*dot(r2M4,jhat);

%FOOT
P25 = m5*g*dot(r2M5,jhat);
dist2 = r2_heela-r2I;
P2a = 1/2*k*((dist2-l_anklerest)'*(dist2-l_anklerest)); %break down spring force into x and y components

% Compute entire system potential energy 
V1 = P1+P2+P3+P4+P5+Pa;
V2 = P21+P22+P23+P24+P25+P2a;
V = V1+V2;

% ---------------------------------- %

% Find Generalized forces
Q_tau1 = M2Q(tau1*khat,omega1*khat);
Q_tau2 = M2Q(tau2*khat,omega2*khat);
Q_tau2R= M2Q(-tau2*khat,omega1*khat);

%foot
Q_tau3 = M2Q(tau3*khat,omega3*khat);
Q_tau3R= M2Q(-tau3*khat,omega1*khat); %WHAT IS THIS

%LEG 2
Q_tau21 = M2Q(tau21*khat,omega21*khat);
Q_tau22 = M2Q(tau22*khat,omega22*khat);
Q_tau22R= M2Q(-tau22*khat,omega21*khat);

%foot
Q_tau23 = M2Q(tau23*khat,omega23*khat);
Q_tau23R= M2Q(-tau23*khat,omega21*khat); %WHAT IS THIS

Q1 = Q_tau1 + Q_tau2 + Q_tau3 + Q_tau2R + Q_tau3R;
Q2 = Q_tau21 + Q_tau22 + Q_tau23 + Q_tau22R + Q_tau23R;
Q = Q1+Q2;

% Assemble the array of cartesian coordinates of the key points
keypoints = [r0(1:2) rA(1:2) rB(1:2) rC(1:2) rD(1:2) rE(1:2) rI(1:2) rH(1:2) r_heela(1:2) ...
             r2A(1:2) r2B(1:2) r2C(1:2) r2D(1:2) r2E(1:2) r2I(1:2) r2H(1:2) r2_heela(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(eom,ddq);
b = A*ddq - eom;

% Equations of motion are
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(V, q)');
Corr_Joint_Sp = simplify( eom + Q - Grav_Joint_Sp - A*ddq);


% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(Grav_Joint_Sp ,'file', ['Grav_leg'] ,'vars',{z p});
matlabFunction(Corr_Joint_Sp ,'file', ['Corr_leg']     ,'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});

%FOOT 1
%toe 
Jtoe = jacobian(rH,q);
dJtoe= reshape( ddt(Jtoe(:)) , size(Jtoe) );

Jtoe  = Jtoe(1:2,:);
dJtoe = dJtoe(1:2,:);

matlabFunction(rH,'file',['position_toe'],'vars',{z p}); 
matlabFunction(drH,'file',['velocity_toe'],'vars',{z p});
matlabFunction(Jtoe,'file',['jacobian_toe'],'vars',{z p});
matlabFunction(dJtoe ,'file',['jacobian_dot_toe'],'vars',{z p});

%heel 
Jheel = jacobian(rI,q);
dJheel= reshape( ddt(Jheel(:)) , size(Jheel) );

Jheel  = Jheel(1:2,:);
dJheel = dJheel(1:2,:);

matlabFunction(rI,'file',['position_heel'],'vars',{z p});
matlabFunction(drI,'file',['velocity_heel'],'vars',{z p});
matlabFunction(Jheel,'file',['jacobian_heel'],'vars',{z p});
matlabFunction(dJheel ,'file',['jacobian_dot_heel'],'vars',{z p});

%heel attachment 
matlabFunction(r_heela,'file',['position_heel_attachment'],'vars',{z p});

%ankle 
Jankle = jacobian(rE,q);
dJankle= reshape( ddt(Jankle(:)) , size(Jankle) );

Jankle  = Jankle(1:2,:);
dJankle = dJankle(1:2,:);

matlabFunction(rE,'file',['position_ankle'],'vars',{z p});
matlabFunction(drE,'file',['velocity_ankle'],'vars',{z p});
matlabFunction(Jankle,'file',['jacobian_ankle'],'vars',{z p});
matlabFunction(dJankle,'file',['jacobian_dot_ankle'],'vars',{z p});

%FOOT 2
%toe 
Jtoe2 = jacobian(r2H,q);
dJtoe2= reshape( ddt(Jtoe2(:)) , size(Jtoe2) );

Jtoe2  = Jtoe2(1:2,:);
dJtoe2 = dJtoe2(1:2,:);

matlabFunction(r2H,'file',['position_toe2'],'vars',{z p}); 
matlabFunction(dr2H,'file',['velocity_toe2'],'vars',{z p});
matlabFunction(Jtoe2,'file',['jacobian_toe2'],'vars',{z p});
matlabFunction(dJtoe2 ,'file',['jacobian_dot_toe2'],'vars',{z p});

%heel 
Jheel2 = jacobian(r2I,q);
dJheel2= reshape( ddt(Jheel2(:)) , size(Jheel2) );

Jheel2  = Jheel2(1:2,:);
dJheel2 = dJheel2(1:2,:);

matlabFunction(r2I,'file',['position_heel2'],'vars',{z p});
matlabFunction(dr2I,'file',['velocity_heel2'],'vars',{z p});
matlabFunction(Jheel2,'file',['jacobian_heel2'],'vars',{z p});
matlabFunction(dJheel2 ,'file',['jacobian_dot_heel2'],'vars',{z p});

%heel attachment 
matlabFunction(r2_heela,'file',['position_heel_attachment2'],'vars',{z p});

%ankle 
Jankle2 = jacobian(r2E,q);
dJankle2= reshape( ddt(Jankle2(:)) , size(Jankle2) );

Jankle2  = Jankle2(1:2,:);
dJankle2 = dJankle2(1:2,:);

matlabFunction(r2E,'file',['position_ankle2'],'vars',{z p});
matlabFunction(dr2E,'file',['velocity_ankle2'],'vars',{z p});
matlabFunction(Jankle2,'file',['jacobian_ankle2'],'vars',{z p});
matlabFunction(dJankle2,'file',['jacobian_dot_ankle2f'],'vars',{z p});

matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
