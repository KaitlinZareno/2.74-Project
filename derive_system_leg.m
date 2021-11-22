name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms x y t th1 th2 th3 ths xs ys dx dy dth1 dth2 dth3 dths dxs dys ddx ddy ddth1 ddth2 ddth3 ddths ddxs ddys real
syms m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_C_s ls g k real
syms l_OA l_OB l_AC l_DE l_IG l_GH l_heela real 
syms tau1 tau2 tau3 taus Fx Fy real
syms Ir N real

% Group them
q   = [x; y; th1 ; th2 ; th3; ths];      % generalized coordinates
dq  = [dx; dy; dth1 ; dth2; dth3 ; dths];    % first time derivatives
ddq = [ddx; ddy; ddth1;ddth2; ddth3 ; ddths];  % second time derivatives
u   = [tau1 ; tau2; tau3; taus];     % controls
F   = [Fx ; Fy];

p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB l_AC l_DE l_IG l_GH l_heela ls l_C_s g k]';        % parameters

% Generate Vectors and Derivatives

% [DEFINE HELPFUL UNIT VECTORS HERE] %
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

th1hat = sin(th1)*ihat - cos(th1)*jhat;
th2hat = sin(th2 + th1)*ihat - cos(th2 + th1)*jhat;
th3hat = sin(th3 + th2 + th1)*ihat - cos(th3+ th2 + th1)*jhat;

thshat = sin(ths)*ihat - cos(ths)*jhat;; %WANT SWING LEG TO BE THE MIRROR IMAGE OF TH1
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

%SWING LEG
rMs = r0 + l_C_s*thshat;
rHs = r0 + ls*thshat;

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

%SWING LEG
%velocity
drHs = ddt(rHs);
vS = ddt(rMs);


% ---------------------------------- %

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

% [DEFINE LINEAR AND ANGULAR VELOCITIES HERE] %

T0 = (1/2)*m0*dot(vR0', vR0); %INCLUDE ANOTHER ROTATIONAL INTERTIA THING FOR BOOM?
T1 = (1/2)*m1*dot(vM1', vM1)+ (1/2)*(I1)*dth1^2; 
T2 = (1/2)*m2*dot(vM2', vM2)+ (1/2)*(I2)*dth2^2;
T3 = (1/2)*m3*dot(vM3', vM3)+ (1/2)*(I3)*dth2^2;
T4 = (1/2)*m4*dot(vM4', vM4)+ (1/2)*(I4)*dth1^2;

T5 = (1/2)*m5*dot(vM5', vM5)+ (1/2)*(I5)*dth3^2;


%SWING LEG
Ts = (1/2)*ms*dot(vS', vS)+ (1/2)*(Is)*dths^2;


% ---------------------------------- %

% Compute kinetic energy about rotational movement
% Note: Don't forget computing kinetic energy for rotor inertias 

% Compute potential energy

P1 = m1*g*dot(rM1,jhat);
P2 = m2*g*dot(rM2,jhat);
P3 = m3*g*dot(rM3,jhat);
P4 = m4*g*dot(rM4,jhat);

%FOOT
P5 = m5*g*dot(rM5,jhat);
Pky = 1/2*k*(dot(rI,jhat)^2-dot(r_heela,jhat)^2); %break down spring force into x and y components
Pkx = 1/2*k*(dot(rI,ihat)^2-dot(r_heela,ihat)^2); 

%SWING LEG
Ps = ms*g*dot(rMs,jhat);

% Compute entire system energy
%T = simplify(T1 + T2 + T3 + T4 +T5); 
T = simplify(T0 + T1 + T2 + T3 + T4 + T5 + Ts); %with x movement
%V = P1+P2+P3+P4+P5;
V = P1+P2+P3+P4+P5+Pky+Pkx + Ps;

% Find Generalized forces
Q_tau1 = M2Q(tau1*khat,dth1*khat);
Q_tau2 = M2Q(tau2*khat,dth2*khat);
%foot
Q_tau3 = M2Q(tau3*khat,dth3*khat);
%swing
Q_taus = M2Q(taus*khat,dths*khat);

Q = Q_tau1 + Q_tau2 + Q_tau3 + Q_taus;

% Assemble the array of cartesian coordinates of the key points
keypoints = [r0(1:2) rA(1:2) rB(1:2) rC(1:2) rD(1:2) rE(1:2) rI(1:2) rH(1:2) r_heela(1:2) rHs(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(g,ddq);
b = A*ddq - g;

% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});

%toe 
matlabFunction(rH,'file',['position_toe'],'vars',{z p}); 
matlabFunction(drH,'file',['velocity_toe'],'vars',{z p});
Jtoe = jacobian(rH,q);
matlabFunction(Jtoe,'file',['jacobian_toe'],'vars',{z p});
%heel 
matlabFunction(rI,'file',['position_heel'],'vars',{z p});
matlabFunction(drI,'file',['velocity_heel'],'vars',{z p});
Jheel= jacobian(rI,q);
matlabFunction(Jheel,'file',['jacobian_heel'],'vars',{z p});

%heel attachment 
matlabFunction(r_heela,'file',['position_heel_attachment'],'vars',{z p});

%ankle 
matlabFunction(rE,'file',['position_ankle'],'vars',{z p});
matlabFunction(drE,'file',['velocity_ankle'],'vars',{z p});
Jankle= jacobian(rE,q);
matlabFunction(Jankle,'file',['jacobian_ankle'],'vars',{z p});

%Swing leg toe
matlabFunction(rHs,'file',['position_toe_swing'],'vars',{z p}); 
matlabFunction(drHs,'file',['velocity_toe_swing'],'vars',{z p});
Jtoes = jacobian(rHs,q);
matlabFunction(Jtoes,'file',['jacobian_toe_swing'],'vars',{z p});

matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
