name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms x t th1 th2 th3 dx dth1 dth2 dth3 ddx ddth1 ddth2 ddth3 real
syms m0 m1 m2 m3 m4 m5 I1 I2 I3 I4 I5 l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 g k real
syms l_OA l_OB l_AC l_DE l_IG l_GH l_heela real 
syms tau1 tau2 tau3 Fx Fy real
syms Ir N real

% Group them
q   = [x; th1 ; th2 ; th3];      % generalized coordinates
dq  = [dx; dth1 ; dth2; dth3];    % first time derivatives
ddq = [ddx; ddth1;ddth2; ddth3];  % second time derivatives
u   = [tau1 ; tau2; tau3];     % controls
F   = [Fx ; Fy];

p   = [m0 m1 m2 m3 m4 m5 I1 I2 I3 I4 I5 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB l_AC l_DE l_IG l_GH l_heela g k]';        % parameters

% Generate Vectors and Derivatives

% [DEFINE HELPFUL UNIT VECTORS HERE] %
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

th1hat = sin(th1)*ihat - cos(th1)*jhat;
th2hat = sin(th2 + th1)*ihat - cos(th2 + th1)*jhat;
th3hat = sin(th3 + th2 + th1)*ihat - cos(th3+ th2 + th1)*jhat;
% ---------------------------------- %

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

% [DEFINE KINEMATIC VECTORS HERE] %
r0 = x*ihat;  %ADD X VALUE TO INDICATE MOVEMTN

rA = r0 + l_OA*th1hat;  %ADD X VALUE TO INDICATE MOVEMTN
rM1 = r0 + l_O_m1*th1hat;
rB = r0 + l_OB*th1hat;

rM3 = rA+ l_A_m3 * th2hat;
rM2 = rB + l_B_m2 * th2hat;

rC = rA + l_AC * th2hat;
rD = rB + l_AC*th2hat;
rE = rD + l_DE*th1hat;

rM4 = rC + l_C_m4*th1hat;

%FOOT
r_heela = rC + l_heela*th1hat;
rH = rE + l_GH*th3hat; %two havles of the same foot
rI = rE-l_IG*th3hat; %second half MAKE NEGATIVE TO FLIP DIRECTION OF VECTOR, heel position
rM5 = rE + c5*th3hat;

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

% ---------------------------------- %

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

% [DEFINE LINEAR AND ANGULAR VELOCITIES HERE] %

T0 = (1/2)*m0*dot(vR0', vR0); %INCLUDE ANOTHER ROTATIONAL INTERTIA THING FOR BOOM?
T1 = (1/2)*m1*dot(vM1', vM1)+ (1/2)*(I1)*dth1^2; %TRANSPOSE rM1?? dot(rM1',rM1)
T2 = (1/2)*m2*dot(vM2', vM2)+ (1/2)*(I2)*dth2^2;
T3 = (1/2)*m3*dot(vM3', vM3)+ (1/2)*(I3)*dth2^2;
T4 = (1/2)*m4*dot(vM4', vM4)+ (1/2)*(I4)*dth1^2;

%FOOT
T5 = (1/2)*m5*dot(vM5', vM5)+ (1/2)*(I5)*dth3^2;

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

% Compute entire system energy
%T = simplify(T1 + T2 + T3 + T4 +T5); 
T = simplify(T0 + T1 + T2 + T3 + T4 +T5); %with x movement
%V = P1+P2+P3+P4+P5;
V = P1+P2+P3+P4+P5+Pky+Pkx;

% Find Generalized forces
Q_tau1 = M2Q(tau1*khat,dth1*khat);
Q_tau2 = M2Q(tau2*khat,dth2*khat);
%foot
Q_tau3 = M2Q(tau3*khat,dth3*khat);
Q = Q_tau1 + Q_tau2 + Q_tau3;

% Assemble the array of cartesian coordinates of the key points
keypoints = [r0(1:2) rA(1:2) rB(1:2) rC(1:2) rD(1:2) rE(1:2) rI(1:2) rH(1:2) r_heela(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(g,ddq);
b = A*ddq - g;

% Compute foot jacobian
Jleg = jacobian(rH,q);%IS IT RH __ TIP OF FOOT

% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(rH,'file',['position_foot'],'vars',{z p}); %FOOT CHANGES
matlabFunction(drH,'file',['velocity_foot'],'vars',{z p});
matlabFunction(Jleg,'file',['jacobian_foot'],'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
