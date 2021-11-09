name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms x t th0 ths th1 th2 th3 dx dth0 dths dth1 dth2 dth3 ddx ddth0 ddths ddth1 ddth2 ddth3 real
syms mc m0 ms m1 m2 m3 m4 m5  Ic I0 Is I1 I2 I3 I4 I5 c_00 c_AM c_B1 c_C2 c_D3 c_E4 c_G5  g real
syms l_0A l_AB l_BC l_BD l_CE l_DF l_EF l_EG l_IG l_GH l_heela real 
syms tau tau0 taus tau1 tau2 tau3 Fx Fy real
syms Ir N k real %added spring constant stuff here 

% Group them
q   = [x ; th0 ; ths ; th1  ; th2 ; th3];      % generalized coordinates
dq  = [dx ; dth0 ; dths ; dth1 ; dth2 ; dth3];    % first time derivatives
ddq = [ddx ; ddth0 ; ddths ; ddth1 ; ddth2 ; ddth3];  % second time derivatives
u   = [tau0 ; taus ; tau1 ; tau2 ; tau3];     % controls
F   = [Fx ; Fy];

p   = [mc m0 ms m1 m2 m3 m4 m5  Ic I0 Is I1 I2 I3 I4 I5 c_00 c_AM c_B1 c_C2 c_D3 c_E4 c_G5 l_0A l_AB l_BC l_BD l_CE l_DF l_EF l_EG l_IG l_GH l_heela g k]';        % parameters

% Generate Vectors and Derivatives

% [DEFINE HELPFUL UNIT VECTORS HERE] %
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

th0hat = -cos(th0)*ihat - sin(th0)*jhat;
thshat = (-cos(th0)+sin(ths))*ihat - (sin(th0)-cos(ths))*jhat;
th1hat = (-cos(th0)+sin(ths+th1))*ihat - (sin(th0)-cos(ths+th1))*jhat;
th2hat = (-cos(th0)+sin(ths+th1+th2))*ihat - (sin(th0)-cos(ths+th1+th2))*jhat;
% %foot angle
% th3hat = (-cos(th0)+sin(ths+th1+th2+th3))*ihat - (sin(th0)-cos(ths+th1+th2+th3))*jhat;
% ---------------------------------- %

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

% [DEFINE KINEMATIC VECTORS HERE] %
r_cart = x*ihat;
r_cm0 = c_00*th0hat;
rA = l_0A*th0hat;
r_cms = rA + c_AM*thshat;
r_AB = rA + l_AB*thshat;
r_cm1 = r_AB + c_B1*th1hat;
r_BC = r_AB + l_BC*th1hat;
r_BD = r_AB + l_BD*th1hat;
r_cm2 = r_BC + c_C2*th2hat;
r_CE = r_BC + l_CE*th2hat;
r_cm3 = r_BD + c_D3*th2hat;
r_DF = r_BD + l_DF*th2hat;

r_cm4 = r_CE + c_E4*th1hat;
r_EG = r_CE + l_EG*th1hat;

%foot stuff
%heel attachment length lies on the RB r_EG
% r_heela = r_EG + l_heela*th1hat;
% 
% r_GH = r_EG + l_GH*th3hat; %two havles of the same foot
% %This vaue is equivalent to rh, gives position of heel
% r_IG = r_EG-l_IG*th3hat; %second half MAKE NEGATIVE TO FLIP DIRECTION OF VECTOR 
% r_cm5 = r_EG + c_G5*th3hat;

%heel attachment and heel points
% r_heela = ; %lies on the shin r_EG
% r_heel = ; %opposite of r_GH -- lie on a straight line

v_cart = ddt(r_cart);
v_cm0 = ddt(r_cm0);
v_cms = ddt(r_cms);
v_cm1 = ddt(r_cm1);
v_cm2 = ddt(r_cm2);
v_cm3 = ddt(r_cm3);
v_cm4 = ddt(r_cm4);

dr_EG = ddt(r_EG);

% v_cm5 = ddt(r_cm5);
% dr_GH = ddt(r_GH); %tip of foot
% dr_IG = ddt(r_IG);%heel of foot

% ---------------------------------- %

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

% [DEFINE LINEAR AND ANGULAR VELOCITIES HERE] %

T0 = (1/2)*m0*dot(v_cm0', v_cm0)+ (1/2)*(I0)*dth0^2; 
Ts = (1/2)*ms*dot(v_cms', v_cms)+ (1/2)*(Is)*dths^2;
T1 = (1/2)*m1*dot(v_cm1', v_cm1)+ (1/2)*(I1)*dth1^2;
T2 = (1/2)*m2*dot(v_cm2', v_cm2)+ (1/2)*(I2)*dth2^2;
T3 = (1/2)*m1*dot(v_cm3', v_cm3)+ (1/2)*(I3)*dth2^2;
T4 = (1/2)*m2*dot(v_cm4', v_cm4)+ (1/2)*(I4)*dth1^2;

Tc = (1/2)*mc*dot(v_cart', v_cart);
%foot
% T5 = (1/2)*m1*dot(v_cm5', v_cm5)+ (1/2)*(I5)*dth3^2; 

% ---------------------------------- %

% Compute kinetic energy about rotational movement
% Note: Don't forget computing kinetic energy for rotor inertias 

% Compute potential energy
P0 = m0*g*dot(r_cm0,jhat);
Ps = ms*g*dot(r_cms,jhat);
P1 = m1*g*dot(r_cm1,jhat);
P2 = m2*g*dot(r_cm2,jhat);
P3 = m3*g*dot(r_cm3,jhat);
P4 = m4*g*dot(r_cm4,jhat);

%foot
% P5 = m5*g*dot(r_cm5,jhat);
% Pkx = 1/2*k*(dot(r_IG,jhat)^2-dot(r_heela,jhat)^2); %break down spring force into x and y components
% Pky = 1/2*k*(dot(r_IG,ihat)^2-dot(r_heela,ihat)^2); 

% Compute entire system energy
T = simplify(T0+ Ts + T1 + T2 + T3 + T4 + Tc); 
V = P0+Ps+P1+P2+P3+P4;

% Find Generalized forces
Q_tau0 = M2Q(tau0*khat,dth0*khat);
Q_taus = M2Q(taus*khat,dths*khat);
Q_tau1 = M2Q(tau1*khat,dth1*khat);
Q_tau2 = M2Q(tau2*khat,dth2*khat);
%foot torque
% Q_tau3 = M2Q(tau3*khat,dth3*khat);
Q = Q_tau0+ Q_taus + Q_tau1 + Q_tau2;

% Assemble the array of cartesian coordinates of the key points
keypoints = [r_cart(1:2) rA(1:2) r_AB(1:2) r_BC(1:2) r_BD(1:2) r_CE(1:2) r_DF(1:2) r_EG(1:2)]
% keypoints with foot
% keypoints = [r_cart(1:2) rA(1:2) r_AB(1:2) r_BC(1:2) r_BD(1:2) r_CE(1:2) r_DF(1:2) r_EG(1:2) r_GH(1:2) r_IG(1:2) r_heela(1:2)]; %endpoints, ankle attachment point (anywhere you need lines) 

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(g,ddq);
b = A*ddq - g;

% Compute foot jacobian
Jleg = jacobian(r_EG,q); %leg until ankle
% Jleg = jacobian(r_GH,q); %calculate tip of foot jacobian

% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(r_EG,'file',['position_foot'],'vars',{z p}); %position is r_EG
matlabFunction(dr_EG,'file',['velocity_foot'],'vars',{z p}); %velocity is dr_EG
matlabFunction(Jleg,'file',['jacobian_foot'],'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
