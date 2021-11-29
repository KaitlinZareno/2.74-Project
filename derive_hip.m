name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t th1 th2 ths dth1 dth2 dths ddth1 ddth2 ddths real
syms m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is l_m1 l_m2 l_m3 l_m4 l_ms ls g k real
syms l_OA l_OB l_AC l_DE l_CE l_IG l_GH l_heela real 
syms tau1 tau2 taus Fx Fy real
syms Ir N real

% Group them
q   = [th1 ; th2 ; ths];      % generalized coordinates
dq  = [dth1 ; dth2; dths];    % first time derivatives
ddq = [ddth1;ddth2; ddths];  % second time derivatives
u   = [tau1 ; tau2; taus];     % controls
F   = [Fx ; Fy];

p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N g k l_m1 l_m2 l_m3 l_m4 l_ms ls l_OA l_OB l_AC l_DE l_CE l_IG l_GH l_heela]';        % parameters

% Generate Vectors and Derivatives

% [DEFINE HELPFUL UNIT VECTORS HERE] %
ihat = [1; 0; 0];
jhat = [0; -1; 0];
khat = cross(ihat,jhat);

th1hat = cos(th1)*ihat + sin(th1)*jhat;
th2hat = cos(th1 + th2)*ihat + sin(th1 + th2)*jhat;
% th3hat = cos(pi-th1)*ihat + sin(pi-th1)*jhat;

thshat = cos(ths)*ihat + sin(ths)*jhat; %WANT SWING LEG TO BE THE MIRROR IMAGE OF TH1
% ---------------------------------- %

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

% [DEFINE KINEMATIC VECTORS HERE] %
%Define from ankle
%rE = 0 0

%foot
% rH = l_GH*ihat;
rI = 0+l_IG*-ihat;
r_heela = l_heela*th1hat;

%leg
rD = l_DE*th1hat;
rC = l_CE*th1hat;
r_m1 = l_m1 *th1hat;

rB = rD - l_AC*th2hat;
rA = rC - l_AC*th2hat;
r_m2 = rC - l_m2 *th2hat;
r_m3 = rD - l_m3 *th2hat;

r0 = rB + l_OB*th1hat;  
r_m4 = rB + l_m4 *th1hat;

%SWING LEG
rMs = r0 + l_ms*thshat;
rHs = r0 + ls*thshat;

%CALCULATE DERIVATIVES
%velocity
vR0 = ddt(r0);
vM1 = ddt(r_m1);
vM2 = ddt(r_m2);
vM3 = ddt(r_m3);
vM4 = ddt(r_m4);

%SWING LEG
%velocity
drHs = ddt(rHs);
vS = ddt(rMs);


% ---------------------------------- %

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

%KINETIC ENERGY
T0 = (1/2)*m0*dot(vR0', vR0); %INCLUDE ANOTHER ROTATIONAL INTERTIA THING FOR BOOM?
T1 = (1/2)*m1*dot(vM1', vM1)+ (1/2)*(I1)*dth1^2; 
T2 = (1/2)*m2*dot(vM2', vM2)+ (1/2)*(I2)*dth2^2;
T3 = (1/2)*m3*dot(vM3', vM3)+ (1/2)*(I3)*dth2^2;
T4 = (1/2)*m4*dot(vM4', vM4)+ (1/2)*(I4)*dth1^2;

%SWING LEG
Ts = (1/2)*ms*dot(vS', vS)+ (1/2)*(Is)*dths^2;

%ROTOR INERTIA ENERGY
T1r = (1/2)*Ir*(N*dth1)^2;
T2r = (1/2)*Ir*(dth1 + N*dth2)^2;
Tsr = (1/2)*Ir*(dth1 + dth2 + N*dths)^2;

% ---------------------------------- %

%POTENTIAL ENERGY
P1 = m1*g*dot(r_m1,jhat);
P2 = m2*g*dot(r_m2,jhat);
P3 = m3*g*dot(r_m3,jhat);
P4 = m4*g*dot(r_m4,jhat);

%FOOT
Pky = 1/2*k*(dot(rI,jhat)^2-dot(r_heela,jhat)^2); %break down spring force into x and y components
Pkx = 1/2*k*(dot(rI,ihat)^2-dot(r_heela,ihat)^2); 

%SWING LEG
Ps = ms*g*dot(rMs,jhat);

% Compute entire system energy 
T = simplify(T0 + T1 + T2 + T3 + T4 + Ts + T1r + T2r + Tsr); %with x movement
V = P1+P2+P3+P4+Pky+Pkx + Ps;

% Find Generalized forces
Q_tau1 = M2Q(tau1*khat,dth1*khat);
Q_tau2 = M2Q(tau2*khat,dth2*khat);

%swing
Q_taus = M2Q(taus*khat,dths*khat);
% Q_tau2R= M2Q(-tau2*khat,omega1*khat); %WHAT IS THIS

Q = Q_tau1 + Q_tau2 + Q_taus;

% Assemble the array of cartesian coordinates of the key points
keypoints = [r0(1:2) rA(1:2) rB(1:2) rC(1:2) rD(1:2) rI(1:2) r_heela(1:2) rHs(1:2)];

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

%hip position
Jhip = jacobian(r0,q);
dJhip= reshape( ddt(Jhip(:)) , size(Jhip) );

matlabFunction(r0,'file',['position_hip'],'vars',{z p}); 
matlabFunction(vR0,'file',['velocity_hip'],'vars',{z p});
matlabFunction(Jhip,'file',['jacobian_hip'],'vars',{z p});
matlabFunction(dJhip ,'file',['jacobian_dot_hip'],'vars',{z p});

matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});


%Swing leg 
Jswing = jacobian(rHs,q);
dJswing= reshape( ddt(Jswing(:)) , size(Jswing) );

% Jswing  = Jswing(1:2,1:6);
% dJswing = dJswing(1:2,1:6);

matlabFunction(rHs,'file',['position_toe_swing'],'vars',{z p}); 
matlabFunction(drHs,'file',['velocity_toe_swing'],'vars',{z p});
matlabFunction(Jswing,'file',['jacobian_toe_swing'],'vars',{z p});
matlabFunction(dJswing ,'file',['jacobian_dot_swing'],'vars',{z p});

matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
