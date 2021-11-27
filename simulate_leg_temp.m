function simulate_leg()
    %% Definte fixed paramters
    m0 = 0.03; %MASS OF BOOM, ETC
    m1 =.0393 + .2;         m2 =.0368; 
    m3 = .00783;            m4 = .0155;
    m5 = .00283;
    
    I1 = 25.1 * 10^-6;      I2 = 53.5 * 10^-6;
    I3 = 9.25 * 10^-6;      I4 = 22.176 * 10^-6;
    I5 = 1876 * 10^-9;
    
    l_OA=.011;              l_OB=.042; 
    l_AC=.096;              l_DE=.091;
    l_IG = 0.01;            l_GH = 0.03;            l_heela = 0.09;
    l_O_m1=0.032;           l_B_m2=0.0344; 
    l_A_m3=0.0622;          l_C_m4=0.0610;
    c5=0.01;
    
    %SWING LEG
    ms = .05;
    Is = 25.1 * 10^-6;
    ls = .13;
    l_C_s = .06;
    
    N = 18.75;
    Ir = 0.0035/N^2;
    g = 9.81;
    %=0;
    k=1;   
    
    restitution_coeff = 0.5;
    friction_coeff = 0.01;
    ground_height = 0;
    %% Parameter vector
    p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB l_AC l_DE l_IG l_GH l_heela ls l_C_s g k]';        % parameters
    
    
    %% Simulation Parameters Set 2 -- Operational Space Control
    p_traj.omega = 15;
    p_traj.x_0   = -0.02;
    p_traj.y_0   = 0;
    p_traj.r     = 0.025;
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 5;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [0; 0.2; -pi/3; pi/2; 0; 0; 0; 0; 0; 0; 0; 0];
    z_out = zeros(12,num_step);
    z_out(:,1) = z0;
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p, p_traj);
        % Velocity update with dynamics
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        
        % constraint handling (Velocity update)
        %z_out(3:4,i+1) = joint_limit_constraint(z_out(:,i+1),p);
        z_out(7:12,i+1) = discrete_impact_contact(z_out(:,i+1), p, restitution_coeff, friction_coeff, ground_height);
        
        % Position update
        z_out(1:6,i+1) = z_out(1:6,i) + z_out(7:12,i+1)*dt;
    end
    
    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time
    rE = zeros(3,length(tspan));
    rH = zeros(3,length(tspan));
    rI = zeros(3,length(tspan));
    vE = zeros(3,length(tspan));
    vH = zeros(3,length(tspan));
    vI = zeros(3,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_ankle(z_out(:,i),p);
        rH(:,i) = position_toe(z_out(:,i),p);
        rI(:,i) = position_heel(z_out(:,i),p);
        vE(:,i) = velocity_ankle(z_out(:,i),p);
        vH(:,i) = velocity_toe(z_out(:,i),p);
        vI(:,i) = velocity_heel(z_out(:,i),p);
    end
    
    %PLOT DESIRED POSITION VS ACTUAL ANKLE POSITION
    figure(2); clf;
    plot(tspan,rE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,p_traj.x_0 + p_traj.r * cos(p_traj.omega*tspan) ,'r--');
    plot(tspan,rE(2,:),'b','LineWidth',2)
    plot(tspan,p_traj.y_0 + p_traj.r * sin(p_traj.omega*tspan) ,'b--');
    
    
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','x_d','y','y_d'});

    %PLOT TIME VS VELOCITY
    figure(3); clf;
    plot(tspan,vE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,vE(2,:),'b','LineWidth',2)
    
    xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x','vel_y'});
    
    %PLOT ANGLE 
    figure(4)
    plot(tspan,z_out(1:2,:)*180/pi)
    legend('q1','q2');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    
    %PLOT ANGULAR VELOCITY
    figure(5)
    plot(tspan,z_out(3:4,:)*180/pi)
    legend('q1dot','q2dot');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    
    %% Animate Solution
    figure(6); clf;
    hold on
   
    %% plot foot target information

    % Target traj
    TH = 0:.1:2*pi;
    plot( p_traj.x_0 + p_traj.r * cos(TH), ...
          p_traj.y_0 + p_traj.r * sin(TH),'k--'); 
    
    % Ground Q2.3
    plot([-.2 .2],[ground_height ground_height],'k'); 
    
    animateSol(tspan, z_out,p);
end

function tau = control_law(t, z, p, p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 150.; % Spring stiffness X
    K_y = 150.; % Spring stiffness Y
    D_x = 10.;  % Damping X
    D_y = 10.;  % Damping Y

    % Desired position of foot is a circle
    % ONLY CONTROLLING HIP AND ANKLE
    omega_swing = p_traj.omega;
    rEd = [p_traj.x_0 p_traj.y_0 0]' + ...
            p_traj.r*[cos(omega_swing*t) sin(omega_swing*t) 0]';
    % Compute desired velocity of foot
    vEd = p_traj.r*[-sin(omega_swing*t)*omega_swing    ...
                     cos(omega_swing*t)*omega_swing   0]';
    % Desired acceleration
    aEd = p_traj.r*[-cos(omega_swing*t)*omega_swing^2 ...
                    -sin(omega_swing*t)*omega_swing^2 0]';
    
    % Actual position and velocity 
    rE = position_foot(z,p);
    vE = velocity_foot(z,p);
    
    % Jacobian matrix \partial r_E / \partial q
    J  = jacobian_foot(z,p);
    dJ = jacobian_dot_foot(z,p);
    dq = z(7:12);

    % Compute virtual foce for Question 1.4 and 1.5
    f  = [K_x * (rEd(1) - rE(1) ) + D_x * (vEd(1) - vE(1) ) ;
          K_y * (rEd(2) - rE(2) ) + D_y * (vEd(2) - vE(2) ) ; 0 ; 0 ];
    
    %% Task-space compensation and feed forward for Question 1.8
    % Get joint space components of equations of motion
    Mass_Joint_Sp = A_leg(z,p);
    Grav_Joint_Sp = Grav_leg(z,p);
    Corr_Joint_Sp = Corr_leg(z,p);

    Mass_Joint_Sp_inv = inv(Mass_Joint_Sp);
    % Task-space mass matrix (Equaiton 51 in Khatib's paper)
    Lambda = inv(J * Mass_Joint_Sp_inv * J');
    
    % Coriolis force in task-space (Equation 51)
    mu     = Lambda*J*Mass_Joint_Sp_inv* Corr_Joint_Sp - Lambda * dJ * dq;
    
    % Gravity force in task-space (Equation 51)
    rho    = Lambda*J*Mass_Joint_Sp_inv * Grav_Joint_Sp; 
    
    % Add task-space acceleration force feedforward, coriolis, and gravity compensation 
    f(1:6) = Lambda*(aEd(1:2) + f(1:2)) + mu + rho; % OSC
%     f(1:2) = Lambda*(aEd(1:2) + f(1:2)) + rho; % OSC w/o mu (coriolis)
%     f(1:2) = Lambda*(aEd(1:2) + f(1:2)) + mu; % OSC w/o rho (gravity)
    
    % Map to joint torques  
    tau = J' * f;
end


function dz = dynamics(t,z,p,p_traj)
    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p,p_traj);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:6) = z(7:12);
    dz(7:12) = qdd;
end


function qdot = joint_limit_constraint(z,p)
    %% 
    q1_min = -50 * pi/ 180;
    C = z(1) - q1_min; % C gives distance away from constraint
    dC= z(3);
    qdot = z(3:4);
    J = [1 0];
    A = A_leg(z,p);

    if (C < 0 && dC <0)% if constraint is violated
        lambda = A(2,2);
        F_c = lambda * (0 - dC);
        qdot = qdot + inv(A)*J.'*F_c;        
    end
end

function qdot = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)

    qdot = z(7:12);
    rToe = position_toe(z,p); 
    rHeel = position_heel(z,p);
    vToe  = velocity_toe(z,p);  
    vHeel = velocity_heel(z,p); 

    %TOE CONTACT
    if(rToe(2)-yC < 0 && vToe(2) < 0)
      J  = jacobian_toe(z,p);
      A = A_leg(z,p);
      Ainv = inv(A);
      
      J_z = J(2,:);
      lambda_z = 1/(J_z * Ainv * J_z.');
      F_z = lambda_z*(-rest_coeff*vToe(2) - J_z*qdot);
      qdot = qdot + Ainv*J_z.'*F_z;
      
      % horizontal
      J_x = J(1,:);
      lambda_x = 1/(J_x * Ainv * J_x.');
      F_x = lambda_x * (0 - J_x * qdot);
      if( abs(F_x) > fric_coeff*F_z)
          F_x = sign(F_x)*F_z*fric_coeff;
      end
      qdot = qdot + Ainv*J_x.'*F_x;
    z_test = z;
    z_test(7:12) = qdot;
    vToe = velocity_toe(z_test, p);
    end
    
    %HEEL CONTACT
    if(rHeel(2)-yC < 0 && vHeel(2) < 0)
      J  = jacobian_heel(z,p);
      A = A_leg(z,p);
      Ainv = inv(A);
      
      J_z = J(2,:);
      lambda_z = 1/(J_z * Ainv * J_z.');
      F_z = lambda_z*(-rest_coeff*vHeel(2) - J_z*qdot);
      qdot = qdot + Ainv*J_z.'*F_z;
      
      % horizontal
      J_x = J(1,:);
      lambda_x = 1/(J_x * Ainv * J_x.');
      F_x = lambda_x * (0 - J_x * qdot);
      if( abs(F_x) > fric_coeff*F_z)
          F_x = sign(F_x)*F_z*fric_coeff;
      end
      qdot = qdot + Ainv*J_x.'*F_x;
    z_test = z;
    z_test(7:12) = qdot;
    vToe = velocity_toe(z_test, p);
    end
end

function animateSol(tspan, x,p)
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_IG = plot([0],[0],'LineWidth',2);
    h_GH = plot([0],[0],'LineWidth',2);
    %THIS IS OUR ANKLE SPRING
    h_ankle = plot([0],[0],'LineWidth',0.5);
    
    %swing leg
    h_swing = plot([0],[0],'LineWidth',0.5);

    %line to indicate the ground
    ground  = plot([0],[0],'LineWidth',0.3);
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 0.5 -.2 .5]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,30)
            continue
        end
        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_leg(z,p);

        r0 = keypoints(:,1);
        rA = keypoints(:,2); % Vector to base of cart
        rB = keypoints(:,3);
        rC = keypoints(:,4); % Vector to tip of pendulum
        rD = keypoints(:,5);
        rE = keypoints(:,6);
        rI = keypoints(:,7);
        rH = keypoints(:,8);
        r_heela = keypoints(:,9);
        r_swing = keypoints(:,10);

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[r0(1) rB(1)]); %should be x rB
        set(h_OB,'YData',[r0(2) rB(2)]); %DO WE HAVE Y MOOVEMENT? HIP MIGHT NATURALLY RISE AND FALL
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);
        
        %FOOT
         set(h_IG,'XData',[rE(1) rI(1)]);
         set(h_IG,'YData',[rE(2) rI(2)]);
        
        set(h_GH,'XData',[rE(1) rH(1)]);
        set(h_GH,'YData',[rE(2) rH(2)]);
        
        %ANKLE
        set(h_ankle, 'XData' , [r_heela(1) rI(1)] );
        set(h_ankle, 'YData' , [r_heela(2) rI(2)] );
        
        %SWING LEG
        set(h_swing, 'XData' , [r0(1) r_swing(1)] );
        set(h_swing, 'YData' , [r0(2) r_swing(2)] );

        %Ground
        set(ground, 'XData' , [-2 2] );
        set(ground, 'YData' , [0 0] );
        pause(.05)
    end
end