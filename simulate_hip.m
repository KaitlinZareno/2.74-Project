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
    l_AC=.096;                l_DE=.0897;
    l_CE=.122;
    l_IG = 0.01;            l_GH = 0.03;            l_heela = 0.03;
    l_m1=0.032;             l_m2=0.0344; 
    l_m3=0.0622;            l_m4=0.0610;
    c5=0.01;
    
    %SWING LEG
    ms = .05;
    Is = 25.1 * 10^-6;
    ls = .13;
    l_ms = .06;
    
    N = 18.75;
    Ir = 0.0035/N^2;
    g = -9.81;
    %=0;
    k=1;   
    
    restitution_coeff = 0.5;
    friction_coeff = 0.01;
    ground_height = 0;
    %% Parameter vector
    p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N g k l_m1 l_m2 l_m3 l_m4 l_ms ls l_OA l_OB l_AC l_DE l_CE l_IG l_GH l_heela]';        % parameters
    
    
    %% Simulation Parameters Set 2 -- Operational Space Control
    p_traj.omega = 15;
    p_traj.x_0   = 0.0;
    p_traj.y_0   = .17;
    p_traj.rx    = 0.1;
    p_traj.ry    = 0.015;
    
    %% Perform Dynamic simulation
    dt = 0.0001;
    tf = 5;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [-pi/2; pi/2; 3*pi/4; 0; 0; 0];
    z_out = zeros(6,num_step);
    z_out(:,1) = z0;
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p, p_traj);
        % Velocity update with dynamics
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        
        % constraint handling (Velocity update)
        %z_out(3:4,i+1) = joint_limit_constraint(z_out(:,i+1),p);
        %z_out(7:12,i+1) = discrete_impact_contact(z_out(:,i+1), p, restitution_coeff, friction_coeff, ground_height);
        
        % Position update
        z_out(1:3,i+1) = z_out(1:3,i) + z_out(4:6,i+1)*dt;
    end
    
    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time
    r0 = zeros(3,length(tspan));
    v0 = zeros(3,length(tspan));
    for i = 1:length(tspan)
        r0(:,i) = position_hip(z_out(:,i),p);
        v0(:,i) = velocity_hip(z_out(:,i),p);
    end
    
    %PLOT DESIRED HIP POSITION VS ACTUAL 
    figure(2); clf;
    plot(tspan,r0(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,p_traj.x_0 + p_traj.rx * cos(p_traj.omega*tspan) ,'r--');
    plot(tspan,r0(2,:),'b','LineWidth',2)
    plot(tspan,p_traj.y_0 + p_traj.ry * sin(p_traj.omega*tspan) ,'b--');
    
    
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','x_d','y','y_d'});

    %PLOT TIME VS VELOCITY
    figure(3); clf;
    plot(tspan,v0(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,v0(2,:),'b','LineWidth',2)
    
    xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x','vel_y'});
    
    %PLOT ANGLE 
    figure(4)
    plot(tspan,z_out(1:2,:)*180/pi)
    legend('q1','q2');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    
    %PLOT ANGULAR VELOCITY
    figure(5)
    plot(tspan,z_out(3:4,:)*180/pi) %LINE MIGHT BE 4:6
    legend('q1dot','q2dot');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    
    %% Animate Solution
    figure(6); clf;
    hold on
   
    %% plot foot target information

    % Target traj
    TH = 0:.1:2*pi;
    plot( p_traj.x_0 + p_traj.rx * cos(TH), ...
          p_traj.y_0 + p_traj.ry * sin(TH),'k--'); 
    
    % Ground Q2.3
    plot([-.2 .2],[ground_height ground_height],'k'); 
    
    animateSol(tspan, z_out,p);
end

function tau = control_law(t, z, p, p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 150.; % Spring stiffness X
    K_y = 150.; % Spring stiffness Y
    K_s = 200.;
    D_x = 10.;  % Damping X
    D_y = 10.;  % Damping Y
    D_s = 30.;
    
    pos_hip = position_hip(z,p);
    v_hip = velocity_hip(z,p);
    
    v1 = norm(pos_hip(1:2) - [0 0]);
    
    v2 = norm(pos_hip(2));
    desired_angle = acos(v2/v1);

    % Desired position of foot is a circle
    % ONLY CONTROLLING HIP AND ANKLE
    omega_swing = p_traj.omega;
    r0d = [p_traj.x_0 + p_traj.rx*cos(omega_swing*t) ...
           p_traj.y_0 + p_traj.ry*sin(omega_swing*t) ...
           0];
%     % Compute desired velocity of foot
    v0d = [p_traj.rx*-sin(omega_swing*t)*omega_swing ...
           p_traj.ry* cos(omega_swing*t)*omega_swing   0]';

    r0 = position_hip(z,p);
    v0 = velocity_hip(z,p);
    
    if pos_hip(1) < 0
        tds = 3*pi/2-desired_angle; %WHY THESE VALUES??
    else
        tds = 0-desired_angle;
    end
    
    if abs(pos_hip(1)) > 0.09
        wd= 0;
    else
        wd = v0d(2);
    end
    
    ts = z(3);
    vs = z(6);
    
    tau = [K_x * (r0d(1) - r0(1) ) + D_x * (v0d(1) - v0(1) ) ;
           K_y * (r0d(2) - r0(2) ) + D_y * (v0d(2) - v0(2) ); K_s * (tds - ts) + D_s * (wd-vs)];
      
end


function dz = dynamics(t,z,p,p_traj)
    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p,p_traj);
%     tau = [0;0;0];
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:3) = z(4:6);
    dz(4:6) = qdd;
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
    z_test(1:12) = z
    vToe = velocity_toe(z_test, p);
    
    %HEEL CONTACT
    vHeel = velocity_heel(z_test, p);
end

function animateSol(tspan, x,p)
    % Prepare plot handles
    hold on
    
    h_IG = plot([0],[0],'LineWidth',2);
    h_GH = plot([0],[0],'LineWidth',2);
    %THIS IS OUR ANKLE SPRING
    h_ankle = plot([0],[0],'LineWidth',0.5);
      
    h_EC = plot([0],[0],'LineWidth',2);
    h_CA = plot([0],[0],'LineWidth',2);
    h_DB = plot([0],[0],'LineWidth',2);
    h_B0 = plot([0],[0],'LineWidth',2);
    
    %swing leg
    h_swing = plot([0],[0],'LineWidth',0.5);

    %line to indicate the ground
    ground  = plot([0],[0],'LineWidth',0.3);
   
    plot(0.05,0.17,'o')
    
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
        rI = keypoints(:,6);
        %rH = keypoints(:,7);
        r_heela = keypoints(:,7);
        r_swing = keypoints(:,8);

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        %FOOT
        %ankle at position(0,0)
        set(h_IG,'XData',[0 -0.01]);
        set(h_IG,'YData',[0 0]);
        
        set(h_GH,'XData',[0 0.03]);
        set(h_GH,'YData',[0 0]);
        
        %ANKLE
        set(h_ankle, 'XData' , [-0.01 r_heela(1)] );
        set(h_ankle, 'YData' , [0 r_heela(2)] );
       
        set(h_EC,'XData',[0 rC(1)]); %should be x rB
        set(h_EC,'YData',[0 rC(2)]); %DO WE HAVE Y MOOVEMENT? HIP MIGHT NATURALLY RISE AND FALL
        
        set(h_DB,'XData',[rD(1) rB(1)]);
        set(h_DB,'YData',[rD(2) rB(2)]);
        
        set(h_CA,'XData',[rC(1) rA(1)]);
        set(h_CA,'YData',[rC(2) rA(2)]);
        
        set(h_B0,'XData',[rB(1) r0(1)]);
        set(h_B0,'YData',[rB(2) r0(2)]);
        
        
        %SWING LEG
        set(h_swing, 'XData' , [r0(1) r_swing(1)] );
        set(h_swing, 'YData' , [r0(2) r_swing(2)] );

        %Ground
         set(ground, 'XData' , [-2 2] );
         set(ground, 'YData' , [0 0] );
         pause(.02)
    end
end