function simulate_system_leg()

    %% Definte fixed paramters
    m0 = 0.03; %MASS OF BOOM, ETC
    m1 =.0393 + .2;         m2 =.0368; 
    m3 = .00783;            m4 = .0155;
    m5 = .00283;
    
    I1 = 25.1 * 10^-6;      I2 = 53.5 * 10^-6;
    I3 = 9.25 * 10^-6;      I4 = 22.176 * 10^-6;
    I5 = 2.25 * 10^-6;
    
    l_OA=.011;              l_OB=.042; 
    l_AC=.096;              l_DE=.091;
    l_IG = 0.01;            l_GH = 0.03;            l_heela = 0.09;
    l_O_m1=0.032;           l_B_m2=0.0344; 
    l_A_m3=0.0622;          l_C_m4=0.0610;
    c5=0.01;
    
    N = 18.75;
    Ir = 0.0035/N^2;
    g = 9.81;
    k=0.001;
    
    %% Parameter vector
    p   = [m0 m1 m2 m3 m4 m5 I1 I2 I3 I4 I5 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB l_AC l_DE l_IG l_GH l_heela g k]';        % parameters
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 1;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [0; -pi/4; pi/2; 0; 0; 0; 0; 0];
    z_out = zeros(8,num_step);
    z_out(:,1) = z0;
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p);
        z_out(3:4,i+1) = z_out(3:4,i) + dz(3:4)*dt;
        z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt;
    end

    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time
    rH = zeros(3,length(tspan));
    for i = 1:length(tspan)
        rH(:,i) = position_foot(z_out(:,i),p);
    end
    figure(2); clf;
    plot(tspan,rH(1:2,:))
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','y'});
    
    %hold on
    %w = 3;
    %Ed = [0.025*cos(w*tspan); -0.125+0.025*sin(w*tspan)];
    %plot(tspan,Ed(1:2,:))
    %% Animate Solution
    figure(3); clf;
    hold on
   
    %% Optional, plot foot target information
    
    % target position Q1.4
    plot(.025 , -.125, 'r.','MarkerSize',6); 
    
    % Target traj. Q 1.6
    %plot( .025*cos(0:.01:2*pi), -.125+.025*sin(0:.01:2*pi),'k--'); 
    
    % Ground Q2.3
    %plot([-.2 .2],[-.125 -.125],'k'); 
    
    animateSol(tspan,z_out,p);
end

function tau = control_law(t,z,p)
    % Controller gains, Update as necessary for Problem 1
    K_x = 40; % Spring stiffness X
    K_y = 40; % Spring stiffness Y
    D_x = 4;  % Damping X
    D_y = 4;  % Damping Y
    rHd = [.025 -.2 0]'; % Desired position of foot
    w = 3;
    
    %% STEPS TO COMPLETE PROBLEM 1.3
    % a. Compute r_E
   rH = position_foot(z,p);
    % b. Compute J, the jacobian of r_E
   jH = jacobian_foot(z,p);
    % c. Use these results to compute \tau as specified in the write-up
   V_virt= 1/2*K_x*(rH(1)-rHd(1))^2+1/2*K_y*(rH(2)-rHd(2))^2; %DESIRED TOE POSITION VS ACTUAL
    % Hint: Some of the functions generated by derive_leg will help with 
    % steps a and b.
   vF = velocity_foot(z,p);
   rHd = [0.025*cos(w*t); -0.2+0.025*sin(w*t)];
   
   tau = -transpose(jH)*[K_x*(rH(1)-rHd(1))+D_x*vF(1);K_y*(rH(2)-rHd(2))+D_y*vF(2);0];  %DAMPING
end


function Fc = contact_force(z,p)

    %% Fixed parameters for contact
    K_c = 100;
    D_c = 2;
    yC  = -.2;
    
    %% STEPS TO COMPLETE PROBLEM 2.1
    % a. Compute constraint C which gives height of foot relative to ground
    % b. Compute constraint rate, \dot{C}
    % c. Set Fc based on compliant contact model
    % d. If foot is above the ground, or Fc<0, set Fc = 0
    
    % Hint: Some of the functions generated by derive_leg will help with 
    % steps a and b.
    rH = position_foot(z,p);
    vF = velocity_foot(z,p);
    
    C = rH(2)-yC;
    C_dot = vF(2);
   
    %If foot in contact with the ground
    if C>0
        Fc = 0;
    else
        Fc = -K_c*C-D_c*C_dot; % replace this line using steps a-d to compute Fc
    end
end


function Tauc = joint_limit_torque(z,p)
    %% Fixed parameters for rotational spring damper at joint limit contact
    Kappa_c = 10;
    Dampa_c = 0.2;

    Tauc = 0;
end


function dz = dynamics(t,z,p)
    x = z(1);        th1 = z(2);     
    th2 = z(3);      th3 = z(4);
    
    dx = z(5);       dth1 = z(6);     
    dth2= z(7);      dth3 = z(8);
    
    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p);
%     tau = [0;0;0;0]; 
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
  
    % Compute the contact force (used for problem 2)
    Fc = contact_force(z,p);
   
    J = jacobian_foot(z,p);
    % Compute the contribution of the contact force to the generalied force
    QFc= transpose(J)*[0;Fc;0];  %% YOUR CODE HERE for Q2.2
%     QFc = zeros([4,1]);

    % Compute the contact force (used for problem 2.5)
    Tauc = joint_limit_torque(z,p);
    QTauc= [0; 0; 0; 0]; %column vector
    
    % Solve for qdd.
    qdd = A\(b + QFc + QTauc);
    dz = 0*z;
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
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
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 1 -.3 .2]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
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

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[r0(1) rB(1)]); %should be x rB
        set(h_OB,'YData',[0 rB(2)]); %DO WE HAVE Y MOOVEMENT? HIP MIGHT NATURALLY RISE AND FALL
        
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

        pause(.01)
    end
end