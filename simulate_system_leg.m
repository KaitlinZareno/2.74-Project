function simulate_system_leg()

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
    %g=0;
    k=1;
    
    restitution_coeff = 0.;
    friction_coeff = 1000;
    ground_height = 0;
    
    %% Parameter vector
    p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB l_AC l_DE l_IG l_GH l_heela ls l_C_s g k]';        % parameters

    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 3;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [0; .25; -pi/4; pi/2; 0; 0; 0; 0; 0; 0; 0; 0]; %moves with x velocity!
    z_out = zeros(12,num_step);
    z_out(:,1) = z0;
    %debugging
    vx = zeros(1,num_step); 
    condition = zeros(1,num_step); 
    
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p);
        z_out(:,i+1) = z_out(:,i) + dz*dt; 
        %discrete impact
         [z_out(7:12,i+1), condition(i)] = discrete_impact_contact(z_out(:,i+1),p, restitution_coeff, friction_coeff, ground_height);
        z_out(1:6,i+1) = z_out(1:6,i) + z_out(7:12,i+1)*dt; %CHECK
        v = velocity_toe(z_out(:,i),p); 
        
        %debugging statements
        vx(i) = v(1); 
        disp(num2str(i*dt)); 
    end

    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Check forces and vel over time
    figure(6); clf;
    plot(tspan, vx);
    xlabel('time'); 
    ylabel('VxToe'); 
    
    %% Check the contact condition over time
    figure(7); clf;
    plot(tspan, condition);
    xlabel('time'); 
    ylabel('Contact Condition'); 
    
    %% Compute foo, heel and ankle position over time
    rH = zeros(3,length(tspan));
    rI = zeros(3,length(tspan));
    rE = zeros(3,length(tspan));
    for i = 1:length(tspan)
        rH(:,i) = position_toe(z_out(:,i),p);
        rI(:,i) = position_heel(z_out(:,i),p);
        rE(:,i) = position_ankle(z_out(:,i),p);
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

function [qdot, c] = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)

    %who is in contact w/ ground?
    %TOE CONTACT
    y = position_toe(z,p); 
    y = -y(2);
    y_dot = velocity_toe(z,p);
    Cy = y-yC;
    Cy_dot = y_dot(2);
    
    %HEEL CONTACT
    y2 = position_heel(z,p); 
    y2 = -y2(2);
    y_dot2 = velocity_heel(z,p);
    Cy2 = y2-yC;
    Cy_dot2 = y_dot2(2); 
  
        
    if (Cy <= 0 && Cy_dot <= 0)
        c =1; 
        disp('toe touching');
        
        J = jacobian_toe(z,p);
        Mm = A_leg(z,p);

        lamY = inv(J(2,:)*inv(Mm)*J(2,:)');
        lamX = inv(J(1,:)*inv(Mm)*J(1,:)');
        
        %update y dir
        F_cy = lamY*(-rest_coeff*Cy_dot-(J(2,:)*z(7:12)));
        z(7:12) = z(7:12) + inv(Mm)*J(2,:)'*F_cy;
        %update x dir
        F_cx = lamX*(0-J(1,:)*z(7:12));
        %truncate if necessary
        if (F_cx > fric_coeff*F_cy)
            F_cx = fric_coeff*F_cy;
        elseif (-F_cx < -fric_coeff*F_cy)
            F_cx = -fric_coeff*F_cy;
        end
        z(7:12) = z(7:12) + inv(Mm)*J(1,:)'*F_cx; 
        qdot = z(7:12);
    elseif (Cy2 <= 0 && Cy_dot2 <= 0)
        c =2;
        disp('heel touching');
        J = jacobian_heel(z,p);
        Mm = A_leg(z,p);

        lamY = inv(J(2,:)*inv(Mm)*J(2,:)');
        lamX = inv(J(1,:)*inv(Mm)*J(1,:)');
        
        %update y dir
        F_cy = lamY*(-rest_coeff*Cy_dot2-(J(2,:)*z(7:12)));
        z(7:12) = z(7:12) + inv(Mm)*J(2,:)'*F_cy;
        %update x dir
        F_cx = lamX*(0-J(1,:)*z(7:12));
        %truncate if necessary
        if (F_cx > fric_coeff*F_cy)
            F_cx = fric_coeff*F_cy;
        elseif (-F_cx < -fric_coeff*F_cy)
            F_cx = -fric_coeff*F_cy;
        end
        z(7:12) = z(7:12) + inv(Mm)*J(1,:)'*F_cx; 
        qdot = z(7:12);
    else
        c=3; 
        disp('no contact'); 
        qdot = z(7:12); 
    end
end

function tau = control_law(t,z,p)
    % Controller gains, Update as necessary for Problem 1
    K_x = 20; % Spring stiffness X
    K_y = 20; % Spring stiffness Y
    D_x = 5;  % Damping X
    D_y = 5;  % Damping Y
    rEd = [0 0 0]'; % Desired position of ANKLE CHANGW
    w = 3;
    
    %% STEPS TO COMPLETE PROBLEM 1.3
    % a. Compute r_E
    
    %ONLY CONTROLLING ANKLE
   rE = position_ankle(z,p); %THIS VALUE NOT HAPPY
    % b. Compute J, the jacobian of r_E
   jE = jacobian_ankle(z,p); %THIS VALUE NOT HAPPY
   jE = jE(:, 4:end); % # OF CONTROLS?
   vF = velocity_ankle(z,p);
   rEd = [0.05; 0.05; 0 ; 0]; %TRACK SINUSOID , make point g move in an ellipse

   tau = -transpose(jE)*[K_x*(rE(1)-rEd(1))+D_x*vF(1);K_y*(rE(2)-rEd(2))+D_y*vF(2); 0 ; 0];  %DAMPING
end


function Fc = contact_force_toe(z,p)

    %% Fixed parameters for contact
    K_c = 2*10^4;
    D_c = 5;
    yC  = 0;
    
    %% STEPS TO COMPLETE PROBLEM 2.1
    % a. Compute constraint C which gives height of foot relative to ground
    % b. Compute constraint rate, \dot{C}
    % c. Set Fc based on compliant contact model
    % d. If foot is above the ground, or Fc<0, set Fc = 0
    
    % Hint: Some of the functions generated by derive_leg will help with 
    % steps a and b.
    rH = position_toe(z,p);         toe = rH(2)-yC;
    vH = velocity_toe(z,p);         toe_dot = vH(2);  
   
    %If toe in contact with the ground
    if toe>0
        Fc = 0;
    else
        Fc = -K_c*toe - D_c*toe_dot; % replace this line using steps a-d to compute Fc
    end
end

function Fh = contact_force_heel(z,p)

    %% Fixed parameters for contact
    K_c = 2*10^4;
    D_c = 5;
    yC  = 0;
    
    %% STEPS TO COMPLETE PROBLEM 2.1
    % a. Compute constraint C which gives height of foot relative to ground
    % b. Compute constraint rate, \dot{C}
    % c. Set Fc based on compliant contact model
    % d. If foot is above the ground, or Fc<0, set Fc = 0
    
    % Hint: Some of the functions generated by derive_leg will help with 
    % steps a and b.    
    rI = position_heel(z,p);        heel = rI(2)-yC;
    vI = velocity_heel(z,p);        heel_dot = vI(2);
   
    %if heel in contact with the ground
    if heel>0
        Fh = 0;
    else
        Fh = -K_c*heel - D_c*heel_dot; %rcross f
    end 
end

function Ankle_Res = ankle_constraint(z,p)
    %make sure the position of the heel is always to the left of that of
    %the ankle -- can't go past the vector of the ankle

    th3 = z(5);
    th1 = z(3);
    if th3 <= 0
        th3=0.;
    end
    if th3 > pi/8
        th3= pi/8;
    end
	Ankle_Res = th3;
end


function Tauc = joint_limit_torque(z,p)
    %% Fixed parameters for rotational spring damper at joint limit contact
    Kappa_c = 10;
    Dampa_c = 0.2;

    Tauc = 0;
end


function dz = dynamics(t,z,p)
%     x = z(1);        y= z(2)      
%    th1 = z(3);     th2 = z(4);      th3 = z(5);
%     
%     dx = z(5);       dy = z(6);     
%     dth1 = z(7);  dth2= z(8);      dth3 = z(9);
    
    % Get mass matrix
    A = A_leg(z,p)
    
    % Compute Controls
    %tau = control_law(t,z,p); %tau EXPLODING WITH CONTROLS
    tau = [0;0;0;0]; 
    
%     % Get b = Q - V(q,qd) - G(q)
%     b = b_leg(z,tau,p); 
  
    % Compute the contact force (used for problem 2)
    Fc = contact_force_toe(z,p);
    Fh = contact_force_heel(z,p);
   
    J_toe = jacobian_toe(z,p);
    J_heel = jacobian_heel(z,p);
    % Compute the contribution of the contact force to the generalied force
    QFc_t=  transpose(J_toe)*[0;Fc;0];
    QFc_h=  transpose(J_heel)*[0;Fh;0];
    QFc = QFc_t +  QFc_h; 
    
    z(5) = ankle_constraint(z,p);
    %controls
    t1= z(3);
    t2 = z(4);
    td1 = -pi/4;
    td2 = pi/2;
    ts = z(6);
    tsd = 0;
    
    K= 100;
    D=5;
    tau = [K*(td1-t1)+D*-t1;K*(td2-t2)+D*-t2; 0 ; K*(tsd-ts)+D*-ts];
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p)
    
    % Compute the contact force (used for problem 2.5)
    Tauc = joint_limit_torque(z,p);
    QTauc= [0; 0; 0; 0;0;0]; %column vector
    
    % Solve for qdd.
    qdd = A\(b + QFc + QTauc);
    dz = 0*z;
    
    % Form dz
    dz(1:6) = z(7:12);
    dz(7:12) = qdd;
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
    
    h_swing = plot([0],[0],'LineWidth',0.5);
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 0.5 -.2 .5]);

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

        pause(.05)
    end
end