function [CoT] = simulate_system_leg_downloaded2(rx, ry, k)

    global ptraj
    
    %% Definte fixed paramters
    m0 = 1; %MASS OF BOOM, ETC
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
    l_A_m3=0.0622;          l_C_m4=0.0610;          l_anklerest = 0.08;
    c5=0.01;
    
    R_motor = 1.46;
    kt = .18;
    
    %SWING LEG
    ms = .05;
    Is = 25.1 * 10^-6;
    ls = .13;
    l_C_s = .06;
    
    m_total = m0 + 2*(m1+m2+m3+m4+m5);
    
    N = 18.75;
    Ir = 0.0035/N^2;
    g = 9.81;
    k = k;
    
    restitution_coeff = 0.;
    friction_coeff = 10000;
    ground_height = 0;
    
    %% Ellipse param
    ptraj.x0 = 0; 
    ptraj.y0 = .16; 
    ptraj.rx = rx; 
    ptraj.ry = ry;
    ptraj.omega = -pi; 
    ptraj.phase = pi; 
    
    %% Parameter vector
    p   = [m0 m1 m2 m3 m4 m5 ms I1 I2 I3 I4 I5 Is Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 c5 l_OA l_OB ...
       l_AC l_DE l_IG l_GH l_heela ls l_C_s g k l_anklerest]';        % parameters
    
    %% Perform Dynamic simulation - Springy ground
    tf = abs(2*pi/ptraj.omega)*2;
    tspan = [0 tf];
    
    z0 = [0; .186; -pi/4; pi/2; 0;-pi/4; pi/2; 0 ...
        ; 0; 0; 0; 0; 0; 0; 0; 0]; %moves with x velocity!
%     opts = odeset('AbsTol',1e-8,'RelTol',1e-6);
%     sol = ode45(@dynamics,tspan,z0,opts,p);
%     
%     %bc lazy and code is copied from everywhere
%     %change variable names for graphs
%     z_out = sol.y;
%     tspan = sol.x; 
    
    %% Perform Dynamic Sim - discrete contact
    dt = .001;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z_out = zeros(16,num_step);
    z_out(:,1) = z0;
    tau = zeros(6,num_step);
    
    for i=1:num_step-1
        [dz, tau(:,i+1)] = dynamics(tspan(i), z_out(:,i), p);
        % Velocity update with dynamics
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        
        % constraint handling (Velocity update)
        
        z_out(9:16,i+1) = discrete_impact_contact(z_out(:,i+1), p, restitution_coeff, friction_coeff, ground_height);
        
        % Position update
        z_out(1:8,i+1) = z_out(1:8,i) + z_out(9:16,i+1)*dt;
        z_out([5, 8,13, 16],i+1) = ankle_constraint(z_out(:,i+1));
    end
    

    %% Compute Energy
    if 1
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
     
    
    %% Compute foot, heel and ankle position over time
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
%     TH = 0:.1:2*pi;
%     w = 3;
%     xunit = 0.025 * cos(w*TH);
%     yunit = 0.025 * sin(w*TH);
%     plot(xunit, yunit);
%    
    %% Optional, plot foot target information
    
    % target position Q1.4
    %plot(.025 , -.125, 'r.','MarkerSize',6); 
    
    % Target traj. Q 1.6
    %plot( .025*cos(0:.01:2*pi), -.125+.025*sin(0:.01:2*pi),'k--'); 
    
    % Ground Q2.3
    %plot([-.2 .2],[-.125 -.125],'k'); 
    
    animateSol(tspan,z_out,p);
    end
    
    %% calculate cost of transport
    displacement = abs(z_out(1, end)); 
    Power = (tau./kt).^2.*R_motor; 
    Energy = sum(sum(Power*dt)); 
    CoT = Energy/(m_total*g*displacement);

end

function [qdot, c] = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)

    %get the necessary vals
    qdot = z(9:16); 
    rToe = position_toe(z,p); 
    rToe = rToe(2); 
    rHeel = position_heel(z,p);
    rHeel = rHeel(2); 
    vToe  = velocity_toe(z,p);  
    vToex = vToe(1);
    vToe = vToe(2);
    vHeel = velocity_heel(z,p); 
    vHeelx = vHeel(1); 
    vHeel = vHeel(2); 
    
    %for the other leg
    rToe2 = position_toe2(z,p); 
    rToe2 = rToe2(2); 
    rHeel2 = position_heel2(z,p);
    rHeel2 = rHeel2(2); 
    vToe2  = velocity_toe2(z,p);  
    vToex2 = vToe2(1);
    vToe2 = vToe2(2);
    vHeel2 = velocity_heel2(z,p); 
    vHeelx2 = vHeel2(1); 
    vHeel2 = vHeel2(2);

    %who is in contact w/ ground?
    if (rToe < 0 && vToe < 0)
        J  = jacobian_toe(z,p); 
        J2 = J; 
        J = J(2,:); 
        
      M = A_leg(z,p);
      Ainv = inv(M);
      
      lambda_z = 1/(J * Ainv * J.');
        F_z = lambda_z*(0 - J*qdot); 
        qdot = qdot + Ainv*J.'*F_z;
        
        %Horizontal
        J_x = J2(1,:);
        lambda_x = 1/(J_x * Ainv * J_x.');
        F_x = lambda_x * (0 - J_x * qdot);
        qdot = qdot + Ainv*J_x.'*F_x;
    end
        
    if(rHeel < 0 && vHeel < 0)
         J = jacobian_heel(z,p);
         J2 = J; 
         J = J(2,:); 
         
          M = A_leg(z,p);
      Ainv = inv(M);
      
      lambda_z = 1/(J * Ainv * J.');
         
         F_z = lambda_z*(0 - J*qdot); 
         qdot = qdot + Ainv*J.'*F_z;
         
         %Horizontal
        J_x = J2(1,:);
        lambda_x = 1/(J_x * Ainv * J_x.');
        F_x = lambda_x * (0 - J_x * qdot);
        qdot = qdot + Ainv*J_x.'*F_x;
    end
    
    %leg 2
    if (rToe2 < 0 && vToe2 < 0)
        J  = jacobian_toe2(z,p); 
        J2 = J; 
        J = J(2,:); 
        
      M = A_leg(z,p);
      Ainv = inv(M);
      
      lambda_z = 1/(J * Ainv * J.');
        F_z = lambda_z*(0 - J*qdot); 
        qdot = qdot + Ainv*J.'*F_z;
        
        %Horizontal
        J_x = J2(1,:);
        lambda_x = 1/(J_x * Ainv * J_x.');
        F_x = lambda_x * (0 - J_x * qdot);
        qdot = qdot + Ainv*J_x.'*F_x;
    end
        
    if(rHeel2 < 0 && vHeel2 < 0)
         J = jacobian_heel2(z,p);
         J2 = J; 
         J = J(2,:); 
         
          M = A_leg(z,p);
      Ainv = inv(M);
      
      lambda_z = 1/(J * Ainv * J.');
         
         F_z = lambda_z*(0 - J*qdot); 
         qdot = qdot + Ainv*J.'*F_z;
         
         %Horizontal
        J_x = J2(1,:);
        lambda_x = 1/(J_x * Ainv * J_x.');
        F_x = lambda_x * (0 - J_x * qdot);
        qdot = qdot + Ainv*J_x.'*F_x;
    end
    
        
  
    c = 1;
    
    
end

function tau = control_law(t,z,p)
    global ptraj; 
    % Controller gains, Update as necessary for Problem 1
    K_x = 400; % Spring stiffness X
    K_y = 400; % Spring stiffness Y
    D_x = .4;  % Damping X
    D_y = .4;  % Damping Y
    
    K_x2 = 300; %300
    K_y2 = 300; %300
    D_y2 = .17;
    D_x2 = .17;
  
   
   
    
    %% STEPS TO COMPLETE PROBLEM 1.3
    % a. Compute r_E
    
    %ONLY CONTROLLING ANKLE
   rE = position_ankle(z,p); %THIS VALUE NOT HAPPY
    % b. Compute J, the jacobian of r_E
   jE = jacobian_ankle(z,p); %THIS VALUE NOT HAPPY
   jE = jE(:, 3:end);
   vF = velocity_ankle(z,p);
   
   %leg 2
   rE2 = position_ankle2(z,p); %THIS VALUE NOT HAPPY
    % b. Compute J, the jacobian of r_E
   jE2 = jacobian_ankle2(z,p); %THIS VALUE NOT HAPPY
   jE2 = jE2(:, 3:end);
   vF2 = velocity_ankle2(z,p);
   
   % EllipsePositioning
   x1 = ptraj.rx*cos(ptraj.omega*t) + z(1) + ptraj.x0; 
   y1 = ptraj.ry*sin(ptraj.omega*t) + z(2) - ptraj.y0; 
   
   x2 = ptraj.rx*cos(ptraj.omega*t + ptraj.phase) + z(1) + ptraj.x0; 
   y2 = ptraj.ry*sin(ptraj.omega*t + ptraj.phase) + z(2) - ptraj.y0;
   
   rEd = [x1, y1];  %TRACK SINUSOID , make point g move in an ellipse
   rEd2 = [x2, y2];
   %denom = ptraj.y0 - y2;
   %thetaDes = 2*atan(ptraj.rx/ptraj.y0)*cos(ptraj.omega*t+ptraj.phase); 
   
   tau = -transpose(jE)*[K_x*(rE(1)-rEd(1))-D_x*vF(1);K_y*(rE(2)-rEd(2))-D_y*vF(2)];  %DAMPING
   tau = tau -transpose(jE2)*[K_x2*(rE2(1)-rEd2(1))-D_x2*vF2(1);K_y2*(rE2(2)-rEd2(2))-D_y2*vF2(2)];
end


function Fc = contact_force_toe(z,p)

    %% Fixed parameters for contact
    K_c = 500;
    D_c = 50;
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
    K_c = 500; %100000
    D_c = 50;
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

function Ankle_Res = ankle_constraint(z)
    th3 = z(5);
    dth3 = z(13); 
    th32 = z(8); 
    dth32 = z(16);
    
    if th3 <= -pi/6
        th3=-pi/6;
        dth3 = 0; 
    end
    if th3 > pi/8
        th3= pi/8;
        dth3 =  0;  
    end
    
    if th32 <= -pi/6
        th32=-pi/6;
        dth32 = 0; 
    end
    if th32 > pi/8
        th32= pi/8;
        dth32 =  0;  
    end
    
	Ankle_Res = [th3,th32, dth3, dth32];
end


function Tauc = joint_limit_torque(z,p)
    %% Fixed parameters for rotational spring damper at joint limit contact
    Kappa_c = 10;
    Dampa_c = 0.2;

    Tauc = 0;
end


function [dz tau] = dynamics(t,z,p,mu)
%     x = z(1);        y= z(2)      
%    th1 = z(3);     th2 = z(4);      th3 = z(5);
%     
%     dx = z(5);       dy = z(6);     
%     dth1 = z(7);  dth2= z(8);      dth3 = z(9);
    
    %check ankle
%     Ankle_Res = ankle_constraint(z);
%     z(5) = Ankle_Res(1);
%     z(11) = Ankle_Res(2);
    
    %y constraint
    %z([2,8]) = 0; 

    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p);
  
%     % Compute the contact force (used for problem 2)
%     Fc = contact_force_toe(z,p);
%     Fh = contact_force_heel(z,p);
%     
%     %if discrete
%     if 1
%         Fc = 0;
%         Fh = 0;
%     end
%    
%     J_toe = jacobian_toe(z,p);
%     J_heel = jacobian_heel(z,p);
%         
%     % Compute the contribution of the contact force to the generalied force
%     QFc_t=  transpose(J_toe)*[0;Fc];
%     QFc_h=  transpose(J_heel)*[0;Fh];
%     QFc = QFc_t +  QFc_h;   
   
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
    
    % Compute the contact force (used for problem 2.5)
    Tauc = joint_limit_torque(z,p);
    QTauc= [0; 0; 0; 0; 0 ; 0]; %column vector
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:8) = z(9:16);
    dz(9:16) = qdd;
    
    %y constraint again
    %dz([2,8]) = 0;
    
end

function animateSol(tspan, x,p)
    global ptraj
    % Prepare plot handles
    filename = 'testAnimated.gif';
    h=figure(3); 
    clf;
    
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_IG = plot([0],[0],'LineWidth',2);
    h_GH = plot([0],[0],'LineWidth',2);
    %THIS IS OUR ANKLE SPRING
    h_ankle = plot([0],[0],'LineWidth',0.5);

    ground  = plot([0],[0],'LineWidth',0.3);
    
    %ellipse
    ellip = plot([0], [0], '--'); 
    dotellip = plot([0], [0], 'o');
    j = plot([0], [0], 'or', 'MarkerSize', 3);
    
    
    %leg 2
    % Prepare plot handles
    hold on
    h_OB2 = plot([0],[0],'LineWidth',2);
    h_AC2 = plot([0],[0],'LineWidth',2);
    h_BD2 = plot([0],[0],'LineWidth',2);
    h_CE2 = plot([0],[0],'LineWidth',2);
    h_IG2 = plot([0],[0],'LineWidth',2);
    h_GH2 = plot([0],[0],'LineWidth',2);
    %THIS IS OUR ANKLE SPRING
    h_ankle2 = plot([0],[0],'LineWidth',0.5);
    
    

    ground  = plot([0],[0],'LineWidth',0.3);
    
    
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 1 -.2 .5]);
    
    pos = get(gcf, 'Position');
    width = pos(3);
    height = pos(4);

    mov = zeros(height, width, 1, length(tspan), 'uint8');
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
        
        r02 = keypoints(:,1+8);%don't use this
        rA2 = keypoints(:,2+8); % Vector to base of cart
        rB2 = keypoints(:,3+8);
        rC2 = keypoints(:,4+8); % Vector to tip of pendulum
        rD2 = keypoints(:,5+8);
        rE2 = keypoints(:,6+8);
        rI2 = keypoints(:,7+8);
        rH2 = keypoints(:,8+8);
        r_heela2 = keypoints(:,9+8);

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
%         set(h_ankle, 'XData' , [r_heela(1) rI(1)] );
%         set(h_ankle, 'YData' , [r_heela(2) rI(2)] );
        
        %leg2 
        set(h_OB2,'XData',[r0(1) rB2(1)]); %should be x rB
        set(h_OB2,'YData',[r0(2) rB2(2)]); %DO WE HAVE Y MOOVEMENT? HIP MIGHT NATURALLY RISE AND FALL
        
        set(h_AC2,'XData',[rA2(1) rC2(1)]);
        set(h_AC2,'YData',[rA2(2) rC2(2)]);
        
        set(h_BD2,'XData',[rB2(1) rD2(1)]);
        set(h_BD2,'YData',[rB2(2) rD2(2)]);
        
        set(h_CE2,'XData',[rC2(1) rE2(1)]);
        set(h_CE2,'YData',[rC2(2) rE2(2)]);
        
        %FOOT
         set(h_IG2,'XData',[rE2(1) rI2(1)]);
         set(h_IG2,'YData',[rE2(2) rI2(2)]);
        
        set(h_GH2,'XData',[rE2(1) rH2(1)]);
        set(h_GH2,'YData',[rE2(2) rH2(2)]);
        
        %ellipse trajectory
        x1 = ptraj.rx*cos(ptraj.omega*tspan) + z(1) + ptraj.x0; 
        y1 = ptraj.ry*sin(ptraj.omega*tspan) + z(2) - ptraj.y0;

        set(ellip, 'XData', x1); 
        set(ellip, 'YData', y1);
        
         x1 = ptraj.rx*cos(ptraj.omega*t) + z(1) + ptraj.x0; 
        y1 = ptraj.ry*sin(ptraj.omega*t) + z(2) - ptraj.y0;
        
        set(dotellip, 'XData', x1); 
        set(dotellip, 'YData', y1); 
        
        % where the robot thinks the foot is
%         [jankle] = jacobian_ankle(z,p); 
%         jankle = jankle*z(1:6); 
%         set(j, 'XData', jankle(1)); 
%         set(j, 'YData', jankle(2)); 
        
        %Ground
        set(ground, 'XData' , [-2 2] );
        set(ground, 'YData' , [0 0] );
        pause(.08)
        
        % Get frame as an image
       f = getframe(gcf);
       img =  frame2im(f);
       [img,cmap] = rgb2ind(img,256);

       % Create a colormap for the first frame. For the rest of the frames, use
       % the same colormap
       if i == 30
          imwrite(img,cmap,'k=0.gif','gif','LoopCount',Inf,'DelayTime',1);
       else
          imwrite(img,cmap,'k=0.gif','gif','WriteMode','append','DelayTime',1);
       end
    end
    imwrite(mov, cmap, 'k=0.gif', 'DelayTime', 0, 'LoopCount', inf)
end