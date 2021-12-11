clear all; close all; clc;


% We can organize our code by filing things in different folders.  These
% folders need to be added to the Matlab path so that it can run the files
% inside them even when they are not the current folder listed at the top
% of the Matlab window.  For more information about the current folder, see
% http://www.mathworks.com/help/matlab/matlab_env/understanding-file-locations-in-matlab.html
% For more information about the Matlab path, see
% http://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html
%setpath                                     % add AutoDerived, Modeling, and Visualization folders to Matlab path


% Note: 5th state is the integral of torque squared over time
% An equation has been added to dynamics_continuous and dynamics_discrete
% to integrate this new state.

% Set guess


ax = 0.05;
by = 0.025;
th = pi/6;
la = 0.01;
lb = 0.01;
lth = -pi/4;
ua = .1;
ub = .03;
uth = pi/4;

res = zeros(3,3);
% A = [.01,.015,.02,.025,.03,.035,.04,.0.45,.05,.055,.06];
A = [30,50,100];

% % setup and solve nonlinear programming problem
for i = 1:length(A)
     problem.objective = @(x) objective(x,A(i));     % create anonymous function that returns objective
     %problem.nonlcon = @(x) constraints(a,b);     % create anonymous function that returns nonlinear constraints
     problem.x0 = [ax by];                   % initial guess for decision variables
     problem.lb = [la lb];     % lower bound on decision variables
     problem.ub = [ua ub];     % upper bound on decision variables
     problem.Aineq = []; problem.bineq = [];         % no linear inequality constraints
     problem.Aeq = []; problem.beq = [];             % no linear equality constraints
     problem.options = optimset('Display','iter');   % set options
     problem.solver = 'fmincon';                     % required
     x = fmincon(problem);                           % solve nonlinear programming problem
    % Note that once you've solved the optimization problem, you'll need to 
    % re-define tf, tfc, and ctrl here to reflect your solution.
    a_optimized = x(1)
    b_optimized = x(2)
    k_optimized = A(i)
    res(:,i) = [a_optimized,b_optimized,k_optimized];

    %% Plot COM for your submissions
    MCoT = simulate_system_leg_downloaded2(a_optimized, b_optimized, k_optimized);
end
res
% Run the animation
%figure(3)                          
