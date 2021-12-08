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
th = 0;
la = 0.01;
lb = 0.01;
lth = -pi/3;
ua = .1;
ub = .03;
uth = pi/3;


% % setup and solve nonlinear programming problem
 problem.objective = @(x) objective(x);     % create anonymous function that returns objective
 %problem.nonlcon = @(x) constraints(a,b);     % create anonymous function that returns nonlinear constraints
 problem.x0 = [ax by th];                   % initial guess for decision variables
 problem.lb = [la lb lth];     % lower bound on decision variables
 problem.ub = [ua ub uth];     % upper bound on decision variables
 problem.Aineq = []; problem.bineq = [];         % no linear inequality constraints
 problem.Aeq = []; problem.beq = [];             % no linear equality constraints
 problem.options = optimset('Display','iter');   % set options
 problem.solver = 'fmincon';                     % required
 x = fmincon(problem);                           % solve nonlinear programming problem
% Note that once you've solved the optimization problem, you'll need to 
% re-define tf, tfc, and ctrl here to reflect your solution.

a_optimized = x(1)
b_optimized = x(2)
th_optimized = x(3)
k = 0;

%% Plot COM for your submissions
MCoT = cost_of_transport_simulate_hip(a_optimized, b_optimized, th_optimized, k);
%%
% Run the animation
figure(3)                          
