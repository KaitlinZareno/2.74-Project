function f = objective(x)
% Inputs:
%a axis of ellipse in x dir
%b axis of ellipse in y dir

% 
% Outputs:
% f - scalar value of the function (to be minimized) evaluated for the
%     provided values of the decision variables.
%
% Note: fmincon() requires a handle to an objective function that accepts 
% exactly one input, the decision variables 'x', and returns exactly one 
% output, the objective function value 'f'.  It is convenient for this 
% assignment to write an objective function which also accepts z0 and p 
% (because they will be needed to evaluate the objective function).  
% However, fmincon() will only pass in x; z0 and p will have to be
% provided using an anonymous function, just as we use anonymous
% functions with ode45().
   
    
    
     f = cost_of_transport_simulate_hip(x(1),x(2),x(3),0); 
     
end