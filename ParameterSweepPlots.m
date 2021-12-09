%% Author: Laura Schwendeman
%date: 12/5/21
%purpose: cycle through different simulation parameters and make a plot

% chose your simulation parameters

%% vary K
K = logspace(1,2,20); 

rx = .01; 
ry = .03; 

theta = 0;

CostsOfTransport = zeros(1,length(K)); 

for i = 1:length(K)
    
    CostsOfTransport(i) = cost_of_transport_simulate_hip(rx, ry, theta, K(i)) 
    
end

figure; 
plot(K, CostsOfTransport); 
xlabel('Ankle Stiffness (N/m)'); 
ylabel('Cost of Transport'); 


%% Vary a and b
K = 3;
rx = .01:.005:.05; 
ry = .01:.005:.05;

[gridx, gridy] = meshgrid(rx, ry); 

costsOfTransport = zeros(1,numel(gridx)); 

for i = 1:numel(gridx)
    
    costsOfTransport(i) = cost_of_transport_simulate_hip(gridx(i), gridy(i), theta, K);
    
end

%% Plot the surface

z = zeros(size(gridx)); 
z(1:numel(gridx)) = costsOfTransport; 

figure; 
s = surf(gridx, gridy, z); 
xlabel('x-axis'); 
ylabel('y-axis'); 
zlabel('Cost of Transport'); 

s.EdgeColor = 'none';

% get cost of transport for the optimal solution and plot that too; 
hold on; 

rx = .01; 
ry = .03;

cost = cost_of_transport_simulate_hip(rx, ry, theta, K);

scatter3(rx, ry, cost, '.r');  

improvePlot(); 
