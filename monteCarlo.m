function [] = monteCarlo(t, state, parameters, windvector, opts)
%monteCarlo plots the error ellipses from Monte Carlo simulation.
%{
The purpose of this program is to 

Inputs: t       - time
        system  - system parameters
        parameters - environmental parameters

Outputs: ------------------------------------

Created by:     Keith Covington
Created on:     04/25/2017
Last modified:  04/27/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_var = 10; % number of variations of each parameter

% Define error in parameters
m_bottle_var = 2;       % error in dry mass of bottle
theta_var = 1.5;        % error in launch angle
cd_var = 0.05;          % error in coefficient of drag
wind_var = 3;           % error in wind measurement
Isp_var = 0.1;          % error in Isp calculation

% Randomly vary parameters within realm of errors
m_bottle = 79-m_bottle_var + 2*m_bottle_var*rand(num_var,1);
theta = 45-theta_var + 2*theta_var*rand(num_var,1);
cd = 0.3-theta_var + 2*cd_var*rand(num_var,1);
%wind = 0-wind_var + 2*wind_var*rand(num_var,1);
Isp = 1.5-Isp_var + 2*Isp_var*rand(num_var,1);

% Initialize landing coordinates
landings = ones(length(num_var*4),2);	% landing coords in simulation

figure
hold on
grid on
axis equal
%ylim([-10 10]);
title('Rocket Flight - Monte Carlo Sim.')
xlabel('Downrange distance (m)')
ylabel('Crossrange distance (m)')
zlabel('Vertical height (m)')


bottleVar = -0.001*(m_bottle_var+2*m_bottle_var*rand(num_var,1))

% Vary dry mass of bottle
for p1 = 1:length(bottleVar)
	state(2) = state(2)+bottleVar(p1);
	%parameters(19) = m_bottle(p1);
	[t,allStates] = ode45(@(t,state) rocketTrajectory(t,state, ...
		parameters,windvector, 45),t,state,opts);
	plot3(allStates(:,7),allStates(:,8),allStates(:,9)); %plot(x,y,z)
	landings(p1,:) = [allStates(end,7), allStates(end,8)]; % record landing
end

% Vary dry mass of bottle
for p1 = 1:length(bottleVar)
	state(2) = state(2)+bottleVar(p1);
	[t,allStates] = ode45(@(t,state) rocketTrajectory(t,state, ...
		parameters,windvector, 45),t,state,opts);
	plot3(allStates(:,7),allStates(:,8),allStates(:,9)); %plot(x,y,z)
	landings(p1,:) = [allStates(end,7), allStates(end,8)]; % record landing
end




% Plot error elipses
%plotErrorEllipse();

