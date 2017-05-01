function [] = monteCarlo(t, state, parameters, windvector, opts)
%monteCarlo plots the error ellipses from Monte Carlo simulation.
%{
The purpose of this program is to 

Inputs: t		- time
        system		- system parameters
        parameters	- environmental parameters
	windvector	- wind vector defined at wind.m call in Main.m
	opts		- ode45 options used to stop integration

Outputs: --------------------------------------

Created by:     Keith Covington
Created on:     04/25/2017
Last modified:  04/30/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_var = 100; % number of variations of each parameter
rng('shuffle')

% Define error in parameters
m_bottle_var = 0.05;	% error in dry mass of bottle
theta_var = 1;		% error in launch angle
cd_var = 0.0068;	% error in coefficient of drag
%cd_var = 0.068;	% error in coefficient of drag
wind_var = 1.34112;	% error in wind measurement
Isp_var = 0.0348;	% error in Isp calculation

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


wbar = waitbar(0, 'Running Monte Carlo simulation...'); % progress bar

% Run Monte Carlo simulation on thermo model
for var = 1:num_var

	% Generate random number
	%{
	randVar = randn;
	randVar = randVar - floor(randVar);
	if randVar < 0
		randVar = -randVar;
	end
	%}

	% Generate random launch parameters
	bottleVar = 0.001*((-m_bottle_var)+2*m_bottle_var*rand);
	theta = 45-theta_var + 2*theta_var*rand;
	cdVar = -cd_var + 2*cd_var*rand;
	windVarX = -wind_var + 2*wind_var*rand;
	windVarY = -wind_var + 2*wind_var*rand;
	Isp = 1.5-Isp_var + 2*Isp_var*rand;

	% Modify parameters
	state(2) = state(2)+bottleVar;
	parameters(9) = parameters(9) + cdVar;
	windvector = [windvector(1)+windVarX; windvector(2)+windVarY; windvector(3)];

	% Solve thermo system
	[t,allStates] = ode45(@(t,state) rocketTrajectory(t,state, ...
		parameters,windvector, theta),t,state,opts);
	plot3(allStates(:,7),allStates(:,8),allStates(:,9)); %plot(x,y,z)
	landings(var,:) = [allStates(end,7), allStates(end,8)]; % record landing

	waitbar(var/num_var); % update progress bar
end
close(wbar); % close progress bar

% Plot error elipses
plotErrorEllipse(landings(:,1), landings(:,2));

