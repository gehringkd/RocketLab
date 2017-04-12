function [] = varyVolWater(t,state,parameters,windvector)
%varyVolWater Defines a range of initial volumes of water, calculates the 
% trajectory for each case, plots the findings, and displays the optimal 
% initial volume of water.
%{
The purpose of this program is to find the optimal initial volume of water
to achieve the maximum distance.

Inputs:	t	- time
	state	- system parameters
	parameters - environmental parameters

Outputs: ------------------------------------

Created by:	Keith Covington
Created on:	03/21/2017
Last modified:	03/21/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extract necessary parameters and system
rho_water = parameters(3);
R = parameters(4);
T_air_i = parameters(7);
p_0 = parameters(11);
vol_water_i = parameters(12);
vol_bottle = parameters(13);
m_bottle = parameters(19);


%% Vary initial volume of water to find optimum volume

% Range of volumes to test
volume = linspace(0.00001, 0.5*vol_bottle, 1000);

w = waitbar(0,'Varying initial water volume...'); % progress bar

x = zeros(1, length(volume));
z = zeros(1, length(volume));

% Loop through different initial volumes of water
for i = 5:length(volume)

	% Change initial values in 'state' vector
	vol_air = vol_bottle - volume(i);
	m_air_0 = (p_0/R/T_air_i)*(vol_air); %kg
	m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
	state(1) = vol_air;
	state(2) = m_R_0;
	state(3) = m_air_0;

	% Change initial values in 'parameters' vector
	parameters(12) = volume(i);
	parameters(14) = m_air_0;
	parameters(15) = vol_air;

	% Solve system
	[~,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,...
	parameters,windvector),t,state);

	% Find max height and distance achieved
	x(i) = max(dsdt(:,7));
	z(i) = max(dsdt(:,9));

	waitbar(i/length(volume)); % update progress bar
end
close(w); % close progress bar

volume = 1000.*volume; % convert water volume from [m^3] to [liters]

% Remove outliers
volume(x >= 3.1*std(x)) = [];
z(x >= 3.1*std(x)) = [];
x(x >= 3.1*std(x)) = [];

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(volume, x)
title('Variation of Initial Volume of Water')
xlabel('Inital Volume of Water (liters)')
ylabel('Distance Achieved (m)')

% Find optimal initial volume of water and display to command window
optVol = volume(x == max(x)); 
disp(['The optimal initial volume of water is ' num2str(optVol) ...
	' liters.'])

end
