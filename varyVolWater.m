function [] = varyVolWater(t, state, parameters)
%varyVolWater Defines a range of initial volumes of water, calculates the 
% trajectory for each case, plots the findings, and displays the optimal 
% initial volume of water.
%{
The purpose of this program is to find the optimal initial volume of water
to achieve the maximum distance.

Inputs:	t	- time
	system	- system parameters
	parameters - environmental parameters

Outputs: ------------------------------------

Created by:	Keith Covington
Created on:	03/21/2017
Last modified:	04/13/2017
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

windvector = [0; 0; 0];

%% Vary initial volume of water to find optimum volume

% Range of volumes to test
volume = linspace(0.1*vol_bottle, 0.5*vol_bottle, 500);

w = waitbar(0,'Varying initial water volume...'); % progress bar

x = zeros(1, length(volume));
z = zeros(1, length(volume));

% Loop through different initial volumes of water
for i = 5:length(volume)

	% Change initial values in 'system' vector
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
	parameters, windvector, 45),t,state);

	% Find max height and distance achieved
	x(i) = max(dsdt(:,7));
	z(i) = max(dsdt(:,9));

	waitbar(i/length(volume)); % update progress bar
end
close(w); % close progress bar

volume = 1000.*volume; % convert water volume from [m^3] to [liters]

% Remove outliers
a = [0 diff(x)];
b = std(a)*ones(1,length(x));
volume(a >= 3*b) = [];
z(a >= 3*b) = [];
x(a >= 3*b) = [];

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(volume, x)
title('Variation of Initial Volume of Water')
xlabel('Inital Volume of Water (liters)')
ylabel('Distance Achieved (m)')

% Find and plot line of best fit
coefs = polyfit(volume, x, 7);
xFit = linspace(min(volume), max(volume), 1000);
yFit = polyval(coefs, xFit);
plot(xFit, yFit)

% Find optimal initial volume of water and display to command window
optVol = volume(x == max(x)); 
disp(['The optimal initial volume of water is ' num2str(optVol) ...
	' liters.'])

end
