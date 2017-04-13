function [] = varyPressure(t, state, parameters)
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
Last modified:	04/12/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extract necessary parameters and system
rho_water = parameters(3);
R = parameters(4);
T_air_i = parameters(7);
vol_water_i = parameters(12);
vol_bottle = parameters(13);
m_bottle = parameters(19);
p_0 = parameters(11);
p_ambient = parameters(6);

windvector = [0; 0; 0];

%% Vary initial air pressure

% Range of pressures to test
p_gage = linspace(32, 60, 500)*6895;%20-60 psi->Pa
pressure = p_ambient + p_gage;
w = waitbar(0,'Varying initial air preesure...'); % progress bar

x = zeros(1, length(pressure));
z = zeros(1, length(pressure));

% Loop through different initial volumes of water
for i = 5:length(pressure)

	% Change initial values in 'system' vector
    vol_air = vol_bottle - vol_water_i;
	m_air_0 = (pressure(i)/R/T_air_i)*(vol_air); %kg
	m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
	state(2) = m_R_0;
	state(3) = m_air_0;

	% Change initial values in 'parameters' vector
	parameters(11) = pressure(i);
	parameters(14) = m_air_0;

	% Solve system
	[~,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,...
	parameters, windvector, 45),t,state);

	% Find max height and distance achieved
	x(i) = max(dsdt(:,7));
	z(i) = max(dsdt(:,9));

	waitbar(i/length(pressure)); % update progress bar
end
close(w); % close progress bar

% Remove outliers
%a = [0 diff(x)];
%b = std(a)*ones(1,length(x));
%pressure(a >= 3*b) = [];
%z(a >= 3*b) = [];
%x(a >= 3*b) = [];

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(pressure/6895, x)
title('Variation of Initial Air Pressure')
xlabel('Inital Air Pressure (psi)')
ylabel('Distance Achieved (m)')

% Find and plot line of best fit
coefs = polyfit(pressure/6895, x, 3);
xFit = linspace(min(pressure/6895), max(pressure/6895), 1000);
yFit = polyval(coefs, xFit);
plot(xFit, yFit)

% Find optimal initial volume of water and display to command window
optVol = pressure(x == max(x)); 
disp(['The optimal initial air pressure is ' num2str(optVol) ...
	' liters.'])

end
