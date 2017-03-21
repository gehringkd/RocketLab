function [] = varyVolWater(t, system, parameters)
% Does things and stuff
%{
The purpose of this program is to...
Inputs:	
Outputs: 
Assumptions: All processes during flight are adiabatic, and all outside
    forces (such as wind) are negligible. Limits of each variable are also
    assumed:
        - max. allowed gage pressure calculated in
        bottleRocketSetParameters
        - the drag coefficient can only vary between 0.3 and 0.5

Created by:	Keith Covington
Created on:	03/21/2017
Last modified:	03/21/2017
%}

%% Extract necessary parameters and system
%parameters = [g; gamma; rho_water; R; rho_air_ambient; p_ambient; ...
% T_air_i; C_d; C_D; p_gage; p_0; vol_water_i; vol_bottle; ... 
% m_air_0; v_0; A_throat; A_bottle; gage_p_allow; m_bottle];

% Parameters
g = parameters(1); %m/s^2
gamma = parameters(2); %specific heat ratio of air
rho_water = parameters(3);
R = parameters(4);
rho_a = parameters(5); %kg/m^3, density of ambient air
p_ambient = parameters(6);
T_air_i = parameters(7);
C_D = parameters(9);
p_gage = parameters(10);
p_0 = parameters(11);
vol_water_i = parameters(12);
vol_bottle = parameters(13);
m_air_i = parameters(14); %kg, initial mass of air
v_0 = parameters(15);
gage_p_allow = parameters(18);
m_bottle = parameters(19);

% System
%s = [v_0 m_R_0 m_air_0 theta_0 V_0 x_0 z_0];
vol_air = system(1);
m_R = system(2);
m_air = system(3);
theta = system(4);
V = system(5);
x = system(6);
z = system(7);


%% Get initial distance x and height z for comparison
[~,dsdt] = ode45(@(t,system) rocketTrajectory(t,system,...
parameters),t,system);
x = max(dsdt(:,6));
z = max(dsdt(:,7));


%% Vary initial volume of water to find optimum volume

% Range of volumes to test
volume = linspace(0.00001, 0.5*vol_bottle, 1000);


w = waitbar(0,'Varying initial water volume...'); % progress bar

% Loop through different initial volumes of water
for i = 5:length(volume)

	% Change initial values in 'system' vector
	vol_air = vol_bottle - volume(i);
	m_air_0 = (p_0/R/T_air_i)*(vol_air); %kg
	m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
	system(1) = vol_air;
	system(2) = m_R_0;
	system(3) = m_air_0;

	% Change initial values in 'parameters' vector
	parameters(12) = volume(i);
	parameters(14) = m_air_0;
	parameters(15) = vol_air;

	% Solve system
	[~,dsdt] = ode45(@(t,system) rocketTrajectory(t,system,...
	parameters),t,system);

	% Find max height and distance achieved
	x(i) = max(dsdt(:,6));
	z(i) = max(dsdt(:,7));

	waitbar(i/length(volume)); % update progress bar
end
close(w); % close progress bar

volume = 1000.*volume; % convert water volume from [m^3] to [liters]

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(volume, x)
title('Variation of Initial Volume of Water')
xlabel('Inital Volume of Water (liters)')
ylabel('Distance Achieved (m)')

optVol = volume(x == max(x)); % optimal initial volume of water...converted to [liters]
disp(['The optimal initial volume of water is ' num2str(optVol) ' liters.'])

end
