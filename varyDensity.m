function [] = varyDensity(t, state, parameters)
%varyVolWater Defines a range of initial volumes of water, calculates the 
% trajectory for each case, plots the findings, and displays the optimal 
% initial volume of water.
%{
The purpose of this program is to find the optimal initial volume of water
to achieve the maximum distance.

Inputs: t       - time
        state   - system parameters
        parameters - environmental parameters

Outputs: ------------------------------------

Created by:     Keith Covington
Created on:     03/21/2017
Last modified:  04/13/2017 by Kayla Gehring
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extract necessary parameters and system
rho_water = parameters(3);
vol_water_i = parameters(12);
m_bottle = parameters(19);
m_air_0 = state(3);

windvector = [0;0;0];
%% Vary initial volume of water to find optimum volume

% Range of launch angles to test
rho_water = linspace(700, 1300, 100);

w = waitbar(0,'Varying fuel density...'); % progress bar

x = zeros(1, length(rho_water));
z = zeros(1, length(rho_water));

% Loop through different initial volumes of water
for i = 1:length(rho_water)

        % Change initial values in 'system' vector
        m_R_0 = m_bottle + rho_water(i)*vol_water_i + m_air_0; %kg
        parameters(3) = rho_water(i);
        state(2) = m_R_0;

        % Solve system
    	[~,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,...
        parameters, windvector, 45),t,state);

        % Find max height and distance achieved
        x(i) = max(dsdt(:,7));
        z(i) = max(dsdt(:,9));

        waitbar(i/length(rho_water)); % update progress bar
end
close(w); % close progress bar

% Remove outliers
%rho_water(x >= 3*std(x)) = [];
%z(x >= 3*std(x)) = [];
%x(x >= 3*std(x)) = [];

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(rho_water, x)
title('Variation of Liquid Density')
xlabel('Liquid Density (kg/m^3)')
ylabel('Distance Achieved (m)')

% Find and plot line of best fit
coefs = polyfit(rho_water, x, 2);
xFit = linspace(min(rho_water), max(rho_water), 1000);
yFit = polyval(coefs, xFit);
plot(xFit, yFit)


% Find optimal initial volume of water and display to command window
optVol = rho_water(x == max(x)); 
disp(['The optimal fuel density is ' num2str(optVol) ...
        ' kg/m^3.'])

end

