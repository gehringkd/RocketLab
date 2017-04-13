function [] = varyAngle(t, state, parameters)
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
Last modified:  04/12/2017 by Kayla Gehring
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

windVector = [0;0;0];

%% Vary initial volume of water to find optimum volume

% Range of launch angles to test
angle = linspace(35, 55, 100);

w = waitbar(0,'Varying launch angle...'); % progress bar

x = zeros(1, length(angle));
z = zeros(1, length(angle));

% Loop through different initial volumes of water
for i = 1:length(angle)


        % Solve system
        [~,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,...
        parameters, windVector, angle(i)),t,state);

        % Find max height and distance achieved
        x(i) = max(dsdt(:,7));
        z(i) = max(dsdt(:,9));

        waitbar(i/length(angle)); % update progress bar
end
close(w); % close progress bar

% Remove outliers
%angle(x >= 3*std(x)) = [];
%z(x >= 3*std(x)) = [];
%x(x >= 3*std(x)) = [];

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(angle, x)
title('Variation of Launch Angle')
xlabel('Launch Angle (degrees)')
ylabel('Distance Achieved (m)')

% Find optimal initial volume of water and display to command window
optVol = angle(x == max(x)); 
disp(['The optimal initial angle is ' num2str(optVol) ...
        ' degrees.'])

end

