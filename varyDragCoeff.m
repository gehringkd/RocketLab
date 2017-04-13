function [] = varyDragCoeff(t, state, parameters)
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
Last modified:	03/21/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Extract necessary parameters and system
C_D = parameters(9);
windvector = [0 0 0]; %wind not needed for sensitivity analysis

%% Vary initial value of drag coefficient to find optimum volume

% Range of volumes to test
drag = linspace(0.05, .6, 1000);

w = waitbar(0,'Varying drag coefficient...'); % progress bar

x = zeros(1, length(drag));
z = zeros(1, length(drag));

% Loop through different values of drag
for i = 1:length(drag)

	% Change initial values in 'parameters' vector
	parameters(9) = drag(i);

	% Solve system
	[~,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,...
	parameters, windvector),t,state);

	% Find max height and distance achieved
	x(i) = max((dsdt(:,6));
	z(i) = max(dsdt(:,7));

	waitbar(i/length(drag)); % update progress bar
end
close(w); % close progress bar

% Remove outliers
drag(x >= 3*std(x)) = [];
z(x >= 3*std(x)) = [];
x(x >= 3*std(x)) = [];

% Plot various volumes and their respective achieved distance
figure
hold on
scatter(drag, x)
title('Variation of Drag Coefficient')
xlabel('C_D')
ylabel('Distance Achieved (m)')

% Find optimal initial volume of water and display to command window
optVol = drag(x == max(x)); 
disp(['The optimal initial volume of water is ' num2str(optVol) ...
	' liters.'])

end