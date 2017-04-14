function [] = plotTrajectory(fig,t,dsdt)
%plotTrajectory plots the 3D trajectory of the rocket.
%{
The purpose of this function is to properly format the 3D plotting of 
rocket trajectories.

Inputs: fig	- the figure object in which to plot
	t	- time values
	dsdt	- solutions to integration

Outputs: -------------------------------------------

Created by:     Keith Covington
Created on:     04/12/2017
Last modified:  04/12/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(fig); % Find figure

% Format and label plot
grid on
view(-30,10)
title('Verification Case - Bottle Rocket Flight')
xlabel('Downrange distance (m)')
ylabel('Crossrange distance (m)')
zlabel('Vertical height (m)')

% Plot (x,y,z)
plot3(dsdt(:,7),dsdt(:,8),dsdt(:,9));

end
