function [] = plotErrorEllipse(x, y)
%plotErrorEllipse plots the error ellipses from Monte Carlo simulation.
%{
The purpose of this program is to 

Inputs: landings - (x,y) coords of landings
        

Outputs: ------------------------------------

Created by:     Torin Clark
Created on:     04/25/2017
Last modified:  04/25/2017 by Keith Covington
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
plot(x,y,'k.','markersize',6)
axis equal; grid on;
title('Monte Carlo - Landing Sites');
xlabel('Downrange Distance [m]');
ylabel('Crossrange Distance [m]');
hold on;

% Calculate covariance matrix
P = cov(x,y);
mean_x = mean(x);
mean_y = mean(y);

% Calculate and define the error ellipses
n = 100; % number of points around ellipse
p = 0:pi/n:2*pi; % angles around a circle

[eigvec,eigval] = eig(P); % compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);

% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b')
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g')
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r')


%% Manually plot our actual landing coordinates

plot(64.3128,0, 'r*','markersize',50)



