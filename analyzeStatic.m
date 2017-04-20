function [thrust, peak, delta_t, Isp] = analyzeStatic(file)
%analyzeStatic reads in and analyzes a file from static test stand tests.
%{
The purpose of this function is to analyze a single test stand file.

Inputs:	file	- name of file to be analyzed

Outputs: thrust	- vector of all total thrust values

Created by:     Keith Covington
Created on:     04/13/2017
Last modified:  04/14/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('Static Test Stand Data');	% add relative path to directory

%% Load data
data = load(file);	% get data from file
thrust = data(:,3);	% get total recorded thrust

% Define time step
timestep = 1/1.652/1000; %1.652 Hz to s

%%  Trim data
% Based off peak data to avoid those annoying offsets in some data sets
peak = max(thrust);
peakInd = find(thrust==peak,1);
startInd = peakInd - round(.03/timestep); %determined by lowering time period until all data appeared to be ~ correct
endInd = peakInd + round(.3/timestep); %make sure there's plenty of room after thrust
thrust = thrust(startInd:endInd);

% Find actual beginning of useful data
startInd = find(thrust>=3,1);	% first value of thrust

% Find end of useful data
thrust = flip(thrust);		% flip the thrust vector around
endInd = find(thrust<=-3,1);	% last-ish value of thrust
endInd = length(thrust)-endInd;	% accound for flipped vector
thrust = flip(thrust);		% flip thrust vector back
%Group18 has a blip of data row 3515 - that's one of the offsets
%Other group is 23 or 24 maybe

% Trim thrust
thrust = thrust(startInd:endInd);

%% Create time array and get total time
time = timestep*[1:length(thrust)]';
delta_t = time(end);


%% Calculate specific Impulse
mass = 1; %kg of water
Isp = calcImpulse(thrust, time, mass);


