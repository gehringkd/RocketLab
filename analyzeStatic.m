function [stats, time, thrust] = analyzeStatic(file)
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
time = data(:,1);	% get time
thrust = data(:,3);	% get total recorded thrust


%%  Trim data

% Go ahead and get rid of some of the trailing data
thrust(end-1000:end) = [];

% Find beginning of useful data
startInd = find(thrust>=3,5);	% first value of thrust

% Find end of useful data
thrust = flip(thrust);		% flip the thrust vector around
endInd = find(thrust<=-3,1);	% last-ish value of thrust
endInd = length(thrust)-endInd;	% accound for flipped vector
thrust = flip(thrust);		% flip thrust vector back

% Trim it
thrust = thrust(startInd:endInd);
time = time(startInd:endInd);

%% Find Peak, Avg, and standard deviation of thrust values
stats = [max(thrust), mean(thrust), std(thrust)]; 

%% Calculate specific Impulse
% moved this to a separate function, but haven't done much with it yet

%mass = 1 L * 1000 kg/L^3 or something
Isp = calcImpulse(thrust, time, mass);

%% Return statement

