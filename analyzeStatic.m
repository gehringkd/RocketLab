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
time = data(:,1);	% get time
thrust = data(:,3);	% get total recorded thrust
time = time.*(1/1.652*1000);	% correct time from 1.652 kHz to seconds


%%  Trim data

% Go ahead and get rid of some of the trailing data
thrust(end-1000:end) = [];

% Find beginning of useful data
startInd = find(thrust>=3,1);	% first value of thrust

% Find end of useful data
thrust = flip(thrust);		% flip the thrust vector around
endInd = find(thrust<=-3,1);	% last-ish value of thrust
endInd = length(thrust)-endInd;	% accound for flipped vector
thrust = flip(thrust);		% flip thrust vector back

% Trim it
thrust = thrust(startInd:endInd);


%% Find Peak thrust and total thrust time
peak = max(thrust);
delta_t = time(endInd) - time(startInd);
time = time(1:endInd-startInd+1);

%% Calculate specific Impulse
% moved this to a separate function, but haven't done much with it yet

%mass = 1 L * 1000 kg/L^3 or something
mass = 1;
Isp = calcImpulse(thrust, time, mass);

