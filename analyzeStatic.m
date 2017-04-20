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
%Group18 has a blip of data row 3515 - that's one of the offsets
%Other group is 23 or 24 maybe

% Trim it
thrust = thrust(startInd:endInd);

% Get time
time = (1/1.652/1000)*[1:length(thrust)]';	% correct time from 1.652 kHz to seconds

%% Find Peak thrust and total thrust time
peak = max(thrust);
delta_t = time(end) - time(1);

%% Calculate specific Impulse

mass = 1; %kg
Isp = calcImpulse(thrust, time, mass);


