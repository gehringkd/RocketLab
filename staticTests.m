function [t, thrust] = staticTests()
%staticTests analyzes all static test stand data and produces data of 
% thrust vs. time.
%{
The purpose of this program analyize the test stand data from all groups 
 in order to produce reliable average static test stand data.

Inputs:	-------------------------------------

Outputs: t	- time vector
	 thrust	- vector of thrust values

Created by:     Keith Covington
Created on:     04/13/2017
Last modified:  04/14/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

% Get names of files to analyze
addpath('Static Test Stand Data');

direc = dir('Static Test Stand Data/G*');	% get directory where files are
numFiles = length(direc);		% number of files
files = cell(numFiles, 1);		% cell array for filenames
thrust = cell(numFiles, 1);		% cell array for thrust data

% Create data table for peak Thrust, Isp, and total time of each data set
thrust = zeros(3, length(direc)); %for Results, Q3
Isp = zeros(1, length(direc)); %for Results, Q2
timeT = zeros(1, length(direc)); %for Results, Q4

% Try plotting
figure
hold on

w = waitbar(0,'Analyzing static test stand data...'); % progress bar
% Loop though files
for i = 1:numFiles
	files{i} = direc(i).name;

	% remove files we don't want
	if ~contains(files{i}, 'Group')
		files{i} = [];		% delete cell/filename entry
	end

	% Analyze files
	%try
        %note from Kayla: analyze Static should return Peak, Avg, and
        %standard deviation of thrust for each data set. You can graph them
        %all if you'd like, but I think we need only one graph.
		thrust{i} = analyzeStatic(files{i});	% analyize file
	%catch
		%disp(['Could not read file: ' files{i} '. Skipping...']);
	%end

	plot(thrust{i})
	waitbar(i/numFiles); % update progress bar
end

close(w); % close progress bar




% dummy return statement for now
t = 0;
thrust = 0;
