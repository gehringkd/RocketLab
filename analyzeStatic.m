function [thrust] = analyzeStatic(file)
%analyzeStatic reads in and analyzes a file from static test stand tests.
%{
The purpose of this function is to analyze a single test stand file.

Inputs:	file	- name of file to be analyzed

Outputs: thrust	- vector of all total thrust values

Created by:     Keith Covington
Created on:     04/13/2017
Last modified:  04/13/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('Static Test Stand Data');	% add relative path to directory

data = load(file);			% get data from file
thrust(:,1) = data(:,3);		% get total recorded thrust

% Trim data
%find(thrust==0)



