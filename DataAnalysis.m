function [ SEM ] = DataAnalysis(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Calculate Standard Deviation of data
sd = std(data);

%Calculate Mean of data
xbar = mean(data);

% divide by number of samples
SEM = sd/sqrt(length(data));

% Calculate the 95%, 97.5%, and 99% confidence intervals
% produces mean +- z*SEM. z values for each CI:
CI = [95,; 97.5; 99];
z = [1.96; 2.24; 2.58];
CI_low = xbar - z*SEM;
CI_high = xbar + z*SEM;

% Create CI matrix
%Column 1: CI, Column 2: lower bound of CI, Column 3: upper bound of CI
CI = [CI, CI_low, CI_high];

end

