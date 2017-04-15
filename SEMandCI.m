function [ CI ] = SEMandCI( data )
%SEM Creates a plot of standard error of the mean vs N samples of data and
%prints the final SEM for the full data set given. Also calculates
%confidence intervals(95%,97.5%,99%) of data.
%
% Inputs: An array of data that SEM Analysis is desired for
%
% Outputs: A graph of SEM vs N, printed statement of final SEM, an excel file
% of confidence interval ranges [CI percentage, lower bound, upper bound]
%
% Created by: Kayla Gehring
% Created on: 4/14
% Last Editted: 4/14

%% Create table of SEM vs N as N grows
SEM = zeros(1, length(data));

for N = 1:length(data)
    SEM(N, :) = [std(data(1:N))/sqrt(N), N];
end

%Create Plot
plot(SEM(2,:), SEM(1,:));

%State the final SEM
disp(['Standard Error of the Mean: ', num2str(SEM(1,end)])
disp(['Number of data points: ', num2str(SEM(2,end)])

%% Calculate the 95%, 97.5%, and 99% confidence intervals
% produces mean +- z*SEM. z values for each CI:
CI = [95,; 97.5; 99];
z = [1.96; 2.24; 2.58];
CI_low = xbar - z*SEM;
CI_high = xbar + z*SEM;

% Create CI matrix
%Column 1: CI, Column 2: lower bound of CI, Column 3: upper bound of CI
CI = [CI, CI_low, CI_high];

% Export to CSV for ease of putting in LaTeX
%I don't know the command 
end

