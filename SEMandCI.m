function SEMandCI( data )
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
SEM = zeros(length(data),2);

for N = 1:length(data)
    SEM(N,:) = [std(data(1:N))/sqrt(N), N];
end

%Create Plot
figure
scatter(SEM(:,2), SEM(:,1));
xlabel('Number of Samples')
ylabel('\sigma_{I_{sp}}, s')
title('\sigma vs N')

%State the final SEM
disp(['Standard Error of the Mean, Isp: ', num2str(SEM(end,1))])
disp(['Number of data points: ', num2str(SEM(end,2))])

%% Calculate the 95%, 97.5%, and 99% confidence intervals
% produces mean +- z*SEM. z values for each CI:
CI = [95, 95, 97.5, 97.5, 99, 99];
z = [1.96, 1.96; 2.24, 2.24; 2.58, 2.58];
xbar = mean(data);
    CI = [xbar - z(:,1).*SEM(end,1), xbar + z(:,2).*SEM(end,1)];


% Create CI matrix
%Column 1: CI, Column 2: lower bound of CI, Column 3: upper bound of CI
Rows = {'95.0%';'97.5%';'99.0%'};
T = table(CI, 'RowNames', Rows);

% Export
writetable(T,'ConfidenceIntervals.csv','WriteRowNames', true);
display('ConfidenceIntervals.csv')
end

