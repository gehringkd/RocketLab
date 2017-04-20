function [Thrust_v_time] = group24Thrust()
%This function returns the useful portion of the thrust, and a vector of
%times, of group 24's thrust data.
%
% Inputs: none
% Outputs: a matrix with thrust in the first column and time in the second
%
% Created by Kayla Gehring, 4/20

    % Load our thrust data
    data = load('Group24_02PM_Statictest1');	% get data from file
    thrust = data(:,3).*4.44822;	% get total recorded thrust in N
    
    % Extract useful data   
    thrust = thrust(2460:2945);

    % Create matching time array
    timestep = 1/1.652/1000; %1.652 Hz to s
    time = timestep*[1:length(thrust)]';

    %Adjust zero-line to account for load cells
    xfit = [time(1) time(end)];
    yfit = [thrust(1) thrust(end)];
    coefs = polyfit(xfit,yfit,1);
    y = polyval(coefs, time);
    thrust = thrust-y;

    %Create output matrix
    Thrust_v_time = [thrust,time];
end
