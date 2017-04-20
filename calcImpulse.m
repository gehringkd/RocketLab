function Isp = calcImpulse(thrust, time,  mass)
%getISP calculates the specific impulse (Isp) of a propellant given the 
% force, time, and initial mass of the propellant
%
% Inputs: thrust and time data for integration, intitial mass of propellant
%% Calculate total impulse
% Total impulse will be the area under the thrust curve


% Account for recovery of dynamic load cell


% Integrate to find total impulse
impulse = trapz(time,thrust);

Isp = impulse/(mass*9.81); %s

end
