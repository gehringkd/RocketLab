function Isp = calcImpulse(thrust, time,  mass)
%getISP calculates the specific impulse (Isp) of a propellant given the 
% force, time, and initial mass of the propellant
%
% Inputs: thrust and time data for integration, intitial mass of propellant


%% Account for recovery of dynamic load cell

% Make linear baseline for thrust curve
xfit = [time(1) time(end)];
yfit = [thrust(1) thrust(end)];
coefs = polyfit(xfit,yfit,1);
y = polyval(coefs, time);
plot(time, y)


%% Calculate impulse
% Integrate to find total impulse
thrust = thrust - y;
impulse = trapz(time,thrust);

Isp = impulse/(mass*9.81); %s

end
