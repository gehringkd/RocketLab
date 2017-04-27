function [dsdt] = rocketTrajectory(t,state,parameters,wind,angle,V1)
%BOTTLEROCKETTRAJECTORY is a system of equations meant to be solved using
%ode45 to plot the trajectory of a bottle rocket.
%{
The system of equations works through 3 phases of flight - 2 thrust
phases and the ballistic phase. 

Inputs: To run this function, use the input arrays t, s, and parameters,
    which are created and may be edited in bottleRocketSetParameters.m.

Outputs: This system should be run with ode45 to produce an output array
    dsdt of the changes in the system s with time, as well as an array of 
    times the differential equations were calculated for.

Assumptions: All calculations assume adiabatic processes and negligible
ooutside influences (such as wind).

Suggested function call:
[t,dsdt] = ode45(@(t,system) bottleRocketTrajectory(t,system,parameters) ...
    ,t,system);

Created by:	Kayla Gehring
Modified by:	Keith Covington
Created:	11/22/16
Modified:	04/04/17
%}


%% Extract needed parameters
rho_a = parameters(5);	%kg/m^3, density of ambient air
c_D = parameters(9);	%drag coefficent
A_B = parameters(17);	%m^2, area of bottle
V_wind = wind; % velocity of wind [m/s]
angle = angle*pi/180; %angle of launcher
m_R = parameters(19);

%% Make array s (inputs of rocket system) useable

Vx = state(1);	% x-component of velocity [m/s]
Vy = state(2);	% y-component of velocity [m/s]
Vz = state(3);	% z-component of velocity [m/s]
x = state(4);	% x-position [m]
y = state(5);	% y-position [m]
z = state(6);	% z-position [m]

V = [Vx; Vy; Vz];

	% Displacement of rocket
	displacement = norm([x y z]);

	Vrel = V - V_wind; % Vrel due to wind [m/s]
	Vrel_mag = norm(Vrel); % magnitude of relative velocity [m/s]
	headVec = Vrel/Vrel_mag; % heading vector (unit vector of Vrel)
%{
	if displacement >= 1
		g = [0; 0; -9.81]; % g into vector form
	elseif displacement < 1
		g = [0; 0; 0];
		headVec = [cos(angle); 0; sin(angle)];
        headVec = headVec/norm(headVec);
	else
		disp('Wuuuuuut');
	end
%}
    g = [0; 0; -9.81];
    
	%Equation for drag
	D = (rho_a/2)*(Vrel_mag^2)*c_D*A_B;  % [N]
	D = D.*headVec; % drag vector

	sumF = -D+g*m_R; % vector of sum of forces [N]

	%derivatives (all vectors [x;y;z];
	dVdt = sumF./m_R; % [m/s^2]

	% Components of velocity vector
	dxdt = V(1);
	dydt = V(2);
	dzdt = V(3);

%% Feed new system back out as a useable array
dsdt(1) = dVdt(1);
dsdt(2) = dVdt(2);
dsdt(3) = dVdt(3);
dsdt(4) = dxdt;
dsdt(5) = dydt;
dsdt(6) = dzdt;

dsdt = dsdt';
end
