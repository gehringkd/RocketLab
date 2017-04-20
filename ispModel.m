function [dsdt] = ispModel(t,state,parameters,wind,angle,thrustData)
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
[t,dsdt] = ode45(@(t,system) ispModel(t,system,parameters) ...
    ,t,system);

Created by:	Kayla Gehring
Modified by:	Keith Covington
Created:	04/20/17
Modified:	04/20/17
%}


%% Extract needed parameters
%parameters = [g; gamma; rho_water; R; rho_air_ambient; p_ambient; ...
 %       T_air_i; C_d; C_D; p_gage; p_0; vol_water_i; vol_bottle; ... 
 %       m_air_0; v_0; A_throat; A_bottle; m_bottle];
g = parameters(1);	%m/s^2
gamma = parameters(2);	%specific heat ratio of air
rho_w = parameters(3);	%kg/m^3, density of water
R = parameters(4);      % J/kg/K
rho_a = parameters(5);	%kg/m^3, density of ambient air
p_a = parameters(6);	%Pa
c_d = parameters(8);	%discharge coefficient
c_D = parameters(9);	%drag coefficent
p_0 = parameters(11);	%Pa
Vol_B = parameters(13);	%m^3, volume of the bottle
m_air_i = parameters(14); %kg, initial mass of air
v_0 = parameters(15);	%m^3, initial volume of air
A_t = parameters(16);	%m^2, area of throat
A_B = parameters(17);	%m^2, area of bottle
V_wind = wind; % velocity of wind [m/s]
angle = angle*pi/180; %angle of launcher


%% Make array s (inputs of rocket system) useable

v = state(1);	% volume of air in tank [m^3]
m_R = state(2);	% mass of rocket [kg]
m_air = state(3);	% mass of air in tank [kg]
Vx = state(4);	% x-component of velocity [m/s]
Vy = state(5);	% y-component of velocity [m/s]
Vz = state(6);	% z-component of velocity [m/s]
x = state(7);	% x-position [m]
y = state(8);	% y-position [m]
z = state(9);	% z-position [m]


%% Get state parameters passed into function
% Define velocity vector and magnitude of velocity
V = [Vx; Vy; Vz];

time = thrustData(:,1);
thrust = thrustData(:,2);


%% Propulsion Phase

if t <= time(end)
	approxInd = find(t>=time,1);
	F = thrust(approxInd);


%% Third Phase - Ballistic
%This phase continues from the end of the second phase until the rocket
%hits the ground
else
    %There is no longer any thrust being produced and the mass is 
    %approximately equal to the mass of the bottle and is no longer
    %changing
    dvdt = 0;
    F = 0;
    dm_airdt = 0;
    dm_Rdt = 0;  
end
end
    if v == Vol_B && p <= p_a
        error('there is something wrong with pressure calculations.')
    end
%% General Rocket Trajectory
%These equations apply to all 3 phases

%The ground is z = 0 so no calculations needed:
if z <= 0
    dm_airdt=0;
    dvdt = 0;
    dm_Rdt = 0;
    dVdt = [0 0 0];
    dxdt = 0;
    dydt = 0;
    dzdt = 0;

else

	% Displacement of rocket
	displacement = norm([x y z]);

	Vrel = V - V_wind;		% Vrel due to wind [m/s]
	Vrel_mag = norm(Vrel);		% magnitude of relative velocity [m/s]
	headVec = Vrel/Vrel_mag;	% heading vector (unit vector of Vrel)

	if displacement >= 1
		g = [0; 0; -9.81]; % g into vector form
	elseif displacement < 1
		g = [0; 0; 0];
		headVec = [cos(angle); 0; sin(angle)];
	else
		disp('Wuuuuuut');
	end

	%Equation for drag
	D = (rho_a/2)*(Vrel_mag^2)*c_D*A_B;  % [N]
	D = D.*headVec; % drag vector
	F = F.*headVec; % thrust vector

	%Given an initial velocity V, the change of V with respect to time is dVdt
	%below, with m being the mass of the rocket and th being theta
	%dVdt = (F-D-m_R*g*sin(th))/m_R; %dVdt = m/s^2, V=m/s

	sumF = F-D+g*m_R; % vector of sum of forces [N]

	%derivatives (all vectors [x;y;z];
	dVdt = sumF./m_R; % [m/s^2]


	%Given an initial angle theta, the change of theta with respect to time is:
	%if Vrel_mag < 1 %V == 0 causes an immediately negative trajectory
	%	dthdt = 0;
	%else
	%	dthdt = (-g(3)*cos(th))/Vrel_mag;
	%end

	%Given an initial horizontal distance x, the change of x with respect to
	%time:
	%dxdt = V*cos(th);

	% Components of velocity vector
	dxdt = V(1);
	dydt = V(2);
	dzdt = V(3);

	%Given an initial height z, the change of z with respect to time:
	%dzdt = V*sin(th);
	    
end


%% Feed new system back out as a useable array
dsdt(1) = dvdt;
dsdt(2) = dm_Rdt;
dsdt(3) = dm_airdt;
dsdt(4) = dVdt(1);
dsdt(5) = dVdt(2);
dsdt(6) = dVdt(3);
dsdt(7) = dxdt;
dsdt(8) = dydt;
dsdt(9) = dzdt;

dsdt = dsdt';
end
