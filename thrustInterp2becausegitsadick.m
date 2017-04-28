function [dsdt] = thrustInterp(t,state,parameters,wind,angle,thrustData)
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

m_R = state(1)	% mass of rocket [kg]
Vx = state(2);	% x-component of velocity [m/s]
Vy = state(3);	% y-component of velocity [m/s]
Vz = state(4);	% z-component of velocity [m/s]
x = state(5);	% x-position [m]
y = state(6);	% y-position [m]
z = state(7);	% z-position [m]


%% Get state parameters passed into function
% Define velocity vector and magnitude of velocity
V = [Vx; Vy; Vz];

thrust = thrustData(:,1);
time = thrustData(:,2);


%% Propulsion Phase

if t <= time(end) && m_R >= 0
	if t<time(1,1)
        x1 = 0;
        x2 = time(1,1);
        y1 = 0;
        y2 = thrust(1,1);
    else
        ind1 = find(time>=t,1)-1;
        ind2 = ind1+1;

        x1 = time(ind1);
        x2 = time(ind2);
        y1 = thrust(ind1);
        y2 = thrust(ind2);
    end

	F = y1 + (t-x1)*(y2-y1)/(x2-x1);
<<<<<<< HEAD
	F = abs(F);
    
    if m_R >= .144+m_air_i %mass bottle + air, just troubleshooting right now
        V_e = sqrt(F/rho_w/A_t);
        dm_Rdt = -rho_w*A_t*V_e;
    elseif m_R >=.144
        m_air = m_R - m_air_i;
        rho_a = m_air/Vol_B;
        V_e = sqrt(F/rho_a/A_t);
        dm_Rdt = -rho_a*A_t*V_e;
    else
        dm_Rdt = 0;
        disp('hm');
    end
            
=======
	V_e = sqrt(F/rho_w/A_t);
	dm_Rdt = rho_w*A_t*V_e;
>>>>>>> ebdcd1c52de83e9540c8a61398096ed3538558a5

disp(['Phase 1  t=' num2str(t) ' s']);

%{
dsdt(1) = dvdt;
dsdt(2) = dm_Rdt;
dsdt(3) = dm_airdt;
dsdt(4) = dVdt(1);
dsdt(5) = dVdt(2);
dsdt(6) = dVdt(3);
dsdt(7) = dxdt;
dsdt(8) = dydt;
dsdt(9) = dzdt;
%}


%% Third Phase - Ballistic
%This phase continues from the end of the second phase until the rocket
%hits the ground
else
    %There is no longer any thrust being produced and the mass is 
    %approximately equal to the mass of the bottle and is no longer
    %changing
    F = 0;
    dm_Rdt = 0;  

disp(['Phase 2  t=' num2str(t) ' s']);
end



%% General Rocket Trajectory
%These equations apply to all 3 phases

%The ground is z = 0 so no calculations needed:
if z <= 0

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

	% Components of velocity vector
	dxdt = V(1);
	dydt = V(2);
	dzdt = V(3);

end


%% Feed new system back out as a useable array
dsdt(1) = dm_Rdt;
dsdt(2) = dVdt(1);
dsdt(3) = dVdt(2);
dsdt(4) = dVdt(3);
dsdt(5) = dxdt;
dsdt(6) = dydt;
dsdt(7) = dzdt;

dsdt = dsdt';
end
