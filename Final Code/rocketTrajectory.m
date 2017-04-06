function [dsdt] = rocketTrajectory(t,state,parameters)
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
%parameters = [g; gamma; rho_water; R; rho_air_ambient; p_ambient; ...
 %       T_air_i; C_d; C_D; p_gage; p_0; vol_water_i; vol_bottle; ... 
 %       m_air_0; v_0; A_throat; A_bottle; m_bottle];
g = parameters(1);	%m/s^2
gamma = parameters(2);	%specific heat ratio of air
rho_w = parameters(3);	%kg/m^3, density of water
R = parameters(4);	%kg/m^3, density of water; % J/kg/K
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
V_wind = parameters(end); % velocity of wind [m/s]


%% Make array s (inputs of rocket system) useable

v = state(1);	% volume of air in tank [m^3]
m_R = state(2);	% mass of rocket [kg]
m_air = state(3);	% mass of air in tank [kg]
th = state(4);	% theta - angle of rocket from trajectory [rad]
Vx = state(5);	% x-component of velocity [m/s]
Vy = state(6);	% y-component of velocity [m/s]
Vz = state(7);	% z-component of velocity [m/s]
x = state(8);	% x-position [m]
y = state(9);	% y-position [m]
z = state(10);	% z-position [m]


%% Get state parameters passed into function
% Define velocity vector and magnitude of velocity
V = [Vx; Vy; Vz]


%% Take-off Phase
%{

Maybe we should consider adding a phase for when the rocket is sliding past
the launch rails...

We could:
- Subtract the frictional force from the thrust force
- Add a normal force from the rails to make sure bottle doesn't spike into
  the ground.
- Define direction of V as a static value

%}


%% First Phase - Thrust generation before water is exhausted
% This phase continues until the volume of the air is equal to the volume of
% the bottle

if v < Vol_B

    disp('1st Stage');

    %during this phase, mass of the air remains constant, so pressure is
    %changing and can be found as a ratio to volume against intitial values
    p = p_0*((v_0/v)^gamma); %Pa
    
    %volume of air increases with time -> rate of change of vol. of air, dvdt:
        %v initially v_air_i, comes from dvdt otherwise
    dvdt = c_d*A_t*sqrt((2/rho_w)*(p-p_a)); %dvdt = m^3/s, v = m^3
    
    %The thrust, F, is a function of mass flow: F = m_dot*V_e
    F = 2*c_d*(p-p_a)*A_t;
    
    %and mass is decreasing with time according to:
    %dmdt = -c_d*A_t*sqrt(2*rho_w*(p-p_a)); %dmdt = kg/s, m = kg
    dm_airdt = 0;
    
    dm_Rdt = -c_d*A_t*sqrt(2*rho_w*(p-p_a));
    dm_Rdt2 = -c_d*A_t*rho_w*sqrt(2*(p-p_a)/rho_w);
    


%% Second Phase - Thrust generation after water is exhausted
%This phase continues until the difference of air pressure inside the
%bottle and the air pressure outside the bottle is 0 -> (p-p_a) = 0

elseif v >= Vol_B

    disp('2nd Stage');

    %p_end is the pressure of the air at the end of phase 1, which is when
    %the volume of the air is equal to the volume of the bottle
    p_end = p_0*(v_0/Vol_B)^gamma;
    
    %during this phase, volume of air remains constant but the pressure
    %continues to fall and can be calculated as a ratio of mass to pressure
    p = p_end*(m_air/m_air_i)^gamma;
    
if p > p_a
    dvdt = 0;
    rho = m_air/Vol_B;
    T = p/(rho*R);
    
    %Define critical pressure
    p_c = p*(2/(gamma+1))^(gamma/(gamma-1));

    if p_c > p_a %the flow is choked, exit Mach number M_e is 1
        T_e = 2*T/(gamma+1);
        p_e = p_c;
        M_e = 1;
    else %p_c <= p_a, the flow is not choked
        p_e = p_a;
        M_e = sqrt(((p/p_a)^((gamma-1)/gamma)-1)/((gamma-1)/2));
        T_e = T*(1+((gamma-1)/2)*(M_e^2));  
    end
    V_e = M_e*sqrt(gamma*R*T_e);
    rho_e = p_e/(R*T_e);
    
    dm_airdt = -c_d*rho_e*A_t*V_e;
    F = -dm_airdt*V_e+(p_e-p_a)*A_t;
    
    dm_Rdt = dm_airdt;

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
    dVdt = 0;
    dthdt = 0;
    dxdt = 0;
    dzdt = 0;

else


%{
How to adapt this to 2D version with wind consideration:

Vrel = Vrocket - V_wind
h = Vrel/|Vrel|   (vector)
T = Thrust*h
D = Drag*h
g = [0; -9.81]
%}

%{
After adapting this to the 2D version with wind, we can add in a y-component to the vectors to make it 3D.
%}

	% Displacement of rocket
	displacement = norm([x y z]);

	Vrel = V - V_wind % Vrel due to wind [m/s]
	Vrel_mag = norm(Vrel); % magnitude of relative velocity [m/s]
	headVec = Vrel/Vrel_mag; % heading vector (unit vector of Vrel)

	if displacement >= 1
		g = [0; 0; -9.81]; % g into vector form
	elseif displacement < 1
		g = [0; 0; 0];
		headVec = [sqrt(2); 0; sqrt(2)];
	else
		disp('Wuuuuuut');
	end

	%Equation for drag
	D = (rho_a/2)*(Vrel_mag^2)*c_D*A_B;  % [N]
	D = D.*headVec % drag vector
	F = F.*headVec % thrust vector

	%Given an initial velocity V, the change of V with respect to time is dVdt
	%below, with m being the mass of the rocket and th being theta
	%dVdt = (F-D-m_R*g*sin(th))/m_R; %dVdt = m/s^2, V=m/s

	sumF = F-D+g % vector of sum of forces [N]

	%derivatives (all vectors [x;y;z];
	dVdt = sumF./m_R % [m/s^2]


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
%dsdt(4) = dthdt;
dsdt(5) = dVdt(1);
dsdt(6) = dVdt(2);
dsdt(7) = dVdt(3);
dsdt(8) = dxdt;
dsdt(9) = dydt;
dsdt(10) = dzdt;

dsdt = dsdt';
end
