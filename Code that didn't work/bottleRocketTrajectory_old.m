function [dsdt] = bottleRocketTrajectory(t,s)
%bottleRocketTrajectory A system of equations meant to be solved using
%ode45 to plot the trajectory of a bottle rocket.
%   The system of equations works through 3 phases of flight - 2 thrust
%   phases and the ballistic phase. To run it, use the following input
%   arrays:
% t = [t0,tf], t0 is the initial time, tf is final time
% s = [v,m,th,V,x,z]
    %v = initial volume of air
    %m_R = initial mass of rocket
    %m_air = initial mass of air
    %th = theta, angle of rocket launch in radians
    %V = initial velocity, probably 0
    %x_0 = initial horizontal distance, probably 0
    %z_0 = initial height, make the ground 0! Calculations stop when z<0
%suggested function call: [dsdt] = ode45('bottleRocketTrajectory',t,s);
%The output will be an array of the values of s as they change with time.

%% Knowns that are unlikely to change 
gamma = 1.4; %specific heat ratio of air
rho_w = 1000; %kg/m^3, density of water
R = 287; % J/kg/K
g = 9.81; %m/s^2

%% Knowns likely to change - look here to change parameters of experiment
%Bottle parameters
d_t = .021; %m, diameter of throat of bottle
d_B = .105; %m, diameter of bottle
p_gage = 50*6894.76; %Pa, initial gage pressure of air in bottle
Vol_B = 0.002; %m^3, volume of the bottle
Vol_w_i = 0.001; %m^3, initial volume of water
m_B = 0.07; %kg, mass of empty bottle

%Flight parameters
%c_D = 0.5; %drag coefficent
c_d = 0.8; %discharge coefficient

%Atmosphere paramters
rho_a = 0.961; %kg/m^3, density of ambient air
p_a = 12.03*6894.76; %psi to Pa, ambient air pressure
T_air_i = 300; %K, initial temperature of air

%% Make array s (inputs of rocket system) useable
%s = [v,m,th,V,x,z]
s = s';
v = s(1);
m_R = s(2);
m_air = s(3);
th = s(4);
V = s(5);
x = s(6);
z = s(7);
p_0 = s(8);
c_D = s(9);

%% Calculations for constants
%valid if verification case data changes
v_0 = Vol_B - Vol_w_i;
%p_0 = p_gage + p_a;
A_t = pi*(d_t/2)^2; %m^2, area of throat
A_B = pi*(d_B/2)^2; %m^2, area of bottle 
%global p_end
%inital mass of entire bottle
    m_air_i = (p_0/(R*T_air_i))*v_0; %kg
    %m_water_i = 1000*Vol_w_i;
    %m_total_i = m_B+rho_w*(Vol_B-v_0)+m_air_i; %kg
    

%% First Phase - Thrust generation before water is exhausted
%This phase continues until the volume of the air is equal to the volume of
%the bottle

if v < Vol_B
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
    %dm_waterdt = -rho_w*c_d*A_t*sqrt(2*(p-p_a)/rho_w);
    
    %dm_Rdt = -c_d*A_t*sqrt(2*rho_w*(p-p_a));
    dm_Rdt = -c_d*A_t*rho_w*sqrt(2*(p-p_a)/rho_w);
    %dm_Rdt = dm_waterdt;

end
%% Constants and calculations for second phase

    %density of air inside bottle is p/(R*T), and mass is density*vol
    

%% Second Phase - Thrust generation after water is exhausted
%This phase continues until the difference of air pressure inside the
%bottle and the air pressure outside the bottle is 0 -> (p-p_a) = 0
 
%during this phase, volume of air remains constant but the pressure
    %continues to fall and can be calculated as a ratio of mass to pressure
        %mass of the air is the total mass minus the mass of the bottle
        %since no water is left
        %T_end = T_air_i*(v_0/Vol_B)^(gamma-1);
        p_end = p_0*(v_0/Vol_B)^gamma;
        p = p_end*(m_air/m_air_i)^gamma;

    
if p > p_a && v >= Vol_B
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
end

%% Third Phase - Ballistic
%This phase continues from the end of the second phase until the rocket
%hits the ground
if p <= p_a && v >= Vol_B && z > 0
    %There is no longer any thrust being produced and the mass is 
    %approximately equal to the mass of the bottle and is no longer
    %changing
    dvdt = 0;
    F = 0;
    dm_airdt = 0;
    %dm_waterdt = 0;
    dm_Rdt = 0;  
end

%% General Rocket Trajectory
%These equations apply to all 3 phases

%The ground is z = 0 so:
if z <= 0
    dm_airdt=0;
    %dm_waterdt=0;
    dvdt = 0;
    dm_Rdt = 0;
    dVdt = 0;
    dthdt = 0;
    dxdt = 0;
    dzdt = 0;
    dp_0dt = 0;
    dc_Ddt = 0;

else
%Equation for drag
D = (rho_a/2)*(V^2)*c_D*A_B; %N

%Given an initial velocity V, the change of V with respect to time is dVdt
%below, with m being the mass of the rocket and th being theta
dVdt = (F-D-m_R*g*sin(th))/m_R; %dVdt = m/s^2, V=m/s

%Given an initial angle theta, the change of theta with respect to time is:
if V < 1 %V == 0 causes an immediately negative trajectory
    dthdt = 0;
else
    dthdt = (-g*cos(th))/V;
end

%Given an initial horizontal distance x, the change of x with respect to
%time:
dxdt = V*cos(th);

%Given an initial height z, the change of z with respect to time:
dzdt = V*sin(th);

%Initial pressure and drag coefficient never change
dp_0dt = 0;
dc_Ddt = 0;

end

%% Feed new system back out as a useable array
%dsdt = [dvdt,dmdt,dthdt,dVdt,dxdt,dzdt]
dsdt(1) = dvdt;
dsdt(2) = dm_Rdt;
dsdt(3) = dm_airdt;
dsdt(4) = dthdt;
dsdt(5) = dVdt;
dsdt(6) = dxdt;
dsdt(7) = dzdt;
dsdt(8) = dp_0dt;
dsdt(9) = dc_Ddt;

dsdt = dsdt';
end
