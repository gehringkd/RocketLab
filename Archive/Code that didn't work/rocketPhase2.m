function [dsdt] = rocketPhase2(t,s)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

global gamma R g Vol_B m_B c_D c_d rho_a p_a T_air_i p_end

%% Make array s (inputs of rocket system) useable
%s = [v,m,th,V,x,z]
s = s';
v = s(1);
m = s(2);
th = s(3);
V = s(4);
x = s(5);
z = s(6);

%% Constants and calculations for second phase

    %density of air inside bottle is p/(R*T), and mass is density*vol
    m_air_i = (p_air_i/(R*T_air_i))*Vol_B;

%% Second Phase - Thrust generation after water is exhausted
%This phase continues until the difference of air pressure inside the
%bottle and the air pressure outside the bottle is 0 -> (p-p_a) = 0
 
%during this phase, volume of air remains constant but the pressure
    %continues to fall and can be calculated as a ratio of mass to pressure
        %mass of the air is the total mass minus the mass of the bottle
        %since no water is left
        m_air = m-m_B;
        p = p_end*(m_air/m_air_i)^gamma;
    
if p > p_a
    dvdt = 0;
    rho = m_air/Vol_B;
    T = p/(rho*R);
    
    %Define critical pressure
    p_c = p*(2/(gamma+1))^(gamma/(gamma-1));

if p_c > p_a %the flow is choked, exit Mach number M_e is 1
    T_e = (2/(gamma+1))*T;
    V_e = sqrt(gamma*R*T_e);
    p_e = p_c;
    rho_e = p_e/(R*T_e);
else %p_c <= p_a, the flow is not choked
    p_e = p_a;
    M_e = sqrt(((p/p_a)^((gamma-1)/gamma)-1)/((gamma-1)/2));
    T_e = T*(1+(gamma-1)/2*M_e^2);
    V_e = M_e*sqrt(gamma*R*T_e);
    rho_e = p_a/(R*T_e);
end
    dm_airdt = c_d*rho_e*A_t*V_e;
    F = dm_airdt*V_e+(p_e-p_a)*A_t;
    
    dmdt = -dm_airdt;
    
    %General trajectory    
D = rho_a/2*V^2*c_D*A_B;

%Given an initial velocity V, the change of V with respect to time is dVdt
%below, with m being the mass of the rocket and th being theta
dVdt = (F-D-m*g*sin(th))/m;

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
else 
    global t_end_phase2
    t_end_phase2 = t;
    dvdt = 0;
    dmdt = 0;
    dthdt = 0;
    dVdt = 0;
    dxdt = 0;
    dzdt = 0;
end

%% Feed new system back out as a useable array
%dsdt = [dvdt,dmdt,dthdt,dVdt,dxdt,dzdt]
dsdt(1) = dvdt;
dsdt(2) = dmdt;
dsdt(3) = dthdt;
dsdt(4) = dVdt;
dsdt(5) = dxdt;
dsdt(6) = dzdt;

dsdt = dsdt';

end

