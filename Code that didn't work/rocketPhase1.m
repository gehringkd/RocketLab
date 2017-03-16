function [dsdt] = rocketPhase1(t,s)

global gamma rho_w g p_gage Vol_B c_D c_d rho_a p_a v_0 A_t A_B

%% Make array s (inputs of rocket system) useable
%s = [v,m,th,V,x,z]
s = s';
v = s(1);
m = s(2);
th = s(3);
V = s(4);
x = s(5);
z = s(6);

%% Calculations for constants
%valid if verification case data changes

p_air_i = p_gage + p_a;

%% First Phase - Thrust generation before water is exhausted
%This phase continues until the volume of the air is equal to the volume of
%the bottle

if v < Vol_B
    %during this phase, mass of the air remains constant, so pressure is
    %changing and can be found as a ratio to volume against intitial values
    p = p_air_i*(v_0/v)^gamma; %Pa

    %volume of air increases with time -> rate of change of vol. of air, dvdt:
        %v initially v_air_i, comes from dvdt otherwise
    dvdt = c_d*A_t*sqrt(2/rho_w*(p-p_a)); %dvdt = m^3/s, v = m^3

    %The thrust, F, is a function of mass flow: F = m_dot*V_e
    F = 2*c_d*(p-p_a)*A_t; %Pa*m^2 = N
    
    %and mass is decreasing with time according to:
    dmdt = -c_d*A_t*sqrt(2*rho_w*(p-p_a)); %dmdt = kg/s, m = kg
    %Equation for drag

%General trajectory    
D = (rho_a/2)*V^2*c_D*A_B;

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
    global p_end
    p_end = p_air_i*(v_0/v)^gamma; %Pa
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

