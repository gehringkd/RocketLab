function [ dsdt] = rocketPhase3(t,s)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

global g c_D rho_a

%% Make array s (inputs of rocket system) useable
%s = [v,m,th,V,x,z]
s = s';
v = s(1);
m = s(2);
th = s(3);
V = s(4);
x = s(5);
z = s(6);
%% Third Phase - Ballistic
%This phase continues from the end of the second phase until the rocket
%hits the ground

    %There is no longer any thrust being produced and the mass is 
    %approximately equal to the mass of the bottle and is no longer
    %changing
    dvdt = 0;
    F = 0;
    dmdt = 0;  


%% General Rocket Trajectory
%These equations apply to all 3 phases

%The ground is z = 0 so:
if z <= 0
    dvdt = 0;
    dmdt = 0;
    dVdt = 0;
    dthdt = 0;
    dxdt = 0;
    dzdt = 0;
else
%Equation for drag
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

