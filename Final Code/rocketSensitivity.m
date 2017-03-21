function [] = bottleRocketSensitivity(t, system, parameters)
%UNTITLED Tests variation of bottle rocket flight as a function of initial
%angle, initial gage pressure inside bottle, initial volume of water, and
%the drag coefficient
%{
The purpose of this program is to test (individually) the effect that each
of 4 parameters has on the distance x and height z of a bottle rocket 
launch. It does this by varying the initial condition and running that
iinitial condition through bottleRocketTrajectory and ode45.

Inputs: This program takes in the arrays t, system, and parameters. See
    bottleRocketSetParameters.m for more detail on what these contain or to 
    change values.
Outputs: This program outputs 3 tables printed to the command window.
Assumptions: All processes during flight are adiabatic, and all outside
    forces (such as wind) are negligible. Limits of each variable are also
    assumed:
        - max. allowed gage pressure calculated in
        bottleRocketSetParameters
        - the drag coefficient can only vary between 0.3 and 0.5

Created by: 4250df2b7e8b
Created on: 11/23/16
Last modified: 12/2/16
%}

%% Extract necessary parameters and system
%parameters = [g; gamma; rho_water; R; rho_air_ambient; p_ambient; ...
% T_air_i; C_d; C_D; p_gage; p_0; vol_water_i; vol_bottle; ... 
% m_air_0; v_0; A_throat; A_bottle; gage_p_allow; m_bottle];

R = parameters(4);
T_air_i = parameters(7);
rho_water = parameters(3);
vol_water_i = parameters(12);
v_0 = parameters(15);
C_D = parameters(9);
p_ambient = parameters(6);
p_gage = parameters(10);
p_0 = parameters(11);
gage_p_allow = parameters(18);
vol_bottle = parameters(13);
m_bottle = parameters(19);

theta = system(4);

%% Get initial distance x and height z for comparison
[~,dsdt] = ode45(@(t,system) rocketTrajectory(t,system,...
parameters),t,system);
x = max(dsdt(:,6));
z = max(dsdt(:,7));

%% Test variation of flight as a function of the drag coefficient, C_D
%Vary C_D by 10% and 15%
C_D_10 = C_D*.9;
C_D_15 = C_D*.85;
C_D_25 = C_D*.75;

%make a modifiable parameters array
parameters_temp = parameters;

%Find final distance x and tallest height z for each C_D
parameters_temp(9) = C_D_10;
[~,dsdt] = ode45(@(t,system) rocketTrajectory(t,system,...
parameters_temp),t,system);
C_D_10_x = max(dsdt(:,6));
C_D_10_z = max(dsdt(:,7));
parameters_temp(9) = C_D_15;
[~,dsdt] = ode45(@(t,system) rocketTrajectory(t,system,...
parameters_temp),t,system);
C_D_15_x = max(dsdt(:,6));
C_D_15_z = max(dsdt(:,7));
parameters_temp(9) = C_D_25;
[~,dsdt] = ode45(@(t,system) rocketTrajectory(t,system,...
parameters_temp),t,system);
C_D_25_x = max(dsdt(:,6));
C_D_25_z = max(dsdt(:,7));

%Compare to original x and z values
C_D_10_x = (C_D_10_x-x)/x*100;
C_D_10_z = (C_D_10_z-z)/z*100;
C_D_15_x = (C_D_15_x-x)/x*100;
C_D_15_z = (C_D_15_z-z)/z*100;
C_D_25_x = (C_D_25_x-x)/x*100;
C_D_25_z = (C_D_25_z-z)/z*100;

%% Test variation of flight as a function of the initial angle, theta_0
%Vary theta by 10% and 15%
th_10 = theta*.9;
th_15 = theta*.85;
th_25 = theta*.75;

%make a modifiable system array
system_temp = system;

%Find final distance x and tallest height z for each C_D
system_temp(4) = th_10;
[~,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,system_temp,...
parameters),t,system_temp);
th_10_x = max(dsdt(:,6));
th_10_z = max(dsdt(:,7));
system_temp(4) = th_15;
[~,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,system_temp,...
parameters),t,system_temp);
th_15_x = max(dsdt(:,6));
th_15_z = max(dsdt(:,7));
system_temp(4) = th_25;
[~,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,system_temp,...
parameters),t,system_temp);
th_25_x = max(dsdt(:,6));
th_25_z = max(dsdt(:,7));

%Compare to original x and z values
th_10_x = (th_10_x-x)/x*100;
th_10_z = (th_10_z-z)/z*100;
th_15_x = (th_15_x-x)/x*100;
th_15_z = (th_15_z-z)/z*100;
th_25_x = (th_25_x-x)/x*100;
th_25_z = (th_25_z-z)/z*100;

%% %% Test variation of flight as a function of the pressure, p_0
%Vary theta by 10% and 15%
p0_10 = p_0*.9;
p0_15 = p_0*.85;
p0_25 = p_0*.75;

%reset modifiable parameters and system arrays
parameters_temp = parameters;
system_temp = system;

%Find final distance x and tallest height z for each C_D
parameters_temp(11) = p0_10;
    m_air_0 = (p0_10/R/T_air_i)*v_0; %kg
    m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
system_temp(1,2) = m_R_0;
system_temp(1,3) = m_air_0;
parameters_temp(14) = m_air_0;
[~,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,...
system_temp,parameters_temp),t,system_temp);
p0_10_x = max(dsdt(:,6));
p0_10_z = max(dsdt(:,7));
parameters_temp(11) = p0_15;
    m_air_0 = (p0_15/R/T_air_i)*v_0; %kg
    m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
system_temp(1,2) = m_R_0;
system_temp(1,3) = m_air_0;
parameters_temp(14) = m_air_0;
[~,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,...
system_temp,parameters_temp),t,system_temp);
p0_15_x = max(dsdt(:,6));
p0_15_z = max(dsdt(:,7));
parameters_temp(11) = p0_25;
    m_air_0 = (p0_25/R/T_air_i)*v_0; %kg
    m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
system_temp(1,2) = m_R_0;
system_temp(1,3) = m_air_0;
parameters_temp(14) = m_air_0;
[~,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,...
system_temp,parameters_temp),t,system_temp);
p0_25_x = max(dsdt(:,6));
p0_25_z = max(dsdt(:,7));


%Compare to original x and z values
p0_10_x = (p0_10_x-x)/x*100;
p0_10_z = (p0_10_z-z)/z*100;
p0_15_x = (p0_15_x-x)/x*100;
p0_15_z = (p0_15_z-z)/z*100;
p0_25_x = (p0_25_x-x)/x*100;
p0_25_z = (p0_25_z-z)/z*100;
%% Test variation of flight as a function of the initial volume of water
%Vary volume by 10% and 20%
v0_w = vol_bottle-v_0;
v0_10 = vol_bottle-v0_w*.9;
v0_15 = vol_bottle-v0_w*.85;
v0_25 = vol_bottle-v0_w*.75;

%reset modifiable parameters and system arrays
parameters_temp = parameters;
system_temp = system;

%Find final distance x and tallest height z for each C_D
parameters_temp(15) = v0_10;
    system_temp(1,1) = v0_10;
    m_air_0 = (p_0/R/T_air_i)*v0_10; %kg
    m_R_0 = m_bottle + rho_water*v0_w*.9 + m_air_0; %kg
        system_temp(1,2) = m_R_0;
        system_temp(1,3) = m_air_0;
        parameters_temp(14) = m_air_0;
    [t,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,...
    system_temp,parameters_temp),t,system_temp);
v0_10_x = max(dsdt(:,6));
v0_10_z = max(dsdt(:,7));
parameters_temp(15) = v0_15;
    system_temp(1,1) = v0_15;
    m_air_0 = (p_0/R/T_air_i)*v0_15; %kg
    m_R_0 = m_bottle + rho_water*v0_w*.85 + m_air_0; %kg
        system_temp(1,2) = m_R_0;
        system_temp(1,3) = m_air_0;
        parameters_temp(14) = m_air_0;
    [t,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,...
    system_temp,parameters_temp),t,system_temp);
v0_15_x = max(dsdt(:,6));
v0_15_z = max(dsdt(:,7));
parameters_temp(15) = v0_25;
    system_temp(1,1) = v0_25;
    m_air_0 = (p_0/R/T_air_i)*v0_25; %kg
    m_R_0 = m_bottle + rho_water*v0_w*.75 + m_air_0; %kg
        system_temp(1,2) = m_R_0;
        system_temp(1,3) = m_air_0;
        parameters_temp(14) = m_air_0;
    [t,dsdt] = ode45(@(t,system_temp) rocketTrajectory(t,...
    system_temp,parameters_temp),t,system_temp);
v0_25_x = max(dsdt(:,6));
v0_25_z = max(dsdt(:,7));


%Compare to original x and z values
v0_10_x = (v0_10_x-x)/x*100;
v0_10_z = (v0_10_z-z)/z*100;
v0_15_x = (v0_15_x-x)/x*100;
v0_15_z = (v0_15_z-z)/z*100;
v0_25_x = (v0_25_x-x)/x*100;
v0_25_z = (v0_25_z-z)/z*100;

%% Create tables
Variables = {'Inital Pressure';'Initial Volume Water';'Launch Angle';'Drag Coefficient'};
percent10x = [p0_10_x;v0_10_x;th_10_x;C_D_10_x];
percent10z = [p0_10_z;v0_10_z;th_10_z;C_D_10_z];

percent15x = [p0_15_x;v0_15_x;th_15_x;C_D_15_x];
percent15z = [p0_15_z;v0_15_z;th_15_z;C_D_15_z];

percent25x = [p0_25_x;v0_25_x;th_25_x;C_D_25_x];
percent25z = [p0_25_z;v0_25_z;th_25_z;C_D_25_z];

T1 = table(percent10x,percent10z,'RowNames',Variables)
T2 = table(percent15x,percent15z,'RowNames',Variables)
T3 = table(percent25x,percent25z,'RowNames',Variables)

%Parameter with most sensitivity in x is either the drag coefficient or
%initial pressure
%Parameter with most sensitivity in y is either the initial pressure or
%launch angle
end

