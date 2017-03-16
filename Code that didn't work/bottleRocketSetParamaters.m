function [parameters,system,t] = bottleRocketSetParameters()
%BOTTLEROCKETSETPARAMATERS  sends all parameters for a bottle rocket to 
%bottleRocketmain.m.
%   To change any parameters, just edit this function.
% Inputs: none
% Outputs: all paramaters as an array
% Assumptions: none apply to this function, see main for assumptions in
% problem
%
% Created by UID!!!
% Created on 11/27/16
% Last edited 11/27/16

%Name all parameters for ease of changing them later
    %Earth constants - unlikely to change
    g = 9.81; %m/s^2
    gamma = 1.4; %no units; ratio of specific heats for air
    rho_water = 1000; %kg/m^3
    R = 287; %J/kg/K, gas constant for air
    
    
    %Atmospheric constants - somewhat likely to change by small amounts
    rho_air_ambient = 0.961; %kg/m^3
    p_ambient = 12.03; %psi
        p_ambient = p_ambient*6894.76; %Pa
    T_air_i = 300; %K
    
    %Flight parameters - likely to change
    C_d = 0.8; %discharge coefficient
    C_D = 0.5; %drag coefficient
    p_gage = 50; %psi, initial gage pressure
        p_gage = p_gage*6894.76; %Pa
        p_0 = p_ambient + p_gage; %Pa
    vol_water_i = 0.001; %m^3, initial volume of water
    
    %Parameters specific to bottle - could change
    vol_bottle = 0.002; %m^3 = 2 L
    d_throat = .021; %m
        A_throat = (d_throat/2)^2*pi; %m^2
    d_bottle = .105; %m
        A_bottle = (d_bottle/2)^2*pi; %m^2
    m_bottle = 0.07; %kg; mass of empty bottle


%Put all useful parameters in outgoing array (some were used to calculate
%other paramters and were not included)
    parameters = [g gamma rho_water R rho_air_ambient T_air_i C_d C_D ...
    p_gage p_0 vol_water_i ti tf vol_bottle A_throat ...
    A_bottle m_bottle];

%% Initial conditions of equations calculated by ode45
    v_0 = vol_bottle - vol_water_i; %m^3
    m_air_0 = (p_0/R/T_air_i)*v_0; %kg
    m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
    theta_0 = 45; %degrees
        theta_0 = theta_0*pi/180; %radians
    V_0 = 0; %m/s
    x_0 = 0; %m
    z_0 = 0.1; %m

%Put initial system conditions in outgoing array
    system = [v_0 m_R_0 m_air_0 theta_0 V_0 x_0 z_0];

%% Time range
    ti = 0;
    tf = 5;

%Put time range in outgoing array
    t = [ti tf];
end

