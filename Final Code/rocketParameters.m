function [parameters,system,state,t] = rocketParameters()
%BOTTLEROCKETSETPARAMATERS  sends all parameters for a bottle rocket to 
%bottleRocketmain.m.
%   To change any parameters, just edit this function.
% Inputs: none
% Outputs: all paramaters as an array
% Assumptions: A 2L bottle can safely withstand 180psi, which was a number 
% pulled off a google search.
%
% Created by 4250df2b7e8b
% Created on 11/27/16
% Last edited 11/27/16

%% Name all parameters for ease of changing them later
    %Earth constants - unlikely to change
    g = 9.81; %m/s^2
    gamma = 1.4; %no units; ratio of specific heats for air
    rho_water = 1000; %kg/m^3
    R = 287; %J/kg/K, gas constant for air
    
    %Atmospheric constants - somewhat likely to change by small amounts
    rho_air_ambient = 0.961; %kg/m^3
    p_ambient = 12.03; %psi
        p_ambient = p_ambient*6895; %Pa
    T_air_i = 300; %K

    % --------------------------- Wind ------------------------------
    %V_wind = [0; 0]; % wind velocity components in [x; z] [m/s]

    % For when we have 3D working...
    V_wind = [0; 0; 0]; % wind velocity components in [x; y; z] [m/s]
    % ---------------------------------------------------------------
    
    %Flight parameters - likely to change
    C_d = 0.8; %discharge coefficient
    C_D = 0.5; %drag coefficient
    p_gage = 50; %psi, initial gage pressure
        p_gage = p_gage*6895; %Pa
        p_0 = p_ambient + p_gage; %Pa
    vol_water_i = 0.001; %m^3, initial volume of water
    
    %Parameters specific to bottle - could change
    vol_bottle = 0.002; %m^3 = 2 L
    d_throat = .021; %m
        A_throat = (d_throat/2)^2*pi; %m^2
    d_bottle = .105; %m
        A_bottle = (d_bottle/2)^2*pi; %m^2
    m_bottle = 0.07; %kg; mass of empty bottle
    
    %Factor of safety and max pressure
    max_gage_p = 150; %psi
        max_gage_p = max_gage_p*6895; %Pa
    FOS = 1.3;
        gage_p_allow = max_gage_p/FOS;
        
%% Initial conditions of equations calculated by ode45
    v_0 = vol_bottle - vol_water_i; %m^3
    m_air_0 = (p_0/R/T_air_i)*(v_0); %kg
    m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
    theta_0 = 45; %degrees
        theta_0 = theta_0*pi/180; %radians
    V_0 = [0.001; 0; 0.001]; %m/s
    x_0 = 0; %m
    y_0 = 0; %m
    z_0 = 0.1; %m

%% Time range
    ti = 0;
    tf = 5;
    
%% Create arrays
%Put initial system conditions in outgoing array
    system = [v_0 m_R_0 m_air_0 theta_0];

 %Put initial state conditions in outgoing array
    state = [V_0(1) V_0(2) V_0(3) x_0 y_0 z_0];
    
%Put all useful parameters in outgoing array (some were used to calculate
%other paramters and were not included)
    parameters = [g; gamma; rho_water; R; rho_air_ambient; p_ambient; ...
        T_air_i; C_d; C_D; p_gage; p_0; vol_water_i; vol_bottle; ... 
        m_air_0; v_0; A_throat; A_bottle; gage_p_allow; m_bottle; V_wind];
%Put time range in outgoing array
    t = [ti tf];
end

