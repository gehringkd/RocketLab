function [] = bottleRocketSensitivity(t, system, parameters)
%UNTITLED Tests variation of bottle rocket flight as a function of initial
%angle, initial gage pressure inside bottle, initial volume of water, and
%the drag coefficient
%{
The purpose of this program is to test (individuall) the effect that each
of 4 parameters has on the distance x and height z of a bottle rocket 
launch. It does this by varying the initial condition and running that
iinitial condition through bottleRocketTrajectory and ode45.

Inputs: This program takes in the arrays t, system, and parameters. See
    bottleRocketSetParameters.m for more detail on what these contain or to 
    change values.
Outputs: This program outputs a graph???
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

%% Test variation of flight as a function of the drag coefficient, C_D
%C_D can only vary between 0.3 and 0.5, check all
C_D_range = linspace(0.3,0.5,20);
x_C_D = size(C_D_range);
y_C_D = size(C_D_range);

%make a modifiable parameters array
parameters_temp = parameters;

%Find final distance x and tallest height z for each C_D
    for i=1:size(C_D_range,2)
        parameters_temp(9) = C_D_range(i);
        [t,dsdt] = ode45(@(t,system) bottleRocketTrajectory(t,system,...
            parameters_temp),t,system);
        x_C_D(1,i) = max(dsdt(:,6));
        y_C_D(1,i) = max(dsdt(:,7));
    end

%Find best fit lines for comparisons against other 3 variables
p_x_C_D = polyfit(C_D_range,x_C_D,2);
p_y_C_D = polyfit(C_D_range,y_C_D,2);
x1_C_D = linspace(.3,.5);
slope_x_C_D = polyval(p_x_C_D,x1_C_D);
slope_y_C_D = polyval(p_y_C_D,x1_C_D);
    figure(2)

%Create plots
% Distance x vs. varying parameters
figure(2)
    scatter(C_D_range,x_C_D)
    hold on
    plot(x1_C_D,slope_x_C_D)
    xlabel('Value of C_D')
    ylabel('Distance Travelled (m)')
    title('Distance Variation as a Function of the Drag Coefficient')
figure(3)
    scatter(C_D_range,y_C_D)
    hold on
    plot(x1_C_D,slope_y_C_D)
    xlabel('Value of C_D')
    ylabel('Height Reached (m)')
    title('Height Variation as a Function of the Drag Coefficient')

%% Test variation of flight as a function of the initial angle, theta_0
%Theta can vary from 0 to 90 degrees, but 30 to 60 degrees is a reasonable
%angle
theta_range = linspace(30,60,30); %degrees, only 30 points to reduce run time
    theta_range = theta_range*pi/180; %radians
x_theta = size(theta_range);
y_theta = size(theta_range);

%make a modifiable system array
system_temp = system;

%Find final distance x and tallest height z for each C_D
    for i=1:size(theta_range,2)
        system_temp(1,4) = theta_range(i);
        [t,dsdt] = ode45(@(t,system_temp) bottleRocketTrajectory(t,...
            system_temp,parameters),t,system_temp);
        x_theta(1,i) = max(dsdt(:,6));
        y_theta(1,i) = max(dsdt(:,7));
    end

%Find best fit lines for comparisons against other 3 variables
p_x_theta = polyfit(theta_range, x_theta,4);
p_y_theta = polyfit(theta_range,y_theta,3);
theta_range2 = linspace(30,60)*pi/180; %wider range for polyval
slope_x_theta = polyval(p_x_theta,theta_range2);
slope_y_theta = polyval(p_y_theta,theta_range2);

%Create plots
% Distance x vs. varying parameters
figure(4)
    scatter(theta_range,x_theta)
    hold on
    plot(theta_range2,slope_x_theta)
    xlabel('Angle (radians)')
    ylabel('Distance Travelled (m)')
    title('Distance Variation as a Function of Initial Launch Angle')
figure(5)
    scatter(theta_range,y_theta)
    hold on
    plot(theta_range2,slope_y_theta)
    xlabel('Angle (radians)')
    ylabel('Height Reached (m)')
    title('Height Variation as a Function of Initial Launch Angle')

%% Test variation of flight as a function of the initial gage pressure
%p_gage can technically vary from 0 to gage_p_allow, but anything under 
%~2.5e5 Pa causes imaginary numbers and/or negative distances
gage_p_min = 35*6895; %psi to Pa
p_gage_range = linspace(gage_p_min,gage_p_allow,13); %35-100 psi-ish, test every 5 psi
    p_range_Pa = p_gage_range+p_ambient; %Pa
x_p = size(p_range_Pa);
y_p = size(p_range_Pa);

%reset modifiable parameters and system arrays
parameters_temp = parameters;
system_temp = system;

%Find final distance x and tallest height z for each C_D
    for i=1:size(p_range_Pa,2)
        parameters_temp(11) = p_range_Pa(i);
            m_air_0 = (p_range_Pa(i)/R/T_air_i)*v_0; %kg
            m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
        system_temp(1,2) = m_R_0;
        system_temp(1,3) = m_air_0;
        parameters_temp(14) = m_air_0;
        [t,dsdt] = ode45(@(t,system_temp) bottleRocketTrajectory(t,...
            system_temp,parameters_temp),t,system_temp);
        x_p(1,i) = max(dsdt(:,6));
        y_p(1,i) = max(dsdt(:,7));
    end

%Find best fit lines for comparisons against other 3 variables
[p_x_p,~,mu1] = polyfit(p_range_Pa, x_p,2);
[p_y_p,~,mu2] = polyfit(p_range_Pa,y_p,2);
p_range_psi =  (linspace(gage_p_min,gage_p_allow)+p_ambient)/6895;
slope_x_p = polyval(p_x_p,p_range_psi,[],mu1);
slope_y_p = polyval(p_y_p,p_range_psi,[],mu2);

%Create plots
% Distance x vs. varying parameters
figure(6)
    scatter(p_range_Pa,x_p)
    hold on
    plot(p_range_psi,slope_x_p)
    xlabel('Initial Air Pressure Inside Bottle (Pa)')
    ylabel('Distance Travelled (m)')
    title('Distance Variation as a Function of Initial Air Pressure')
figure(7)
    scatter(p_range_Pa,y_p)
    hold on
    plot(p_range_psi,slope_y_p)
    xlabel('Initial Air Pressure Inside Bottle (Pa)')
    ylabel('Height Reached (m)')
    title('Height Variation as a Function of Initial Air Pressure')
    
%% Test variation of flight as a function of the initial volume of water
%initial volume of water can vary from 0 to the volume of the bottle but
%anything under ~40% or over ~60% is more reasonable
vol_w_range = linspace(.4*vol_bottle,.6*vol_bottle,20); %m^3
    v_0_range = vol_bottle-vol_w_range; %m^3
x_v = size(v_0_range);
y_v = size(v_0_range);

%reset modifiable parameters and system arrays
parameters_temp = parameters;
system_temp = system;

%Find final distance x and tallest height z for each v_0
    for i=1:size(v_0_range,2)
        parameters_temp(15) = v_0_range(i);
        system_temp(1,1) = v_0_range(i);
            m_air_0 = (p_0/R/T_air_i)*v_0_range(i); %kg
            m_R_0 = m_bottle + rho_water*vol_w_range(i) + m_air_0; %kg
        system_temp(1,2) = m_R_0;
        system_temp(1,3) = m_air_0;
        parameters_temp(14) = m_air_0;
        [t,dsdt] = ode45(@(t,system_temp) bottleRocketTrajectory(t,...
            system_temp,parameters_temp),t,system_temp);
        x_v(1,i) = max(dsdt(:,6));
        y_v(1,i) = max(dsdt(:,7));
    end

%Find best fit lines for comparisons against other 3 variables
[v_x_p,~,mu3] = polyfit(v_0_range, x_v,3);
[v_y_p,~,mu4] = polyfit(v_0_range,y_v,3);
v_0_w_range = (linspace(.4*vol_bottle,.6*vol_bottle,20));
slope_x_v = polyval(v_x_p,v_0_w_range,[],mu3);
slope_y_v = polyval(v_y_p,v_0_w_range,[],mu4);

%Create plots
% Distance x vs. varying parameters
figure(8)
    scatter(vol_w_range,x_v)
    hold on
    plot(v_0_w_range,slope_x_v)
    xlabel('Initial Water Volume Inside Bottle (m^3)')
    ylabel('Distance Travelled (m)')
    title('Distance Variation as a Function of Initial Volume of Water')
figure(9)
    scatter(vol_w_range,y_v)
    hold on
    plot(v_0_w_range,slope_y_v)
    xlabel('Initial Water Volume Inside Bottle (m^3)')
    ylabel('Height Reached (m)')
    title('Height Variation as a Function of Initial Volume of Water')
end

