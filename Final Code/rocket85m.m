function rocket85m(t,state,parameters)
%BOTTLEROCKET85m  The purpose of this code is to find a combination of 
%parameters that will allow the rocket to land within 1 m of 85m.

%Inputs: This program takes in the arrays t, system, and parameters. See
%    bottleRocketSetParameters.m for more detail on what these contain or to 
%    change values.
% Outputs: text to screen and a graph
% Assumptions: A 2L bottle can safely withstand 150psi, which was a number 
% pulled off a google search.
%
% Created by 4250df2b7e8b
% Created on 11/27/16
% Last edited 12/2/16
%% Parameters
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

t = [0 8];

%% Get state parameters passed into function
v_0 = state(1);
m_R_0 = state(2);
m_air_0 = state(3);
%th_0 = state(4);
V_0 = state(5);
x_0 = state(6);
z_0 = state(7);
Vx = state(5);  % x-component of velocity [m/s]
Vy = state(6);  % y-component of velocity [m/s]
Vz = state(7);  % z-component of velocity [m/s]
x = state(8);   % x-position [m]
y = state(9);   % y-position [m]
z = state(10);   % z-position [m]



%% Find a set of parameters that works
x = 0;
while x < 84 || x > 86
    
%increase pressure by 80%
p_0 = 1.8*p_0;
    if p_0 > (gage_p_allow+p_ambient)
        p_0 = gage_p_allow+p_ambient;
    end

%change parameters and systems arrays
parameters(11) = p_0;
    m_air_0 = (p_0/R/T_air_i)*v_0; %kg
    m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
state(2) = m_R_0;
state(3) = m_air_0;
parameters(14) = m_air_0;

%Find final distance x
[~,dsdt] = ode45(@(t,state) rocketTrajectory(t,...
state,parameters),t,state);
x = max(dsdt(:,6));

while x<84 || x>86
    if x<84 && C_D > 0.3
        %Reduce the drag coefficient
        C_D = round(.75*C_D*100)/100;
        if C_D < 0.3
           C_D = 0.3;
        end
    
        %Change parameter array
        parameters(9) = C_D;
    
        %Find final distance x
        [~,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,...
        parameters),t,state);
        x = max(dsdt(:,6));
    end
    
    if x > 86 && vol_water_i < 0.68*vol_bottle
        %Increase volume of water
        vol_water_i = vol_water_i+2.5e-5; %increase by 25 mL
            if vol_water_i > 0.68*vol_bottle %causes calculation errors
                vol_water_i = 0.68*vol_bottle
            end
        v_0 = vol_bottle-vol_water_i;
        
        
        %Change parameter and system arrays
        parameters(15) = v_0;
        state(1) = v_0;
        m_air_0 = (p_0/R/T_air_i)*v_0; %kg
        m_R_0 = m_bottle + rho_water*vol_water_i + m_air_0; %kg
        state(2) = m_R_0;
        state(3) = m_air_0;
        parameters(14) = m_air_0;
        [t,dsdt] = ode45(@(t,state) rocketTrajectory(t,...
        state,parameters),t,state);
        x = max(dsdt(:,6));
    end
    
    if vol_water_i == 0.68*vol_bottle && C_D == 0.3 && x > 86 && x < 84 && p_0 == gage_p_allow+p_ambient
        error('Add more parameters to calculations - this one is impossible.')
    end
end

    figure(2)
    plot(dsdt(:,6),dsdt(:,7)); %plot(x,z)
    title('Bottle Rocket Flight - 85m')
    xlabel('Horizontal distance (m)')
    ylabel('Vertical height (m)') 
    axis([0,90,0,30])
    
    fprintf('The rocket reaches %0.2f m when C_D, v_0, p_0, and th_0 are %0.2f,%0.5f,%0.0f,%0.0f respectively.\n',x, C_D,v_0,p_0/6895,th_0/pi*180')
end
end

