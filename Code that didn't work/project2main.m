%% Knowns that are unlikely to change 
global gamma rho_w R g d_B p_gage Vol_B Vol_w_i m_B c_D c_d rho_a p_a T_air_i A_t A_B v_0

gamma = 1.4; %specific heat ratio of air
rho_w = 1000; %kg/m^3, density of water
R = 287; % J/kg/K
g = 9.81; %m/s^2

%% Knowns likely to change - look here to change parameters of experiment
%Bottle parameters
d_B = 10.5/100; %cm to m, diameter of bottle
d_t = 2.1/100; %cm to m, diameter of throat of bottle
p_gage = 50*6894.76; %Pa, initial gage pressure of air in bottle
Vol_B = 0.002; %m^3, volume of the bottle
Vol_w_i = 0.001; %m^3, initial volume of water
m_B = 0.07; %kg, mass of empty bottle

%Flight parameters
c_D = 0.5; %drag coefficent
c_d = 0.8; %discharge coefficient

%Atmosphere paramters
rho_a = 0.961; %kg/m^3, density of ambient air
p_a = 12.03*6894.76; %psi to Pa, ambient air pressure
T_air_i = 300; %K, initial temperature of air
%inital mass of entire bottle
    m_total_i = m_B+rho_w*(Vol_w_i)+(p_air_i/(R*T_air_i))*v_air_i; %kg

%equations based off changeable parameters    
A_t = pi*(d_t/2)^2; %m^2, area of throat
A_B = pi*(d_B/2)^2; %m^2, area of bottle 


% s = [v,m,th,V,x,z]
    %v = initial volume of air
    %m = initial mass of bottle
    %th = theta, angle of rocket launch in radians
    %V = initial velocity, probably 0
    %x = initial horizontal distance, probably 0
    %z = initial height, make the ground 0! Calculations stop when z<0
    global v_o m_o
v_0 = Vol_B - Vol_w_i;
m_0 = m_total_i; %kg
th = 45*(pi/180); %radians
V_0 = 0; %m/s
x_0 = 0; %m
y_0 = 0.1; %m

s1 = [v_0 m_0 th V_0 x_0 y_0];
t_0 = 0; %s
t_f = 5; %s
t = [t_0 t_f];

%% Calculate phase 1
[t1,dsdt1] = ode45('rocketPhase1',t,s1);
for i=1:size(t1);
    if dsdt(i+2,1)-dsdt(i,1)
        ind1 = i;
        break
    end
end
dsdt1 = dsdt1(1:ind1,:);
t1 = t1(1:ind1,:);

t = [t_end_phase1 5];

%% Calculate phase 2
s2 = dsdt1(ind1,:);
[t2,dsdt2] = ode45('rocketPhase2',t,s2);
global t_end_phase2
t = [t_end_phase2 2];
ind2 = find(t==t_end_phase2);
dsdt2 = dsdt2(1:ind2,:);
t2 = t2(2:ind,:);

%% Calculate phase 3
s3 = dsdt2(ind2,:);
[t3,dsdt3] = ode45('rocketPhase3',t,s3);

%% Combine all 3 trajectories
dsdt = [dsdt1;dsdt2;dsdt3];
t = [t1;t2;t3];
x = dsdt(:,5);
y = dsdt(:,6);

plot(x,y)
xlabel('Horizontal Distance (m)')
ylabel('Height (m)')
title('Bottle Rocket Trajectory')