%% Knowns that are unlikely to change 
gamma = 1.4; %specific heat ratio of air
rho_w = 1000; %kg/m^3, density of water
R = 287; % J/kg/K
g = 9.81; %m/s^2

%% Knowns likely to change - look here to change parameters of experiment
p_gage = 50*6894.76; %Pa, initial gage pressure of air in bottle
Vol_B = 0.002; %m^3, volume of the bottle
Vol_w_i = 0.001; %m^3, initial volume of water
m_B = 0.07; %kg, mass of empty bottle
c_D = 0.5; %drag coefficent

%Atmosphere paramters
p_a = 12.03*6894.76; %psi to Pa, ambient air pressure
T_air_i = 300; %K, initial temperature of air

%inital mass of entire bottle
    v_0 = Vol_B - Vol_w_i;
    p_0 = p_gage + p_a;
    m_total_i = m_B+rho_w*(Vol_w_i)+(p_0*v_0/(R*T_air_i)); %kg
    m_air_i = (p_0/(R*T_air_i))*v_0; %kg
    m_water_i = rho_w*Vol_w_i;

%% Set up Verification Case Run
% Pass t = [t0 tf] and s = [v,m,th,V,x,z] to bottleRocketTrajectory.m
m_0 = m_total_i; %kg
th_0 = 45*(pi/180); %radians
V_0 = 0; %m/s
x_0 = 0; %m
y_0 = 0.1; %m
s = [v_0 m_0 m_air_i th_0 V_0 x_0 y_0 p_0 c_D];

t_0 = 0; %s
t_f = 5; %s
t = [t_0 t_f];

%Run bottleRocketTrajectory.m using ode45
[t,dsdt] = ode45('bottleRocketTrajectory',t,s);

%Extract necessary data
x = dsdt(:,6);
y = dsdt(:,7);

%Plot verification case
plot(x,y)
title('Rocket Trajectory')
xlabel('Horizontal Distance (m)')
ylabel('Height (m)')
axis([0,55,0,18])

%% Find Most Sensitive Parameter
%To see which parameter is most sensitive, vary each parameter by 5% and
%test the final distance x to see which makes the greatest change

%Extract x_final from verification case
x_final = max(x);
y_max = max(y);

%Test difference that launch angle makes (th_0)
th_1 = 0.95*th_0; %theta-5%theta
s = [v_0 m_0 m_air_i th_1 V_0 x_0 y_0 p_0 c_D];
[t,dsdt] = ode45('bottleRocketTrajectory',t,s);
x_th = max(dsdt(:,6));
y_th = max(dsdt(:,7));

%Test difference that amount of water makes (v_0)
Vol_w_1 = 0.95*Vol_w_i; %make water volume 110% original volume
v_1 = Vol_B - Vol_w_1;
s = [v_1 m_0 m_air_i th_0 V_0 x_0 y_0 p_0 c_D];
[t,dsdt] = ode45('bottleRocketTrajectory',t,s);
x_v = max(dsdt(:,6));
y_v = max(dsdt(:,7));

%Test difference that initial air pressure (gage) makes (p_0)
max_p_gage = 130; %some random value I found
F_safety = 1.5;
p_gage_allow = F_safety*max_p_gage;
p_allow = p_gage_allow+p_a;

p_gage_1 = 0.95*p_gage;
p_1 = p_gage_1 + p_a;
s = [v_0 m_0 m_air_i th_0 V_0 x_0 y_0 p_1 c_D];
[t,dsdt] = ode45('bottleRocketTrajectory',t,s);
x_p = max(dsdt(:,6));
y_p = max(dsdt(:,7));

    %alternatively, find rate of change that occurs as a function of the
    %drag coefficient
    p_opts = linspace(0,p_allow,30);
    x_p = zeros(1,size(c_D_opts,2));
    y_p = zeros(1,size(c_D_opts,2));
    for i=1:size(c_D_opts,2)
        s = [v_0 m_0 m_air_i th_0 V_0 x_0 y_0 p_opts(i) c_D];
        [t,dsdt] = ode45('bottleRocketTrajectory',t,s);
        x_p(1,i) = max(dsdt(:,6));
        y_p(1,i) = max(dsdt(:,7));
    end
    figure(4)
    p = polyfit(p_opts,x_p,2)
    x1 = linspace(3,5);
    y1 = polyval(p,x1);
    scatter(c_D_opts,x_c_D)
    hold on
    plot(x1,y1)
    xlabel('value of p')
    ylabel('distance travelled')

%Test difference that drag coefficient makes (c_D)
c_D_1 = 0.95*c_D; 
s = [v_0 m_0 m_air_i th_0 V_0 x_0 y_0 p_0 c_D_1];
[t,dsdt] = ode45('bottleRocketTrajectory',t,s);
%x_c_D = max(dsdt(:,6));
%y_c_D = max(dsdt(:,7));

%{
    %alternatively, find rate of change that occurs as a function of the
    %drag coefficient
    c_D_opts = linspace(3,5,30);
    x_c_D = zeros(1,size(c_D_opts,2));
    y_c_D = zeros(1,size(c_D_opts,2));
    for i=1:size(c_D_opts,2)
        s = [v_0 m_0 m_air_i th_0 V_0 x_0 y_0 p_0 c_D_opts(1,i)];
        [t,dsdt] = ode45('bottleRocketTrajectory',t,s);
        x_c_D(1,i) = max(dsdt(:,6));
        y_c_D(1,i) = max(dsdt(:,7));
    end
    p = polyfit(c_D_opts,x_c_D,2)
    x4 = linspace(3,5);
    y4 = polyval(p,x1);
    figure(5)
    scatter(c_D_opts,x_c_D)
    hold on
    plot(x1,y1)
    xlabel('value of C_D')
    ylabel('distance travelled')

%Calculate percent differences
diff_th_x = (x_th - x_final)/x_final*100;
diff_v_x = (x_v - x_final)/x_final*100;
diff_p_x = (x_p - x_final)/x_final*100;
diff_c_D_x = (x_c_D - x_final)/x_final*100;
diff_th_y = (y_th - y_max)/y_max*100;
diff_v_y = (y_v - y_max)/y_max*100;
diff_p_y = (y_p - y_max)/y_max*100;
diff_c_D_y = (y_c_D - y_max)/y_max*100;

%Find greatest difference
fprintf('Decreasing the launch angle 5 percent causes a %0.2f percent change in distance x and a %0.2f percent change in height y.\n', diff_th_x,diff_th_y);
fprintf('Decreasing the initial volume of water 5 percent causes a %0.3f percent change in distance x and a %0.2f percent change in height y..\n', diff_v_x,diff_v_y);
fprintf('Decreasing the initial gage pressure 5 percent causes a %0.3f percent change in distance x and a %0.2f percent change in height y..\n', diff_p_x,diff_p_y);
fprintf('Decreasing the drag coefficient 5 percent causes a %0.3f percent change in distance x and a %0.2f percent change in height y..\n', diff_c_D_x,diff_c_D_y);

%% Make the rocket land at 85m +- 1 m
%Find the percent increase needed to get to 85 m
diff_x_desired = (85-x_final)/x_final*100;

%The greatest possible increase in x comes decreasing the initial volume 
%of water