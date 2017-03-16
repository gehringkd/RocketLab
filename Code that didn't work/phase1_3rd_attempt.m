function [dsdt] = phase1_3rd_attempt(t,s)
%% Constants and equations
gamma = 1.4;
rho_w = 1000;
g = 9.81;

c_d = 0.8;
c_D = 0.5;
A_t = pi*(0.0105^2); %m^2
A_B = pi*(0.0525^2); %m^2
p_a = 12.03*6894.75729; %Pa
p_0 = 50*6894.75729+p_a; %Pa
v_0 = .001; %Change this to be in terms of bottle and water
v_B = .002; %m^3
rho_a = 0.961;

global p_end


%% Input array
s = s';
v = s(1);
m = s(2);
V = s(3);
x = s(4);
z = s(5);
th = s(6);

%% Phase 1

if v < v_B
    dvdt = c_d*A_t*sqrt(2*(p_0*((v_0/v)^gamma)-p_a));
    dmdt = -c_d*A_t*sqrt(2*rho_w*(p_0*((v_0/v)^gamma)-p_a));

    F = 2*c_d*A_t*(p_0*((v_0/v)^gamma)-p_a);
    D = rho_a/2*V^2*c_D*A_B;
    
    dVdt = F-D-m*g*sin(th);
    dxdt = V*cos(th);
    dzdt = V*sin(th);
    if V < 1
        dthdt = 0;
    else
        dthdt = -g*cos(th)/V;
    end
    p_end = (p_0*((v_0/v)^gamma)-p_a);
else
    dvdt = 0;
    dmdt = 0;
    dVdt = 0;
    dxdt = 0;
    dzdt = 0;
    dthdt = 0;
end

dsdt(1) = dvdt;
dsdt(2) = dmdt;
dsdt(3) = dVdt;
dsdt(4) = dxdt;
dsdt(5) = dzdt;
dsdt(6) = dthdt;
dsdt = dsdt';

end