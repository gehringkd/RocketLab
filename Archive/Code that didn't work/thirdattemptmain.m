rho_w = 1000;
R = 287;
T_0 = 300;
p_0 = (50+12.03)*6894.75729;

m_B = 0.07;
v_B = .002;
v_w = .001;
v_0 = v_B-v_w;

m_0 = m_B+rho_w*v_w+p_0/R/T_0*v_0;

V_0 = 0;
x_0 = 0;
z_0 = 0.1;
th_0 = 45*(pi/180);

s = [v_0 m_0 V_0 x_0 z_0 th_0];

ti = 0;
t1 = 0.4;

[t,dsdt] = ode45('phase1_3rd_attempt',[ti 5], s);