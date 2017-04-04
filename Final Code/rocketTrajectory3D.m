function rocketTrajectory3D()
%Changes needed for 3D trajectory

%% Parameters
g = [0;0;9.81]; %[x,y,z] m/s^2

%% Define V_rel so F can be calculated throughout phases
V_wind = somefunction(p1, p2); %[x,y,z] m/s

%V brought in same as before, except now as a vector. Ode45 will need to
%calculate V, and V_rel needs to be used to calculate T and D to calculate
%dVdt
V_rel = V - V_wind;
V_mag = norm(V_rel);

%Unit vector h, direction of V_rel
h = V_rel/V_mag;

%Thrust F calculated through the 3 phases is still valid, it is just  a
%magnitude for now
%% Trajectory, 3D
    
    %Drag equation, magnitude
    D = (rho_a/2)*(V_mag^2)*c_D*A_B;
    
    %Change D and T to vectors instead of magnitudes
    D = D*h;
    F = F*h;
    
    sumF = F-D+g; %vector in x,y,z
    
    %derivatives (all vectors [x;y;z];
    dVdt = sumF/m_R;
    dxdt = V;
    
    
    
end
