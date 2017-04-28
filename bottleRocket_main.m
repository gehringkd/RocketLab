%{
The purpose of this program is to access 4 functions to calculate the 
trajectory of a rocket under given conditions, find which conditions/
parameters have the greatest sensitivity for the flight, and report 
combinations of parameters that allow the bottle rocket to land a 
horizontal distance of 85 +- 1 m from the initial starting point.

Inputs: This program takes no inputs. Please note that if any initial 
    parameters need to be changed, those changes should be written in
    bottleRocketSetParameters.m
Outputs: This program outputs 2 graphs - one for the verification case,
    'verification_case_graph.png,' and one for any flight path that reaches
    83-86m, 'all_85m_flights.png.' It also prints a statement to the screen
    explaining which parameters have the greatest sensitivity for flight
    distance and height.
Assumptions: All processes during flight are adiabatic, and all outside
    forces (such as wind) are negligible. Also, that a 2L bottle can
    withstand 130psi, which was a number pulled off a google search.

Created by:	Kayla Gehring
Modified by:	Keith Covington
Created:	11/23/16
Modified:	04/12/17
%}
%Housekeeping
clear all;close all;clc;

%Get all parameters for flight (all came from verification data)
[parameters,system,state,t] = rocketParameters;

state = [system state];

%Get wind vector
windvector = wind;

%Calculate flight path using ode45
opts = odeset('Events',@stopping_point); % define event to stop ode45
[t,allStates] = ode45(@(t,state) rocketTrajectory(t,state,parameters,windvector, 45) ...
    ,t,state,opts);

% Display expectations
fprintf('\n');
disp(['Downrange distance achieved: ' num2str(allStates(end,7)) ' m.']);
disp(['Crossrange distance achieved: ' num2str(allStates(end,8)) ' m.']);
disp(['Height achieved: ' num2str(max(allStates(:,9))) ' m.']);
disp(['Total distance achieved: ' num2str(norm([allStates(end,7), allStates(end,8)])) ' m.']);
disp(['Angle: ' num2str(atand(allStates(end,8)/allStates(end,7))) ' degrees.']);

% Display results in English units
fprintf('\n');
disp(['Downrange distance achieved: ' num2str(allStates(end,7)*3.28084) ' ft.']);
disp(['Crossrange distance achieved: ' num2str(allStates(end,8)*3.28084) ' ft.']);
disp(['Height achieved: ' num2str(max(allStates(:,9))*3.28084) ' ft.']);
disp(['Total distance achieved: ' num2str(norm([allStates(end,7), allStates(end,8)])*3.28084) ' ft.']);
fprintf('\n');


%Graph flight path
figure(1)
plot3(allStates(:,7),allStates(:,8),allStates(:,9)); %plot(x,y,z)
grid on
axis equal
ylim([-10 10]);
title('Verification Case - Bottle Rocket Flight')
xlabel('Downrange distance (m)')
ylabel('Crossrange distance (m)')
zlabel('Vertical height (m)')


%% Monte Carlo Simulation
monteCarlo(t,state,parameters,windvector,opts);


%% Thrust Interpolation model
% Get static test stand data for Group 24
%thrustData = group24Thrust();

%plot(thrustData(:,2),thrustData(:,1))
%{

%Calculate flight path using ode45
[t,allStates] = ode45(@(t,state) thrustInterp(t,state,parameters,windvector,45,thrustData) ...
    ,t,state,opts);

%Graph flight path
figure
plot3(allStates(:,7),allStates(:,8),allStates(:,9)); %plot(x,y,z)
grid on
axis equal
ylim([-10 10]);
title('Verification Case - Bottle Rocket Flight')
xlabel('Downrange distance (m)')
ylabel('Crossrange distance (m)')
zlabel('Vertical height (m)')

%}

%Sensitivity Analysis
%varyVolWater(t, state, parameters);
%varyDragCoeff(t, state, parameters);
%varyDensity(t, state, parameters);
%varyPressure(t, state, parameters);
%varyAngle(t, state, parameters);

%TODO display delta distance achieved by each optimization

%{
% Analyze static test stand data
Isp_Avg = staticTests();

%Thrust Interpolation Model
T_v_t = group24Thrust();
tspan = [T_v_t(end,1), T_v_t(end,2)];
Isp = calcImpulse(T_v_t(:,1), T_v_t(:,2), 1);

%Rocket Equation/Isp model
    %find initial velocity
    m_initial = state(2);
    m_final = parameters(19); %mass bottle
    V1 = Isp*9.81*log(m_initial/m_final);
    V1 = [cosd(45)*V1, 0, cosd(45)*V1];
    state2 = [V1, 0, 0, 0.1];

opts = odeset('Events',@stopping_point2); % define event to stop ode45
t2 = [0,5];
[t2,allStates2] = ode45(@(t2,state2) rocketTrajectoryIsp(t2,state2,parameters,windvector,41) ...
    ,t2,state2,opts);


figure
plot3(allStates2(:,4),allStates2(:,5),allStates2(:,6)); %plot(x,y,z)
grid on
axis equal
ylim([-10 10]);
title('Bottle Rocket Flight, I_{sp} Model')
xlabel('Downrange distance (m)')
ylabel('Crossrange distance (m)')
zlabel('Vertical height (m)')
%}
