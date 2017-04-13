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

[t,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,parameters,windvector, 45) ...
    ,t,state);

%Graph flight path
figure(1)
plot3(dsdt(:,7),dsdt(:,8),dsdt(:,9)); %plot(x,y,z)
grid on
axis equal
ylim([-10 10]);
title('Verification Case - Bottle Rocket Flight')
xlabel('Downrange distance (m)')
ylabel('Crossrange distance (m)')
zlabel('Vertical height (m)')
   

%Sensitivity Analysis
varyVolWater(t, state, parameters);
%varyDragCoeff(t, state, parameters);
varyDensity(t, state, parameters);
%varyPressure(t, state, parameters);
%varyAngle(t, state, parameters);



