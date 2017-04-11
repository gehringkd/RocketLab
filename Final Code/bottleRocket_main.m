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
Modified:	04/04/17
%}
%Housekeeping
clear all;close all;clc;

%Get all parameters for flight (all came from verification data)
[parameters,system,state,t] = rocketParameters();

state = [system state];

%Calculate flight path using ode45
opts = odeset('Events', @stopping_point);
[t,dsdt] = ode45(@(t,state) rocketTrajectory(t,state,parameters) ...
    ,t,state,opts);

% display x- and y- coordinates
%disp([dsdt(:,8) dsdt(:,9) dsdt(:,10) dsdt(:,5:7)]);
disp([dsdt(:,8) dsdt(:,9) dsdt(:,10)]);


%Graph flight path
figure
plot(dsdt(:,8),dsdt(:,10)); %plot(x,z)
title('Verification Case - Bottle Rocket Flight')
xlabel('Horizontal distance (m)')
ylabel('Vertical height (m)')
%axis([0,55,0,18])

figure
plot(t,dsdt(:,10)); %plot(t,z)


%Find which parameter is most "sensitive" (which changes flight path most)
%rocketSensitivity(t,system,state,parameters);
 
%Find options to get 85m
%rocket85m(t,state,parameters);
   
% Vary the initial volume of water
%varyVolWater(t,state,parameters);


