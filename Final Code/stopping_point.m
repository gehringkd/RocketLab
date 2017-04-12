function [value,isterminal,direction] = stopping_point(t, state)
%stopping_point defines an event for ode45 to stop numerical integration 
% of rocket flight.
%{
The purpose of this program is to define when ode45 should stop running 
 i.e. when the rocket hits the ground.

Inputs: t       - time
        state	- state

Outputs: value		- value
         isterminal	- isterminal
	 direction	- direction

Created by:     Keith Covington
Created on:     04/08/2017
Last modified:  04/12/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value = state(9); % detect when z=state(9) < 0
isterminal = 1; % stop the integration
direction = -1; % negative direction

end
