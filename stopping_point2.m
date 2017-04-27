function [value,isterminal,direction] = stopping_point2(t, state)
%stopping_point Does this and that.
%{
The purpose of this program is to 

Inputs: t       - time
        state	- state

Outputs: value		- value
         isterminal	- isterminal
	 direction	- direction

Created by:     Keith Covington
Created on:     04/08/2017
Last modified:  04/08/2017
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

value = state(6); % detect when z=state(10) < 0
isterminal = 1; % stop the integration
direction = -1; % negative direction

end
