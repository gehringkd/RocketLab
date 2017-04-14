function [windVector] = wind()
%Purpose: This function promped and recieves speed(mph) and direction of 
%the wind then converts to m\s and returns a 3-D vector.

% Define directions
direcs = {'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'};
angles = linspace(210,-127.5,16);

%get inputs
mph = input('Enter a windspeed (mph): ');
if (mph == 0)
	windVector = [0; 0; 0];
else
	direction = input('Enter a direction (ex: nne): ','s');
	direction = upper(direction); % convert to uppercase

	% Check to make sure input direction is valid
	while ~any(strcmp(direcs, direction))
		disp('Not an acceptable direction. Try again...');
		direction = input('Enter a direction (ex: nne): ','s');
		direction = upper(direction); % convert to uppercase
	end

	%convert to m/s
	windSpeed = mph*1609.34/3600;

	%determine wind angle
	switch direction
	    case direcs{1}
		angle = angles(1);
	    case direcs{2}
		angle = angles(2);
	    case direcs{3}
		angle = angles(3);
	    case direcs{4}
		angle = angles(4);
	    case direcs{5}
		angle = angles(5);
	    case direcs{6}
		angle = angles(6);
	    case direcs{7}
		angle = angles(7);
	    case direcs{8}
		angle = angles(8);
	    case direcs{9}
		angle = angles(9);
	    case direcs{10}
		angle = angles(10);
	    case direcs{11}
		angle = angles(11);
	    case direcs{12}
		angle = angles(12);
	    case direcs{13}
		angle = angles(13);
	    case direcs{14}
		angle = angles(14);
	    case direcs{15}
		angle = angles(15);
	    case direcs{16}
		angle = angles(16);
	    otherwise
		error('Not a acceptable wind direction.')
	end

	%correct error in above angles
	angle = angle + 180;

	%create vector for wind
	windVector(1) = windSpeed*cosd(angle);
	windVector(2) = windSpeed*sind(angle);
	windVector(3) = 0;  %Z is always zero
	windVector = windVector.';

	end

end
