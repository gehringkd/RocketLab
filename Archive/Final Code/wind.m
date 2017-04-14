function [windVector] = wind()
%Purpose: This function promped and recieves speed(mph) and direction of 
%the wind then converts to m\s and returns a 3-D vector.

% Define directions
direcs = {'N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW'};


%get inputs
mph = input('Enter a windspeed (mph): ');
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
        angle = 210;
    case direcs{2}
        angle = 187.5;
    case direcs{3}
        angle = 165;
    case direcs{4}
        angle = 142.5;
    case direcs{5}
        angle = 120;
    case direcs{6}
        angle = 97.5;
    case direcs{7}
        angle = 75;
    case direcs{8}
        angle = 52.5;
    case direcs{9}
        angle = 30;
    case direcs{10}
        angle = 7.5;
    case direcs{11}
        angle = -15;
    case direcs{12}
        angle = -37.5;
    case direcs{13}
        angle = -60;
    case direcs{14}
        angle = -82.5;
    case direcs{15}
        angle = -105;
    case direcs{16}
        angle = -127.5;
    otherwise
        error('Not a acceptable wind direction.')
end

%create vector for wind
windVector(1) = windSpeed*cosd(angle);
windVector(2) = windSpeed*sind(angle);
windVector(3) = 0;  %Z is always zero
windVector = windVector.';

end
