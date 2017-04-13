function [windVector] = wind
%Purpose: This function promped and recieves speed(mph) and direction of 
%the wind then converts to m\s and returns a 3-D vector.

%get inputs
mph = input('Enter a windspeed (mph): ');

if (mph == 0)
    windVector = [0; 0; 0];
else
    
direction = input('Enter a direction (lowercase): ','s');

%convert to m/s
windSpeed = mph*1609.34/3600;

%determine wind angle
switch direction
    case 'n'
        angle = 210;
    case 'nne'
        angle = 187.5;
    case 'ne'
        angle = 165;
    case 'ene'
        angle = 142.5;
    case 'e'
        angle = 120;
    case 'ese'
        angle = 97.5;
    case 'se'
        angle = 75;
    case 'sse'
        angle = 52.5;
    case 's'
        angle = 30;
    case 'ssw'
        angle = 7.5;
    case 'sw'
        angle = -15;
    case 'wsw'
        angle = -37.5;
    case 'w'
        angle = -60;
    case 'wnw'
        angle = -82.5;
    case 'nw'
        angle = -105;
    case 'nnw'
        angle = -127.5;
    otherwise
        warning('Not a acceptable wind direction')
end

%create vector for wind
windVector(1) = windSpeed*cosd(angle);
windVector(2) = windSpeed*sind(angle);
windVector(3) = 0;  %Z is always zero
windVector = windVector.';
end

