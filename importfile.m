function [dataOther] = importfile(dataArray)
%IMPORTFILE Import numeric data from an excel file as column vectors.
%Subtract the 0 m/s data from 25 m/s data.
%
% Returns dynamic pressure, velocity, angle of attack, N, A

%Note to self: error correction should probably be added here

% Created by Kayla Gehring, 4/13

%% Open the file and extract required data
dataArray = load(dataArray);

%Extract only required variables
    %To calculate needed coefficients, only need v, p_dyn, AoA L, D    
    %This is columns D, E, and W-Y (4-5, 23-25)
    dataArray = dataArray(:,[4:5,23:25]);
        %round the velocities for ease of subtracting V=0 from other values
        dataArray(:,1) = round(dataArray(:,1));


 %% Take average of data for each velocity
%Average the data
%5 velocities total, 20 samples each
index = size(dataArray,1)/20; %20 data points for each angle
averaged_data = zeros(index,size(dataArray,2));

j=0;
for i=1:index
        %average the 20 samples taken
        averaged_data(i,:) = mean(dataArray(j+1:j+19,:));
        j=j+20;
end


%Separate zero velocity
data0 = averaged_data(1, :);
for i=2:size(averaged_data,2)
    dataOther(i-1,:) = averaged_data(i, :) - data0;
end
        
%Close file
fclose('all');

end


