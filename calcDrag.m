function [ Drag ] = calcDrag( data, S)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Extract/name variables
V = data(:,1);
dyn_p = data(:,2);
AoA = data(:,3)*pi/180; %radians
N = data(:,4);
A = data(:,5);

%% Calculate L and D
D = N.*sin(AoA) + A.*cos(AoA);

%% Calculate the coefficient
C_D = D./dyn_p/S;

%% Pass into D_vs_Aoa
Drag = mean(C_D);

end

