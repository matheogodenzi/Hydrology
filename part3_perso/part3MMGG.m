% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 3: Construction of depth-duration-frequency (DDF) curves 
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1: Calibrate the DDF curve parameters
% useful functions: linspace
% -------------------------------------------------------------------------

% import the data from previous part
load assignment1_output_part2.mat 

% we search for parameters c, e, f, to build the DDF curves given by:
% h = c*D/(D**e + f)
% finding them by brute force, i.e. testing all possibilities and keeping
% the ones that matches the best

%% (1) Computing parameters

param = zeros(3,3);
threshold = [30 70 220];

for k = 1:3 % iterating over return periods T
    for cspace = linspace(0,100,150)
        for fspace = linspace(-1, 1, 200)
            for espace = linspace(0, 1, 200)
                for l = 1:6 % iterating over event durations
                    hcomp = cspace*D(l)/((D(l).^espace)+fspace);
                    hgum = H_Gum(k,l);
                    square = (hcomp-hgum).^2;
                    if l==1           % sum of squarred errors 
                        diff_squared = square;
                    end 
                    diff_squared = diff_squared + square;
                end
                if diff_squared <= threshold(k)
                    threshold(k) = diff_squared;
                    param(k,1) = cspace;
                    param(k,2) = espace;
                    param(k,3) = fspace;
                end 

            end
        end
    end 
end
%% (2) Output table with parameters and errors

%% (3) Plotting DDF curves and Gumbel estimation
D_prime = 0:1:100;
h_prime = zeros(length(D_prime), 3);
for k = 1:length(D_prime)
    for l = 1:3
        h_prime(k,l)= param(l,1)*D_prime(k)/(D_prime(k).^param(l,2)+ param(l,3));
    end 
end 

flipped = H_Gum.';

figure
plot(h_prime)