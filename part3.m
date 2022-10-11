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

for cspace = linspace(0,100,150)
    for fspace = linspace(-1, 1, 200)
        for espace = linspace(0, 1, 200)
        
        
        
        end
    end
end



