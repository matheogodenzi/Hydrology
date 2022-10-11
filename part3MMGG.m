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

parameters = zeros(18,3); %rows : 3 x 6 durations (corresponding to all 
% the durations we have for a given return period) cols: parameters  

threshold = [25 60 120];
for k = 1:3
    for l = 1:6
        min_diff = threshold(k);
        for cspace = linspace(0,100,150)
            for fspace = linspace(-1, 1, 200)
                for espace = linspace(0, 1, 200)
                    hcomp = cspace*D(l)/(D(l).^espace+fspace);
                    hgum = H_Gum(k,l);
                    diff_squared = (hcomp-hgum).^2;
                    if   diff_squared < min_diff
                        min_diff =  diff_squared;
                        parameters(k*l,1) = cspace;
                        parameters(k*l,2) = espace;
                        parameters(k*l,3) = fspace;
                    end 
                end
            end
        end
    end
end 

parameters
