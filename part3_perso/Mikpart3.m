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

% we compare computed h with estimated h by gumbel curves

D = [1,3,6,12,24,48];
% output is 3 matrices of parameters, one for each return period
% for each return period, we have 3 parameters per event duration
paramT10 = zeros(3,6); % rows are parameters, columns are durations
paramT40 = zeros(3,6); % rows: c, e, f
paramT100 = zeros(3,6);

thres = [25,60,120];
for k = 1:3 % iterating over return periods [T]
    min_diff = thres(k);
    for l = 1:6 % iterating over event durations [D]
        for cspace = linspace(0,100,150)
            for fspace = linspace(-1, 1, 200)
                for espace = linspace(0, 1, 200)
                    hcomp = cspace*D(l)/(D(l).^espace + fspace); 
                    hgum = H_Gum(k,l);
                    sqerr = (hcomp - hgum).^2;
                    if sqerr < min_diff
                        min_diff = sqerr;
                        if k == 1
                            paramT10(1, l) = cspace;
                            paramT10(2, l) = espace;
                            paramT10(3, l) = fspace;
                        end
                        if k == 2
                            paramT40(1, l) = cspace;
                            paramT40(2, l) = espace;
                            paramT40(3, l) = fspace;
                        end
                        if k == 3
                            paramT100(1, l) = cspace;
                            paramT100(2, l) = espace;
                            paramT100(3, l) = fspace;
                        end
                   
                    end

                end
            end
        end
    end
end

%% second try

param = zeros(3,3);
threshold = [30 70 220];

for k = 1:3
    for cspace = linspace(0,100,150)
        for fspace = linspace(-1, 1, 200)
            for espace = linspace(0, 1, 200)
                for l = 1:6
                    hcomp = cspace*D(l)/((D(l).^espace)+fspace);
                    hgum = H_Gum(k,l);
                    square = (hcomp-hgum).^2;
                    if l==1
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

