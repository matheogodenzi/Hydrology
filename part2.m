% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 2: Fit a Gumbel distribution and calculate critical rainfall depths
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1: Compute the Weibull plotting position
% -------------------------------------------------------------------------
%% (1) Computing Fh nd Yf
% import the data from Part1 using the function load
load assignment1_output_part1.mat

% sorting data columns by amplitude
sortedAnnualMax = sort(AnnualMax,1, 'descend');
sortedAnnualMax;

% initializing empty vector
Fh = zeros([size(sortedAnnualMax, 1),1]);

% computing frequency : Fh = i/(N+1)
for c = 1:size(Fh,1)
    Fh(c) = (size(Fh,1) - (c-1))/ (size(Fh,1) + 1);
end
Fh;
% computing reduced variable Yf

Yf = zeros([size(Fh, 1),1]);

for c = 1:size(Yf,1)
    Yf(c) = -log(-log(Fh(c)));
end

% Plot of Fh vs h
% for first column
plot(sortedAnnualMax(:,1), Fh, 'o');

%% (2) Fitting Gumbel curve
% method of moments

% Gumbel method
