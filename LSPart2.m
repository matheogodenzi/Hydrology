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

% import the data from Part1 using the function load
load assignment1_output_part1.mat


SortedAnnualMax = sort(AnnualMax, 1, 'ascend');

%SortedAnnualMax

N = length(SortedAnnualMax) ;
Fh = zeros(N,1) ;
RedVar = zeros(N,1) ;

for i = 1:N
    Fh(i,1) = i/(N+1) ;
    RedVar(i,1) = -log10(-log10(Fh(i,1))) ; 
end

SortedAnnualMax(:,1)

figure
plot(SortedAnnualMax(:,i), Fh, "o");


title('Empirical Frequency');
xlabel('Precipitations depths');
ylabel('Empirical Frequency');
