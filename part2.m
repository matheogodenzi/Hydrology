% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 2: Fit a Gumbel distribution and calculate critical rainfall depths
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1: Compute the Weibull plotting position
% -------------------------------------------------------------------------
%% (1) Computing Fh nd YF
% import the data from Part1 using the function load
load assignment1_output_part1.mat

% sorting data columns by amplitude
sortedAnnualMax = sort(AnnualMax,1, 'ascend');
sortedAnnualMax;

%% Computing Weibull coefficients 

% initializing empty vector
Fh = zeros([size(sortedAnnualMax, 1),1]);

% computing frequency : Fh = i/(N+1)
sz = size(Fh,1);
for c = 1:sz
    Fh(c) = c/ (sz + 1);
end
Fh;

%% computing reduced variable YF

YF = -log(-log(Fh));

% setting up colors for all the plots
newcolors = [0.6350 0.0780 0.1840; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
        0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];

% Plot of Fh vs h
figure(21)
colororder(newcolors)
plot(sortedAnnualMax(), Fh, '.');
title('Empirical Frequencies vs Precipitation Depth ')
xlabel('Precipitation Depth [mm]') 
ylabel('Empirical Frequencies [Fh]')
lgd = legend({'1 hour', '3 hours', '6 hours', ...
        '12 hours', '24 hours', '48 hours'}, ...
        "Location", "southeast", "NumColumns",2);
title(lgd, "Annual Max depth over time span [mm]")
print('FhvsD','-vector','-dpdf') % this saves 'my_figure.pdf' (useful for LaTeX)

%% (2) Fitting Gumbel curve
% method of moments
% calculating means and standard deviation for each duration (i.e. each
% column) 
std_vect = std(AnnualMax);
mean_vect = mean(AnnualMax);

%defining variables 
GumbelParMoments = zeros(2,6);
m = length(std_vect);

for k = 1:m %iterating over all durations
    alpha = pi/(std_vect(k)*sqrt(6)); %defining the parameter alpha 
    u = mean_vect(k)- 0.5772/alpha; %defining the parameter u 

    GumbelParMoments(1,k) = alpha;
    GumbelParMoments(2,k) = u;
end 
GumbelParMoments;
%%
% Gumbel method
GumbelPar = zeros(2,6);
m = length(std_vect);

%calculating the means and standard deviations of the reduced variable (Fh)
mu_YF = mean(YF);
sigma_YF = std(YF);

% Defining paramaters 
for k = 1:m
    alphaG = sigma_YF/std_vect(k);
    uG = mean_vect(k)-(mu_YF/sigma_YF)*std_vect(k);
    
    GumbelPar(1,k) = alphaG;
    GumbelPar(2,k) = uG;
end
GumbelPar;

%% (3)Compute analytical Gumbel distributions

n = 0:0.5:110; %initiating a vector of 221 terms             
m = length(n); %storing n's length in one variable
GumbelCompute = zeros(m,6); %saving up some memory for the vector to be built
GumbelComputeMoments = zeros(m,6);

for k = 1:m
    for l = 1:6
       % using the gumbel method paramters 
       GumbelCompute(k,l) = exp(-exp(-GumbelPar(1,l)*(n(k)-GumbelPar(2,l))));
       % using the moments parameters
       GumbelComputeMoments(k,l) = exp(-exp(-GumbelParMoments(1,l)*(n(k)-GumbelParMoments(2,l))));
    end 
end 

%% (4) Plot 

figure(22)
colororder(newcolors)
plot(n, GumbelCompute) %using the Gumbel parameters as they are more precise
title('Gumbel Distributions') 
xlabel('Precipitation Depth h [mm]') 
ylabel('Empirical Frequencies [Fh] and Regression Curves')
hold on 
plot(sortedAnnualMax(), Fh, '.'); % (39x6), (39x1)
lgd = legend({' ',' ',' ',' ',' ',' ','1 hour','3 hours', ...
        '6 hours','12 hours', ...
        '24 hours','48 hours',}, ...
        'Location','southeast', 'NumColumns',2);
title(lgd, "Annual Max depth over time span [mm]")
print('GumbelDist','-vector','-dpdf') % this saves 'my_figure.pdf' (useful for LaTeX)
%% (5) Th and h

% computing return period

q = length(Fh);
Weibull_T = zeros(q,1); % saving space in the memory 

for k = 1:q
    Weibull_T(k) = 1/(1-Fh(k));
end 
Weibull_T;

%% (6) plotting 

% computing h by reverting analytical formula
hrevert = zeros(221,6);
T = 0:0.5:110;
for k = 1:221
    for l = 1:6
        hrevert(k,l) = GumbelPar(2,l) -log(-log(1-1/T(k)))/GumbelPar(1,l);
    end 
end 


% measured (dots) and estimated (smooth lines) rainfall depth vs return period 
figure(23)
n = 1:39;

colororder(newcolors)
plot(Weibull_T, sortedAnnualMax(),'o'); % dimensions : (39x1), (39x6) -> empirical data
title('Rainfall Depth vs Return Period')
xlabel('Return period T [years]') 
ylabel('Rainfall depth [mm]')
axis([0, 60, 0, 140])

hold on 
plot(T,hrevert); % (221x1), (221x6) % Gumbel distribution -> derived values
lgd = legend({' ',' ',' ',' ',' ',' ', '1 hour', '3 hours', '6 hours', ...
        '12 hours', '24 hours', '48 hours'}, ...
        "Location", "southeast", "NumColumns",2);
title(lgd, "Annual Max depth over time span [mm]")

hold off
print('RDvsT','-vector','-dpdf') % this saves 'my_figure.pdf' (useful for LaTeX)
%% (7) Building matrix of predicted depth

H_Gum = zeros(3,6); % initializing the desired matrix
k = 1;
for t = [10, 40, 100]
    s = T == t;
    for u = 1:6
        H_Gum(k,u) = hrevert(s,u);
    end 
    k = k + 1;
end 
H_Gum;

%% (8) saving elements 
T = [10 40 100];
D = [1 3 6 12 24 48];
save('assignment1_output_part2.mat','H_Gum', 'T', 'D');

