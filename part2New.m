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

%% Plotting Weibull coefficients 
figure(21)
stem(Fh)
title('Weibull coefficients')
xlabel('Rank') 
ylabel('Fh') 

%% computing reduced variable YF

YF = -log(-log(Fh));

% Plot of Fh vs h
% for first column
figure(22)
plot(sortedAnnualMax(), Fh, 'o');
title('Empirical Frequencies vs Precipitation Depth ')
xlabel('Precipitation Depth [mm]') 
ylabel('Empirical Frequencies [Fh]')
legend({'Annual Max depth in 1 hour [mm]','Annual Max depth in 3 hours [mm]', ...
      'Annual Max depth in 6 hours [mm]','Annual Max depth in 12 hours [mm]', ...
      'Annual Max depth in 24 hours [mm]','Annual Max depth in 48 hours [mm]',}, ...
      'Location','southeast')

%% (2) Fitting Gumbel curve
% method of moments
% calculating means and standard deviation for each duration 
std_vect = std(AnnualMax);
mean_vect = mean (AnnualMax);

%defining variables 
GumbelParMoments = zeros(2,6);
m = length(std_vect);

for k = 1:m
    alpha = pi/(std_vect(k)*sqrt(6));
    mu = mean_vect(k)-0.5772/alpha;

    GumbelParMoments(1,k) = alpha;
    GumbelParMoments(2,k) = mu;
end 
GumbelParMoments;
%%
% Gumbel method
GumbelPar = zeros(2,6);
m = length(std_vect);

mu_YF = mean(YF);
sigma_YF = std(YF);

% Defining paramaters 
for k = 1:m
    alphaG = sigma_YF/std_vect(k);
    muG = mean_vect(k)-(mu_YF/sigma_YF)*std_vect(k);
    
    GumbelPar(1,k) = alphaG;
    GumbelPar(2,k) = muG;
end
GumbelPar;

%% (3)Compute analytical Gumbel distributions

n = 0:0.5:110;%TO-DO : revoir si on veut pas partir à 1 plutôt...?
m = length(n);
GumbelCompute = zeros(m,6);
GumbelComputeMoments = zeros(m,6);

for k = 1:m
    for l = 1:6
       GumbelCompute(k,l) = exp(-exp(-GumbelPar(1,l)*(n(k)-GumbelPar(2,l))));
       GumbelComputeMoments(k,l) = exp(-exp(-GumbelParMoments(1,l)*(n(k)-GumbelParMoments(2,l))));
    end 
end 


%% (4) Plot 

%Going to set different call to "plot" to set the color right for each
%curve

figure(23)
plot(n, GumbelComputeMoments(:,1),"color",[0.6350 0.0780 0.1840])
hold on
plot(n, GumbelComputeMoments(:,2),"color",[0 0.4470 0.7410])
hold on
plot(n, GumbelComputeMoments(:,3),"color",[0.8500 0.3250 0.0980])
hold on
plot(n, GumbelComputeMoments(:,4),"color", [0.9290 0.6940 0.1250])
hold on 
plot(n, GumbelComputeMoments(:,5),"color",[0.4940 0.1840 0.5560])
hold on 
plot(n, GumbelComputeMoments(:,6),"Color",[0.4660 0.6740 0.1880])
title('Gumbel Distributions')
xlabel('Precipitation Depth h [mm]') 
ylabel('Empirical Frequencies [Fh] and Regression Values')
legend({'Annual Max depth in 1 hour [mm]','Annual Max depth in 3 hours [mm]', ...
        'Annual Max depth in 6 hours [mm]','Annual Max depth in 12 hours [mm]', ...
        'Annual Max depth in 24 hours [mm]','Annual Max depth in 48 hours [mm]',}, ...
        'Location','southeast')

hold on 
plot(sortedAnnualMax(), Fh, '.');
legend({'Annual Max depth in 1 hour [mm]','Annual Max depth in 3 hours [mm]', ...
        'Annual Max depth in 6 hours [mm]','Annual Max depth in 12 hours [mm]', ...
        'Annual Max depth in 24 hours [mm]','Annual Max depth in 48 hours [mm]',}, ...
        'Location','southeast')
% TO-DO : set the same colors to both plot elements 1 to 1 (regressions and
% points) i.e. the points corresponding to a regression curve should have
% the same color


%% (5) Th and h

%computing return period
q = length(Fh);
Weibull_T = zeros(q,1);

for k = 1:q
    Weibull_T(k) = 1/(1-Fh(k));
end 
Weibull_T;


% inverting Gumbel distribution 
T = zeros(221,6);
h = zeros(221,1);
for k = 1:221
    for l = 1:6
        T(k,l) = 1/(1-GumbelCompute(k,l));
        h(k,l) = muG - (1/alphaG)*log(-log(1-1/T(k,l)));
    end 
end 





figure(24)
n = 1:39;
m = 0:0.5:110;
plot(Weibull_T, sortedAnnualMax(),'o')

figure(25)
plot(T,h)

h







