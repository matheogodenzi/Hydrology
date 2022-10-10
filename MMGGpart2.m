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
plot(sortedAnnualMax(), Fh, '.');
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
mean_vect = mean(AnnualMax);

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

n = 0:0.5:110;%TO-DO : revoir si on veut pas partir à 1 plutôt...?  % selon les graphe du cours, mieux de partir à 0 MMG               
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

figure(23)
plot(n, GumbelComputeMoments) % (1x39), (221x6) ??
title('Gumbel Distributions')
xlabel('Precipitation Depth h [mm]') 
ylabel('Empirical Frequencies [Fh] and Regression Values')
legend({'Annual Max depth in 1 hour [mm]','Annual Max depth in 3 hours [mm]', ...
        'Annual Max depth in 6 hours [mm]','Annual Max depth in 12 hours [mm]', ...
        'Annual Max depth in 24 hours [mm]','Annual Max depth in 48 hours [mm]',}, ...
        'Location','southeast')

hold on 
plot(sortedAnnualMax(), Fh, '.'); % (39x6), (39x1)
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


T = zeros(221,6);
for k = 1:221
    for l = 1:6
        T(k,l) = 1/(1-GumbelCompute(k,l));
    end 
end 

%% (6) plotting 

figure(24)

n = 1:1:39;

plot(Weibull_T, sortedAnnualMax(),'.'); % dimensions : (39x1), (39x6)

figure(25)
n = 1:39;
plot(Weibull_T, sortedAnnualMax(),'.'); % dimensions : (39x1), (39x6)
hold on 
% compute h by reverting analytical formula (check variables)
% hrevert = muG - log(-log(1-1/T))/alphaG;
hrevert = zeros(221,6);
Tbis = 0:0.5:110;
for k = 1:221
    for l = 1:6
        hrevert(k,l) = GumbelPar(2,l) -log(-log(1-1/Tbis(k)))/GumbelPar(1,l);
    end 
end 

plot(Tbis,hrevert) % (221x1), (221x6)


%% (7) Building matrix of predicted depth

H_Gum = zeros(3,6);
k = 1;
for t = [10, 40, 100]
    s = Tbis == t;
    for u = 1:6
        H_Gum(k,u) = hrevert(s,u);
    end 
    k = k + 1;
end 
H_Gum

%% (8) saving elements 
T = [10 40 100];
D = [1 3 6 12 24 48];
save('assignment1_output_part2.mat','H_Gum', 'T', 'D');

