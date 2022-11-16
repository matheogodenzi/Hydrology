% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 4: Finding the the 16h window statistics
% -------------------------------------------------------------------------

clear variables %clear the workspace variables
close all %close alla figures
clc %clear the command window

% -------------------------------------------------------------------------
% # 1-3: Data import and cleaning
% useful functions: readtable, isnan, year, month
% -------------------------------------------------------------------------

% import data into table T
T = readtable('data.txt',...
    'HeaderLines', 2,...
    'Format','%s%s%f',... %the format is: text string, text string and float number
    'TreatAsEmpty','-'); %this is how empty data is reported (see legend)

% Create a vector h containing the hourly precipitation and a vector t
% containing the timestamp
h = T.rre150h0;    %rre150h0 is the MeteoSwiss code for hourly rainfall depth [mm]
t = datetime(T.time,'InputFormat','yyyyMMddHH'); %convert to datetime (can be slow)
m = month(t); %gives a value in 1-12 to indicate the month of each date
y = year(t); %gives the year of each date

% fix empty values (which appear as NaN values in the Matlab)
emptyValues = isnan(h); %logical test to tell whether a value is missing or not
h(emptyValues) = 0; %give zero to those values
fprintf('%i empty values\n', sum(emptyValues)); %display how many missing values there are


% -------------------------------------------------------------------------
% # Compute rainfall maxima for 16 consecutive hours over a year
% -------------------------------------------------------------------------
%%
sixteen_hours_max = zeros(39,1); % building a vector 

count = 1;
for k = 1981:2019 %iterationg over the years
    indices = y == k; %building booleen matrix to get a hold of the indiced...
    % corresponding to the right year in the matrix of precipitation
    yearly_prec = h(indices);
    max = 0;
    for l = 1:sum(indices)-15 % summing over the window
        s = sum(yearly_prec(l:l+15));
        if s >= max
            max = s;
        end
    end 
    sixteen_hours_max(count,1) = max;
    count = count + 1; 
end 

sixteen_hours_max;

% sorting data columns by amplitude
SortedSixteenMax = sort(sixteen_hours_max,1, 'ascend');
SortedSixteenMax;

%% Approximating the return period based on empirical values 
T_empirical = (length(SortedSixteenMax)+1)/((length(SortedSixteenMax)+1)-(sum(SortedSixteenMax<80)+1));
%here we set a bondary where rainfall depths is strictly smaller than 80mm.

%it is not the most precise because there's a gap in the data between 76 mm
%and 90 mm, Thereofre we eill proceed to a gumbel fit to obtain a closer
%estimation of the actual return period. 

%% Computing Weibull coefficients 

% initializing empty vector
Fh = zeros([size(SortedSixteenMax, 1),1]);

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
figure(41)
colororder(newcolors)
plot(SortedSixteenMax, Fh, '.');
title('Empirical Frequencies vs Precipitation Depth ')
xlabel('Precipitation Depth [mm]') 
ylabel('Empirical Frequencies [Fh]')
lgd = legend({'16 hour'}, ...
        "Location", "southeast", "NumColumns",2);
title(lgd, "Annual Max depth over time span [mm]")
%print('FhvsD','-vector','-dpdf') % this saves 'my_figure.pdf' (useful for LaTeX)

%% (2) Fitting Gumbel curve
% Gumbel method
std_vect = std(SortedSixteenMax);
mean_vect = mean(SortedSixteenMax);
GumbelPar = zeros(2,1);
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

n = 0:0.2:110; %initiating a vector of 551 terms             
m = length(n); %storing n's length in one variable
GumbelCompute16 = zeros(m,1); %saving up some memory for the vector to be built

for k = 1:m
   % using the gumbel method paramters 
   GumbelCompute16(k,1) = exp(-exp(-GumbelPar(1,1)*(n(k)-GumbelPar(2,1))));
end 
GumbelCompute16;
%% (4) Plot 

figure(42)
colororder(newcolors);
plot(n, GumbelCompute16) %using the Gumbel parameters as they are more precise
title('Gumbel Distributions');
xlabel('Precipitation Depth h [mm]') ;
ylabel('Empirical Frequencies [Fh] and Regression Curves');
hold on 
plot(SortedSixteenMax, Fh, '.'); % (39x1), (39x1)
lgd = legend({'regression ','empirical data',}, ...
        'Location','southeast', 'NumColumns',2);
title(lgd, "Annual Max depth over 16h [mm]");
%print('GumbelDist','-vector','-dpdf') % this saves 'my_figure.pdf' (useful for LaTeX)
%% (5) Computing the return period based on empirical data

%computing return period

q = length(Fh);
Weibull_T = zeros(q,1); %saving some space in the memory 

for k = 1:q
    Weibull_T(k) = 1/(1-Fh(k));
end 
Weibull_T;

%% (6) Computing return period based on gumbel parameters and plotting

%Here's the answer to the question, derived from the Gumbel distribution
%fit.
T16 = 1/(1-exp(-exp(-GumbelPar(1)*(80-GumbelPar(2)))))

% computing h by reverting analytical formula
hrevert = zeros(551,1);
T = 0:0.2:110;
for k = 1:551
    hrevert(k,1) = GumbelPar(2,1) -log(-log(1-1/T(k)))/GumbelPar(1,1);
end 


% measured (dots) and estimated (smooth lines) rainfall depth vs return period 
figure(43)
n = 1:39;

colororder(newcolors)
plot(Weibull_T, SortedSixteenMax,'o'); % dimensions : (39x1), (39x6) -> empirical data
title('Rainfall Depth vs Return Period')
xlabel('Return period T [years]') 
ylabel('Rainfall depth [mm]')
axis([0, 60, 0, 140])

hold on 
plot(T,hrevert); % (221x1), (221x6) % Gumbel distribution -> derived values
lgd = legend({'empirical data', 'regression'}, ...
        "Location", "southeast", "NumColumns",2);
title(lgd, "Annual Max depth over 16 hours [mm]")

hold off
%print('RDvsT','-vector','-dpdf') % this saves 'my_figure.pdf' (useful for LaTeX)