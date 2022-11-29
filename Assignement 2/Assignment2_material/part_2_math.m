
clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

% importing data

load('output_part1.mat');
%% (1)

n_sub = 3; %each hour is broken into n_int intervals
dt = 1/n_sub; % timestep [h]

t_Je= linspace(0,4-1,4);
t_Jedt = linspace(0,4-dt,4/dt);
Jedt = interp1(t_Je, Je, t_Jedt ,'previous','extrap');

%% (2) Watershed IUH
% loading IUH params
load('../Parameters/IUHpars.mat')
% generating watershed IUH gamma distribution
%gammaIUH = zeros(size(Jedt));

% trying on a finite time interval
T = linspace(0,1000, 1000/dt);

IUHW = gampdf(T , par_shape, par_scale);

% verification of gamma curve 
sumIUHW = sum(IUHW*dt);

figure(1)
bar(T*dt,IUHW);

%% testing gammy function

T_test = 0:40;
gammaIUH_test = gampdf(T_test ,7, 1);

figure(2)
bar(gammaIUH_test)

%% Channel IUH

TC = linspace(0,100, 100/dt);
sz = size(TC);
L = 10^4;
D = 10^6;
c = 0.3*3600;
%length((L-c*TC).^2)
%length(4*D*T)
%length(sqrt(4*pi*D)*T.^(3/2))

IUHC = L./sqrt(4*pi*D)*TC.^-(3/2).*exp(-((L-c*TC).^2)./(4*D*TC));
IUHC(1)=0;

figure(3)
bar(TC,IUHC);
sumIUHC = sum(IUHC*dt);

%% (3) IUHW figures 

figure(4)

subplot(2,1,1);
bar(T, IUHW, 'red');
xlabel('step : ' + string(dt) +' h');
ylabel('discretized probability distribution');
xlim([0 150]);
ylim([0 0.2]);
title('IUHW');
legend("Probability Density Curve");

subplot(2,1,2);
bar(TC, IUHC);
xlabel('step : ' + string(dt) +' h');
ylabel('discretized probability distribution');
xlim([0 150]);
ylim([0 0.2]);
title('IUHC');
legend("Probability Density Curve");

%% (4) Convolution 

% Convolution over the effective precipitation falling on the watershed
n = length(Jedt);
M = length(IUHW);
N = n+M-1;
DischargeW = zeros(1,N);

for i = 1:n
    DischargeW(i:i+M-1) = DischargeW(i:i+M-1)+ Je(i)*IUHW*dt;
end 
DischargeW;

%% plotting the results to see the difference between the IUH and the system's response.
figure
subplot(2,1,1);
bar(DischargeW)
xlim([0 500])
subplot(2,1,2);
bar(IUHW);
xlim([0 500]);

%% (5) IUHC convolution 

n = length(Jedt);
M = length(IUHC);
N = n+M-1;
DischargeC = zeros(1,N);

for i = 1:n
    DischargeC(i:i+M-1) = DischargeC(i:i+M-1)+ Je(i)*IUHC*dt;
end 
DischargeW;
%% plotting the results to see the difference between the IUH and the system's response.
figure
subplot(2,1,1);
bar(DischargeC)
xlim([0 500])
ylim([0 0.6]);
subplot(2,1,2);
bar(IUHC);
xlim([0 500]);
ylim([0 0.6]);