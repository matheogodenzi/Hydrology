
clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

% importing data

load('output_part1.mat');
%% (1)

n_sub = 5; %each hour is broken into n_int intervals
dt = 1/n_sub; % timestep [h]

t_Je= linspace(0,4-1,4); % precipitation time 
%t_Jedt = linspace(0,4-dt,4/dt); % precipitation time for dt timesteps starting at zero (first version)
t_Jedt = dt*(dt:4/dt); % precipitation time for dt timesteps starting at dt (second version)
Jedt = interp1(t_Je, Je, t_Jedt ,'previous','extrap'); % effective precipitation extrapolation
Jedt
%% (2) Watershed IUH
% loading IUH params
load('../Parameters/IUHpars.mat')
% generating watershed IUH gamma distribution
%gammaIUH = zeros(size(Jedt));

% defining IUH time interval 
cutoff = 70;
%T = linspace(0,1000, 1000/dt); % first version
t_iuh=(dt:dt:cutoff); %second version 

%IUHW = gampdf(T , par_shape, par_scale);% first version --> less precise 
IUHW = gampdf(t_iuh , par_shape, par_scale); %second version

% verification of gamma curve 
sumIUHW = sum(IUHW*dt)

figure(1)
%bar(T*dt,IUHW); %first version
bar(t_iuh*dt,IUHW); %second version

%% Channel IUH

% defining IUH time interval 
cutoff = 70;
%TC = linspace(0,100, 100/dt); %first version --> less precise 
t_iuh=(dt:dt:cutoff); %second version 
% sz = size(TC);

% defining parameters of the probability density function 
L = 10^4; % [m]
D = 10^6; % [m2/h]
c = 0.3*3600; % [m/h] (3600 is to convert seconds to hours)

% first version 
%IUHC = L./sqrt(4*pi*D)*TC.^-(3/2).*exp(-((L-c*TC).^2)./(4*D*TC));
%IUHC(1)=0;

%second version 
IUHC = L./sqrt(4*pi*D)*t_iuh.^-(3/2).*exp(-((L-c*t_iuh).^2)./(4*D*t_iuh));

figure(3)
%bar(TC,IUHC);% first version
bar(t_iuh,IUHC);% second version
sumIUHC = sum(IUHC*dt)

%% (3) IUHW figures 

figure(4)

subplot(2,1,1);
bar(t_iuh, IUHW, 'red');
xlabel('step : ' + string(dt) +' h');
ylabel('discretized probability distribution');
xlim([0 50]);
ylim([0 0.15]);
title('IUHW');
legend("Watershed Probability Density Curve");

subplot(2,1,2);
bar(t_iuh, IUHC);
xlabel('step : ' + string(dt) +' h');
ylabel('discretized probability distribution');
xlim([0 50]);
ylim([0 0.15]);
title('IUHC');
legend("Channel Probability Density Curve");

%% (4) Convolution 

% Convolution over the effective precipitation falling on the watershed
[nrowsW,ncolsW] = size(Jedt);
M = length(IUHW);
NW = nrowsW+M-1;
DischargeW = zeros(3,NW);

for l = 1:ncolsW
    for i = 1:nrowsW
        DischargeW(l,i:i+M-1) = DischargeW(l,i:i+M-1)+ Jedt(i,l)*IUHW*dt;
    end 
end 
DischargeW = transpose(DischargeW);

%% matlab conv

QW = conv(IUHW, Jedt(:, 3)*dt);

figure
subplot(2,1,1)
bar(DischargeW(:,3))
subplot(2,1,2)
bar(QW)

%% plotting the results 
xW = (dt:dt:NW*dt);
figure
for i = 1:3
    subplot(3,1,i);
    bar(xW,DischargeW(:, i))
    title("event" + string(i))
    xlim([0 40])
    ylim([0 1.1])
end 

%% (5) IUHC convolution 

[nrowsC, ncolsC] = size(DischargeW);
M = length(IUHC);
NC = nrowsC+M-1;
DischargeC = zeros(3,NC);
 
for l = 1:ncolsC % the inverse of what was done in the previous convolution
    for i = 1:nrowsC % see comment above. This is due to the shape of the matrix defined 
        DischargeC(l,i:i+M-1) = DischargeC(l,i:i+M-1)+ DischargeW(i,l)*IUHC*dt;
    end 
end 
DischargeC = transpose(DischargeC);
%%
QC = conv(IUHC, QW*dt);

figure
subplot(2,1,1)
bar(DischargeC(:,3))
subplot(2,1,2)
bar(QC)
%% plotting the results to see the difference between the IUH and the system's response.
figure
plot(DischargeW(:,1))
xlim([0 500])
ylim([0 0.5]);
hold on 
plot(DischargeC(:,1));
xlim([0 500]);
ylim([0 0.5]);
hold off

%% (6) figure containing everything 
[spl_nbr, ev_nbr] = size(Jedt);

figure
for k = 1:ev_nbr
    subplot(ev_nbr, 1, k)
    xW = (dt:dt:NW*dt);
    xC = (dt:dt:NC*dt);
    bar(t_Jedt,Jedt(:,k));
    hold on 
    bar(xW, DischargeW(:,k));
    hold on 
    bar(xC-dt/10, DischargeC(:,k)); %-dt/10 is a way of visualizing the data. This way the histograms are not exactly superposed
    hold off
    xlim([0 35]);
    ylim([0 max(Jedt(:,k)+0.2)]);
    xlabel("time [h]")
    ylabel("intensity [mm/h]")
    if k == 1
        legend("Effective Precipitation Intensity (Je)",  "Watershed Discharge Intensity (QW)",  "Channel Discharge Intensity (QC)")
    end 
    title(" Precipitation Event & Hydrological Response - " + string(k))
end 

%% crash test 

figure
xW = (dt:dt:NW*dt);
xC = (dt:dt:NC*dt);
ax = axes;
set(gca, 'YDir', 'reverse');
bar(t_Jedt,Jedt(:,1));
ax = axes;
set(gca, 'Color', 'none')
hold on
plot(xW, DischargeW(:,1));
hold on 
plot(xC-dt/10, DischargeC(:,1)); %-dt/10 is a way of visualizing the data. This way the histograms are not exactly superposed
hold off
xlim([0 45]);
ylim([0 0.8]);
xlabel("time [h]")
ylabel("intensity [mm/h]")
legend("Effective Precipitation Intensity (Je)",  "Watershed Discharge Intensity (QW)",  "Channel Discharge Intensity (QC)")
title("Precipitation & Hydrological Response")
%% (7) Discharge peaks

% getting max discharge 
[nr, nc] = size(DischargeW);
peakW = zeros(1,nc);
pkW_id = zeros(1,nc);
peakC = zeros(1,nc);
pkC_id = zeros(1,nc);
for i = 1:nc
    [peakW(i), pkW_id(i)] = max(DischargeW(:,i));
    [peakC(i), pkC_id(i)] = max(DischargeC(:,i));
end

% getting corssponding time in hours 
max_time_W = zeros(1, nc);
max_time_C = zeros(1, nc);
for i = 1:nc
    max_time_W(i) = pkW_id(i)*dt; % the index represents the subinterval in terms of time. By multiplying by dt we get the ectual time in hours
    max_time_C(i) = pkC_id(i)*dt;
end
%% (8) saving variables 
%uncomment the next line to save the desired results on your terminal
%save('output_part2.mat','dt','IUHW');

