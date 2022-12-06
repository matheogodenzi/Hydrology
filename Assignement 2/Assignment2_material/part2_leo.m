clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

% importing data
load('output_part1.mat');

% Setting a color code reference 
blue = [0 0.4470 0.7410];   % Used for precipitation values
orange = [0.8500 0.3250 0.0980];   % Used for catchment response
yellow = [0.9290 0.6940 0.1250];   % Used fo channel response 

%% (1)

n_sub = 4; %each hour is broken into n_int intervals
dt = 1/n_sub; % timestep [h]

t_Je= linspace(0,4-1,4); % precipitation time 
t_Jedt = dt*(dt:4/dt); % precipitation time for dt timesteps starting at dt (second version)
Jedt = interp1(t_Je, Je, t_Jedt ,'previous','extrap'); % effective precipitation extrapolation

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

% verification of gamma curve discretization 
sumIUHW = sum(IUHW*dt);

figure(1)
%bar(T*dt,IUHW); %first version
bar(t_iuh,IUHW, 'FaceColor', orange, 'EdgeColor', orange); %second version
xlabel('Timestep: '+ string(dt));
ylabel('Discretized probability distribution');
title('Watershed IUH');

%% Channel IUH

% defining IUH time interval 
cutoff = 70;
%TC = linspace(0,100, 100/dt); %first version --> less precise 
t_iuh=(dt:dt:cutoff); %second version 
% sz = size(TC);

%defining parameters of the probability density function 
L = 10^4;
D = 10^6;
c = 0.3*3600;

% first version 
%IUHC = L./sqrt(4*pi*D)*TC.^-(3/2).*exp(-((L-c*TC).^2)./(4*D*TC));
%IUHC(1)=0;

%second version 
IUHC = L./sqrt(4*pi*D)*t_iuh.^-(3/2).*exp(-((L-c*t_iuh).^2)./(4*D*t_iuh));

figure(3)
%bar(TC,IUHC);% first version
bar(t_iuh,IUHC, 'FaceColor', yellow, 'EdgeColor', yellow);% second version
xlabel('Timestep: '+ string(dt));
ylabel('Discretized probability distribution');
title('Channel IUH');

% verification of the inverse gaussian curve discretization 
sumIUHC = sum(IUHC*dt);

%% (3) IUHW figures 

figure(4)

bar(t_iuh, IUHW, 'FaceColor', orange);
hold on
bar(t_iuh, IUHC, 'FaceColor', yellow, 'FaceAlpha', 0.8);
xlabel('Timestep : ' + string(dt) +' h');
ylabel('Discretized probability distribution');
xlim([0 50]);
ylim([0 0.15]);
title('IUHW and IUHC');
legend("Watershed probability density discretization", "Channel probability density discretization");
%saveas(gcf,'figure3.png') %Saves the figure as "figure3" in png. format


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

% plotting the results 
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

% plotting the results to see the difference between the IUH and the system's response.
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
    bar(t_Jedt,Jedt(:,k), 'EdgeColor', blue, 'FaceColor', blue, 'LineWidth', 2);
    hold on 
    plot(xW, DischargeW(:,k), 'LineWidth', 1.5);
    hold on 
    plot(xC-dt/10, DischargeC(:,k), 'LineWidth', 1.5); %-dt/10 is a way of visualizing the data. This way the histograms are not exactly superposed
    hold off
    xlim([0 35]);
    ylim([0 max(Jedt(:,k)) + 0.2]);
    xlabel("time [h]")
    ylabel("intensity [mm/h]")
    if k == 1
        legend("Effective Precipitation Intensity (Je)",  "Watershed Discharge Intensity (QW)",  "Channel Discharge Intensity (QC)")
    end 
    title(" Precipitation Event & Hydrological Response - " + string(k))
end 
%saveas(gcf,'figure6.png') %Saves the figure as "figure6" in png. format

%% (6) Plotting of the effective precipitation intensity with catchment and chanel discharge

% Méthode qui permet d'afficher les precipitations effectives en inversé
% mais sur des figures différentes pour chaque event

for k = 1:ev_nbr
    fig = figure;
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    
    xW = (dt:dt:NW*dt);
    xC = (dt:dt:NC*dt);
    
    yyaxis left 
    bar(t_Jedt, Jedt(:,k), 'FaceColor', blue, 'EdgeColor', blue, 'LineWidth', 2);
    axis ([0 50 0 max(Jedt(:,k))+2]);
    set(gca, 'Ydir', 'reverse');
    ylabel('Precipitation intetnsity (Je) [mm/h]')
    
    hold on 
    yyaxis right 
    stairs (xW, DischargeW(:,k),'LineWidth', 1.5, 'LineStyle', '-', 'Color', orange);
    stairs (xC-dt/10, DischargeC(:,k),'LineWidth', 1.5, 'LineStyle', '-', 'Color', yellow);
    ylabel('Discharge (Q) [mm/h]');
    axis ([0 50 0 max(DischargeW(:,k))+1]);
    hold off 
    legend ('Effective precipitation intensity (Je)', 'Watershed discharge (Qw)', 'Channel discharge (Qc)');
    xlabel('Time [h]');
    title('Precipitation & Hydrological Response of Event ' +string(ev_nbr));
    %saveas(gcf,'figure6.' + str(k)'.png') %Saves the figure as "figure6.k" in png. format
end


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
save('output_part2.mat','dt','IUHW');
