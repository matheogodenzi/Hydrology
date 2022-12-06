clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

% importing data

%% (1) Prepare precipitation data
load('output_part2.mat');

T = readtable('rainfall_data.txt',...
    'HeaderLines', 2,...
    'Format','%s%s%f',... %the format is: text string, text string and...
    ...                   % float number...
    'TreatAsEmpty','-'); %this is how empty data is reported (see legend)

% Create a vector h containing the hourly precipitation and a vector t
% containing the timestamp
h = T.rre150h0; %rre150h0 is the MeteoSwiss code for hourly rainfall...
                %depth [mm]
t = datetime(T.time,'InputFormat','yyyyMMddHH'); %convert to datetime...
                                                 %(can be slow)
m = month(t); %gives a value in 1-12 to indicate the month of each date
y = year(t); %gives the year of each date

% fix empty values (which appear as NaN values in the Matlab)
emptyValues = isnan(h); %logical test to tell whether a value is missing...
                        % or not
h(emptyValues) = 0; %give zero to those values
fprintf('%i empty values\n', sum(emptyValues)); %display how many...
                                                %missing values there are
%%                                                
% defining sub intervals
[nr,nc]=size(T);
tprime= linspace(0,nr-1,nr);
t_dt = dt*(dt:nr/dt); % precipitation time for dt timesteps starting at dt (second version)
Pdt = interp1(tprime, h, t_dt ,'previous','extrap'); % effective precipitation extrapolation

%% Verifying sub intervals 
figure
new_t = t(1):minutes(dt*60):t(end);
length(new_t)
length(Pdt)
bar(new_t,Pdt(1:length(new_t)))
ylabel('rainfall intensity [mm/h] over timestep dt')
title('total precipitation throughout the year 2018')

%% (2) Separation into effective precipitation and infiltration
Je = 0.3*Pdt;
I = 0.1*Pdt;

%% (3) Surface contributions

% Convolution over the effective precipitation falling on the watershed
n = length(Je);
M = length(IUHW);
NW = n+M-1;
DischargeW = zeros(1,NW);

for i = 1:n
    DischargeW(i:i+M-1) = DischargeW(i:i+M-1)+ Je(i)*IUHW*dt;
end 

DischargeW = transpose(DischargeW);

%% Conv matlab - verification
%intended to verify the handly writeen convolution

QW = conv(Je,IUHW)*dt;

figure
subplot(3,1,1);
bar(QW);
subplot(3,1,2); 
bar(DischargeW);
subplot(3,1,3);
bar(transpose(QW)-DischargeW)

% my convolution seems appropriate

%diff = transpose(QW)-DischargeW

%% (4.1) Subsurface contributions 

load('../Parameters/IUHpars.mat')
par_scale = par_scale*10;

% defining IUH time interval 
cutoff = 500;
t_iuh=(dt:dt:cutoff); % time vector

%IUHW = gampdf(T , par_shape, par_scale);% first version --> less precise 
IUHSW = gampdf(t_iuh , par_shape, par_scale); % IUH sub surface generation 

% verification of gamma curve 
sumIUHSW = sum(IUHSW*dt);

figure
%bar(T*dt,IUHW); %first version
bar(t_iuh*dt,IUHSW); %second version

%% (4.2) subsurface convolution 

% Convolution over the effective precipitation falling on the watershed
n = length(I);
M = length(IUHSW);
NSW = n+M-1;
DischargeSW = zeros(1,NSW);

for i = 1:n
    DischargeSW(i:i+M-1) = DischargeSW(i:i+M-1)+ I(i)*IUHSW*dt;
end 

DischargeSW = transpose(DischargeSW);

%% Conv matlab - verification 
%intended to verify the handly writeen convolution

QSW = conv(I,IUHSW)*dt;

figure
subplot(3,1,1);
bar(QSW);
subplot(3,1,2);

bar(DischargeSW);
subplot(3,1,3);
bar(transpose(QSW)-DischargeSW)

% my convilution seems appropriate

%diff = transpose(QW)-DischargeW


%% (5) Total discharge 

sz = length(t_dt);

Qtot = DischargeW(1:sz, 1)+DischargeSW(1:sz, 1);

%% (6) Figure

% extrapolating month values to the same 
t_dt = dt*(dt:nr/dt); % precipitation time for dt timesteps starting at dt (second version)
mdt = interp1(tprime, m, t_dt ,'previous','extrap'); % effective precipitation extrapolation

% creating relevant november vector with sub interval time
t1 = datetime('01/11/2018 00:00:00','InputFormat','dd/MM/yyyy HH:mm:ss');
t2 = datetime('30/11/2018 23:59:59','InputFormat','dd/MM/yyyy HH:mm:ss');
tnov = t1:minutes(dt*60):t2;

% collecting all timesteps relative to the month of november
november = mdt==11;

figure
subplot(2,1,1)
plot(tnov, Qtot(november), LineWidth=1.5)
hold on 
bar(tnov, DischargeW(november), 'FaceAlpha', .8)
hold on
bar(tnov, DischargeSW(november), 'FaceAlpha', .8)
hold off
ylabel('discharge [mm/h]')
legend("total runoff", "surface runoff", "subsurface runoff")
title('Distinct Runoffs')
subplot(2,1,2)
plot(tnov, DischargeW(november)./DischargeSW(november))
ylabel('percentage [%]')
legend("QW/QSW ratio")
title('November Time Series Ratio')

%% (7) Compute total amount of timsteps 

high = sum(Qtot > 0.2); % boolean mask
low = sum(Qtot < 0.002);

high_percentage = high/length(Qtot)*100;
low_percentage = low/length(Qtot)*100;


%% Bonus - additional graphs 

ti = t(1);
tf = t(length(t));

%defining sub intervals time vector 
tyear = ti:minutes(dt*60):tf;

% plotting Qtot throughout the year
figure 
plot(tyear, Qtot(1:length(tyear)))

% plotting QW and QSW throughout the year 
figure
subplot(2,1,1)
plot(tyear, DischargeW(1:length(tyear), 1))
subplot(2,1,2)
plot(tyear, DischargeSW(1:length(tyear), 1))

%%
figure 
plot(tyear, Qtot(1:length(tyear)))
ylabel('Discharge intensity [mm/h] over time step dt')
title('Total Discharge Throughout the Year')