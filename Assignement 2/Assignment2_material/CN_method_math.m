% Specific goal: Compute effective precipitation Je and infiltration I intensity  
% for three different rainfall events using the Curve Number method


% beginning of the script

clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

%%
% --------------------------------------------------------------------------
% 1) Data import
%--------------------------------------------------------------------------


% importing the data

load('event_1.mat')
load('event_2.mat')

% for T = 100 our DDF had the following parameters: 
c = 40.2685;
e = 0.7035;
f = -0.0452;

% for D = 4h it yields the following
D = 4;
h_4h = c*D/(D.^e+ f);

% assigning mean intensity to each 1h step during 4h event
event_3 = zeros(4,1);
for k = 1:4
    event_3(k)=h_4h/4;
end 

event_3;

% creating a matrix 
event_matrix = [event_1(:), event_2(:), event_3(:)];


%%
% --------------------------------------------------------------------------
% 2) Compute the average CN of the catchment
%--------------------------------------------------------------------------
load('../Parameters/SCSpars.mat')
percentages = reshape(percentages,[1,3]);
% weighted average
Average_CN = sum(CN.*percentages)/sum(percentages);


%%
% --------------------------------------------------------------------------
% 3) Implement the CN method
%--------------------------------------------------------------------------

% Potential maximum soil moisture retention S
S = 25.4*(1000/Average_CN-10);

%SETTING VARIABLE VECTORS FOR EACH PARAMETERS 
% The cumulative precipitation (P [mm])
P = zeros(4,3);
% The initial abstraction (Ia [mm])
Ia = zeros(1,3);
% The cumulative infiltration (Fa [mm]);
Fa = zeros(4,3);
% The cumulative effective precipitation (Pe [mm]);
Pe = zeros(4,3);
% The infiltration intensity (I [mm/h]);
I = zeros(4,3);
% The effective rainfall intensity (Je [mm/h]);
Je = zeros(4,3);

for l = 1:3
    Ia(l)=0.2*S;
    for k = 1:4
        % time step of an hour so we multiply each rainfall intensity by 1 and
        % add it to the previous to get the cumulative rainfall after n hours
        if k == 1
            P(k,l)= event_matrix(k,l);
        else 
            P(k,l)= P(k-1, l) + event_matrix(k,l);
        end 
        
        % computing Pe
        if P(k,l) <= Ia(l)
            Pe(k,l)=0;
        else 
            Pe(k,l) = ((P(k,l)-Ia(l))^2)/(P(k,l)-Ia(l)+S);
        end 

        % computing Fa
        if P(k,l) <= Ia(l)
            Fa(k,l) = 0;
        else 
            Fa(k,l)= (S*Pe(k,l))/(P(k,l)-Ia(l));
        end 

        % computing Je = dPe/dt and I since the time step in 1 we divide by 1
        if k == 1
            Je(k,l) = Pe(k,l);
            I(k,l) = Fa(k,l);
        else
            Je(k,l)= (Pe(k,l)-Pe(k-1,l))/1;
            I(k,l)= (Fa(k,l)-Fa(k-1,l))/1;
        end 
    end 
end 


%%
% --------------------------------------------------------------------------
% 4) Show rainfall separation in a figure
%--------------------------------------------------------------------------

% Make a figure that includes, for each event: the rainfall intensity J [mm/h],
% the infiltration intensity I [mm/h] and the effective precipitation intensity 
% Je [mm/h]. Please use the same limits for the y-axis in the three events. 
% You can either use one figure and three subplots or three different figures 
% for the three events. Line plots with markers or bar plots are equally fine.
t = 1:4;
figure(1)
legend("Precipitation Intensity (J)",  "Infiltration Intensity (I)",  "Effective Precipitation Intensity (Je)");
for m = 1:3
    subplot(3,1,m)
    stem(t, event_matrix(:,m));
    hold on
    plot(t, I(:,m), '-o');
    hold on 
    bar(t, Je(:,m));
    hold off
    xlabel('time [h]')
    ylabel('intensity [mm/h]')
    ylim([0 20])
    title('event ' + string(m))
    legend("J",  "I",  "Je");
end 

%% This is the plot version we kept for the report

figure (2)
for m= 1:3
    subplot(3,1,m)
    bar(t, [event_matrix(:, m) I(:,m) Je(:,m)])
    xlabel('time [h]')
    ylabel('intensity [mm/h]')
    ylim([0 20])
    title('event ' + string(m))
    if m == 1
        legend("Precipitation Intensity (J)",  "Infiltration Intensity (I)",  "Effective Precipitation Intensity (Je)")
    end
end 

% les docs sur le drive ont tous fait full bar plot, ça rend assez bien

%%
% --------------------------------------------------------------------------
% 5) Compute
%--------------------------------------------------------------------------

% since they are all relative to one hour we multiply intensities by 1
Je_tot = sum(1*Je, 1);

%%
% --------------------------------------------------------------------------
% 6) saving result
%--------------------------------------------------------------------------
J = event_matrix;
save('output_part1.mat','J','Je', 'I');