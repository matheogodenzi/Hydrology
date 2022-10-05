% -------------------------------------------------------------------------
% Course: Hydrology for Engineers
% Assignment 1
% Part 1: Process rainfall data from MeteoSwiss
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
% # 4: Plot with annual rainfall over the years
% -------------------------------------------------------------------------
annual_prec = zeros(39,1);

count = 1;
for k = 1981:2019
    indices = y == k;
    annual_prec(count, 1) = sum(h(indices));
    count = count + 1;
end 

annual_prec;
years = 1981:2019;


figure
stem(years, annual_prec,'Color',[0,0.7,0.9]);

title('Annual Precipitations');
xlabel('years');
ylabel('precipitations (mm)');

% -------------------------------------------------------------------------
% # 5: Compute rainfall maxima of a certain duration
% -------------------------------------------------------------------------

max = 0;
for k = 1:length(h)-2
    s = sum(h(k:k+2));
    if s >= max
        max = s;
    end
end 

% -------------------------------------------------------------------------
% # 6: Compute rainfall maxima of a certain duration over a year
% -------------------------------------------------------------------------
%%
three_hours_max = zeros(39,1);

count = 1;
for k = 1981:2019
    indices = y == k;
    yearly_prec = h(indices);
    max = 0;
    for l = 1:sum(indices)-2
        s = sum(yearly_prec(l:l+2));
        if s >= max
            max = s;
        end
    end 
    three_hours_max(count,1) = max;
    count = count + 1;
end 

three_hours_max;

years = 1981:2019;


figure
stem(years, three_hours_max,'Color',[0,0.7,0.9]);

title('Three hours max -- Yearly');
xlabel('years');
ylabel('precipitations (mm)');

% -------------------------------------------------------------------------
% # 7: Extend to multiple years and multiple durations
% -------------------------------------------------------------------------
%%

AnnualMax = zeros(39,6);

count = 1;
for k = 1981:2019
    indices = y == k;
    yearly_prec = h(indices);

    max = 0;
    for l = 1:sum(indices)
        s = yearly_prec(l);
        if s >= max
            max = s;
        end
    end 
    AnnualMax(count,1) = max;
    
    
    max = 0;
    for l = 1:sum(indices)-2
        s = sum(yearly_prec(l:l+2));
        if s >= max
            max = s;
        end
    end 
    AnnualMax(count,2) = max;
    
    max = 0;
    for l = 1:sum(indices)-5
        s = sum(yearly_prec(l:l+5));
        if s >= max
            max = s;
        end
    end
    AnnualMax(count,3) = max;
    
    max = 0;
    for l = 1:sum(indices)-11
        s = sum(yearly_prec(l:l+11));
        if s >= max
            max = s;
        end
    end
    AnnualMax(count,4) = max;

    max = 0;
    for l = 1:sum(indices)-23
        s = sum(yearly_prec(l:l+23));
        if s >= max
            max = s;
        end
    end
    AnnualMax(count,5) = max;
    
    max = 0;
    for l = 1:sum(indices)-47
        s = sum(yearly_prec(l:l+47));
        if s >= max
            max = s;
        end
    end
    AnnualMax(count,6) = max;
    
    count = count + 1;
end 

AnnualMax


% -------------------------------------------------------------------------
% # 7 bis: Extend to multiple years and multiple durations
% -------------------------------------------------------------------------
%%

AnnualMax = zeros(39,6);

count = 1;
for k = 1981:2019
    indices = y == k;
    yearly_prec = h(indices);
    windows = [1 3 6 12 24 48];
    
    ind = 1;
    for m = 1:49
        if ismember(m, windows)
            max = 0;
            for l = 1:sum(indices)+1-m
                s = sum(yearly_prec(l:l+m-1));
                if s >= max
                    max = s;
                end
            end 
            AnnualMax(count,ind) = max;
            ind = ind + 1;
        end
    end 
    count = count + 1;
end 

AnnualMax

%% Meme resultat, je trouve plus simple cette ecriture
%{
D = [1,3,6,12,24,48];

AnnualMax = zeros(39,6);


    count = 1;
    for k = 1981:2019
        indices = y == k;
        yearly_prec = h(indices);
        for f = 1:length(D)
            d = D(f)-1;
            max = 0;
            for l = 1:sum(indices)-d
                s = sum(yearly_prec(l:l+d));
                if s >= max
                    max = s;
                end
            end 
            AnnualMax(count,f) = max;
 
        end 
  count = count + 1;
    
    end
AnnualMax



%}





% -------------------------------------------------------------------------
% # 8: save the output
% useful functions: save
% -------------------------------------------------------------------------
% when you are confident about your results, save the variables AnnualMax
% and D by uncommenting the following line: 
save('assignment1_output_part1.mat','D','H_Gum');
