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


% -------------------------------------------------------------------------
% # 4: Plot with annual rainfall over the years
% -------------------------------------------------------------------------
annual_prec = zeros(39,1);   %Initiating an empty vector to store the...
                             %total annual rainfall for each years
count = 1;
for k = 1981:2019 %To iterate over all years
    indices = y == k; %Logical indexing to select the measurements of...
                      %the indexed year
    annual_prec(count, 1) = sum(h(indices)); %Sum of those measurements
    count = count + 1;
end 

annual_prec;
years = 1981:2019; 

[min_prec, min_year] = min(annual_prec);
min_prec
years(min_year)

figure(11) 
stem(years, annual_prec,'Color',[0,0.7,0.9]);
title('Total Annual Precipitations');
xlabel('years');
ylabel('Total annual precipitations [mm]');
print('figure11','-vector','-dpdf') %Saved to 'figure11.pdf'...

% -------------------------------------------------------------------------
% # 5: Compute rainfall maxima of a certain duration
% -------------------------------------------------------------------------
%%

max = 0;
for k = 1:length(h)-2 %Iteration over the measurements vector
    s = sum(h(k:k+2)); %Sum of the elements in the moving window 
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
    indices = y == k; %Logical indexing with the vector "indices"
    yearly_prec = h(indices); %Selection of the measurment of the...
                              %indexed year
    max = 0;
    for l = 1:sum(indices)-2 %Setting the moving window range to the...
                             %number of measurements in the indexed year 
        s = sum(yearly_prec(l:l+2)); %Summing over the moving window
        if s >= max
            max = s;
        end
    end 
    three_hours_max(count,1) = max;
    count = count + 1;
end 

three_hours_max;


% -------------------------------------------------------------------------
% # 7 bis: Extend to multiple years and multiple durations
% -------------------------------------------------------------------------
%%

AnnualMax = zeros(39,6);

count = 1;
for k = 1981:2019
    indices = y == k; %Logical indexing with the vector "indices"
    yearly_prec = h(indices); %Selection of the measurment of the...
                              %indexed year
    D = [1 3 6 12 24 48];
    
    ind = 1;
    for m = 1:49
        if ismember(m, D) %To select only the desired durations 
            max = 0;
            for l = 1:sum(indices)+1-m %Defining a moving window the...
                                       %length of the duration
                s = sum(yearly_prec(l:l+m-1)); %Summing over the moving...
                                               %window 
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





% -------------------------------------------------------------------------
% # 8: save the output
% useful functions: save
% -------------------------------------------------------------------------
% when you are confident about your results, save the variables AnnualMax
% and D by uncommenting the following line: 
save('assignment1_output_part1.mat','D','AnnualMax');
