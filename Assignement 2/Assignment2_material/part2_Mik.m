
clear variables %clear the workspace variables
close all %close all figures
clc %clear the command window

% importing data

load('output_part1.mat');
%% (1)

n_sub = 5; %each hour is broken into n_int intervals
dt = 1/n_sub; % timestep [h]

t_Je= linspace(0,4-1,4);
t_Jedt = linspace(0,4-dt,4/dt);
Jedt = interp1(t_Je, Je, t_Jedt ,'previous','extrap');

%% (2)
% loading IUH params
load('../Parameters/IUHpars.mat')
% generating watershed IUH gamma distribution
%gammaIUH = zeros(size(Jedt));

gammaIUH = pdf('Gamma',40, par_shape, par_scale);




