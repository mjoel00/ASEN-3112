%% ASEN 3112 Lab 2
% Data Analysis for Loads on a 16-bay Truss
% Matthew Pabin; Jack Davis; Kyler Stirewalt; Alicia Wu; Nathan Braunstein

% Housekeeping
clear; clc; close all;

%==========================================================================
%% Read Data 
data = readtable('Lab');
data = table2array(data);

% Format Data to Variables
loadCase = data(:,1);
F0 = data(:,2);
F1 = data(:,3);
F2 = data(:,4);
F3D = data(:,5);
LVDT = data(:,6);

%==========================================================================
%% Load -vs.- Displacement
% Linear Regression
[linFit,S] = polyfit(loadCase,LVDT,1);
[y,delta] = polyval(linFit,loadCase,S);

% Plotting
plot(loadCase,LVDT,'b.','MarkerSize',15)
xlabel('Load [lbf]')
ylabel('Displacement [in]')
title('Load vs Displacement of Truss')
hold on 
plot(loadCase,y,'g-','LineWidth',1.5)
errorbar(loadCase,y,delta,'r.')
legend('Exp. Value','Linear Fit','Error')
axis padded
