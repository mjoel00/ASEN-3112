% ASEN 3112 - Lab 2
%
% Linear Regression & Uncertainty analysis
% Matthew Pabin, Jack Davis, Alicia Wu, Nathan Braunstein

clear
clc
close all

% Read experimental data
data = readtable('Lab');
data = table2array(data);

loadCase = data(:,1) * 4.448;  % N 
F0 = data(:,2) * 4.448;
F1 = data(:,3) * 4.448;
F2 = data(:,4) * 4.448;      
F3 = data(:,5)* 4.448;     
LVDT = data(:,6) * 25.4; %mm 

% read ANSYS data
Rdata = readtable('Lab_Reac');
Rdata = table2array(Rdata);

%% Load vs Displacement

% Linear Regression
[linFit_disp,S_disp] = polyfit(loadCase,LVDT,1);
[y_disp,delta_disp] = polyval(linFit_disp,loadCase,S_disp);

% Plot
figure
plot(loadCase,LVDT,'b.','MarkerSize',15)
xlabel('Load [N]')
ylabel('Displacement [mm]')
title('Load vs Displacement of Truss')
hold on 
plot(loadCase,y_disp,'g-')
errorbar(loadCase,y_disp,delta_disp,'r.')
legend('Exp. Value','Linear Fit','Error')

%% Load vs Forces

[linFit_F0,S_F0] = polyfit(loadCase,F0,1);
[y_F0,delta_F0] = polyval(linFit_F0,loadCase,S_F0);

[linFit_F1,S_F1] = polyfit(loadCase,F1,1);
[y_F1,delta_F1] = polyval(linFit_F1,loadCase,S_F1);

[linFit_F2,S_F2] = polyfit(loadCase,F2,1);
[y_F2,delta_F2] = polyval(linFit_F2,loadCase,S_F2);

[linFit_F3,S_F3] = polyfit(loadCase,F3,1);
[y_F3,delta_F3] = polyval(linFit_F3,loadCase,S_F3);

figure
plot(loadCase,F0,'ro','MarkerSize',8)
hold on
plot(loadCase,F1,'b.','MarkerSize',15)
hold on 
plot(loadCase,F2,'magx','MarkerSize',15)
xlabel('Load [N]')
ylabel('Reaction Forces [N]')
title('Load vs Reaction Forces of Truss')
hold on 
plot(loadCase,y_F0,'g-')
%errorbar(loadCase,y_F0,delta_F0,'r.')
hold on 
plot(loadCase,y_F1,'g-')
%errorbar(loadCase,y_F1,delta_F1,'r.')
hold on 
plot(loadCase,y_F2,'g-')
%errorbar(loadCase,y_F2,delta_F1,'r.')
legend('F0','F1','F2','Linear Fit')
hold off 

figure
plot(loadCase,F3,'mag.','MarkerSize',15)
xlabel('Load [N]')
ylabel('Internal Force [N]')
title('Load vs Internal Force of Truss')
hold on 
plot(loadCase,y_F3,'g-')
errorbar(loadCase,y_F3,delta_F3,'r.')
legend('F3D','Linear Fit','Error')

%% ANSYS
RloadCase = 4.448 * [0 10 20 30 40 50 40 30 20 10 0];

% Displacement
D = 0.0015439;   %m
RLVDT(1) = 0;
RLVDT(2) = (1/5) * D;
RLVDT(3) = (2/5) * D;
RLVDT(4) = (3/5) * D;
RLVDT(5) = (4/5) * D;
RLVDT(6) = (5/5) * D;
RLVDT(7) = (4/5) * D;
RLVDT(8) = (3/5) * D;
RLVDT(9) = (2/5) * D;
RLVDT(10) = (1/5) * D;
RLVDT(11) = 0;  
RLVDT = RLVDT * (-1000);

[linFit_Rdisp,S_Rdisp] = polyfit(RloadCase,RLVDT,1);
[y_Rdisp,delta_Rdisp] = polyval(linFit_Rdisp,RloadCase,S_Rdisp);

figure
plot(RloadCase,RLVDT,'b.','MarkerSize',15)
xlabel('Load [N]')
ylabel('Displacement [mm]')
title('Load vs Displacement of Truss (ANSYS)')
hold on 
plot(RloadCase,y_Rdisp,'g-')
%errorbar(RloadCase,y_Rdisp,delta_Rdisp,'r.')
legend('Exp. Value','Linear Fit','Error')

% Reac Forces
f0 = sqrt( Rdata(1,2)^2 + Rdata(1,3)^2 + Rdata(1,4)^2);
f1 = sqrt( Rdata(2,2)^2 + Rdata(2,3)^2 + Rdata(2,4)^2);
f2 = sqrt( (Rdata(3,2)+Rdata(4,2))^2 + (Rdata(3,3)+Rdata(4,3))^2 + (Rdata(3,4)+Rdata(4,4))^2);

RF0(1) = 0;
RF0(2) = (1/5) * f0;
RF0(3) = (2/5) * f0;
RF0(4) = (3/5) * f0;
RF0(5) = (4/5) * f0;
RF0(6) = (5/5) * f0;
RF0(7) = (4/5) * f0;
RF0(8) = (3/5) * f0;
RF0(9) = (2/5) * f0;
RF0(10)= (1/5) * f0;
RF0(11) = 0;

RF1(1) = 0;
RF1(2) = (1/5) * f1;
RF1(3) = (2/5) * f1;
RF1(4) = (3/5) * f1;
RF1(5) = (4/5) * f1;
RF1(6) = (5/5) * f1;
RF1(7) = (4/5) * f1;
RF1(8) = (3/5) * f1;
RF1(9) = (2/5) * f1;
RF1(10)= (1/5) * f1;
RF1(11) = 0;


RF2(1) = 0;
RF2(2) = (1/5) * f2;
RF2(3) = (2/5) * f2;
RF2(4) = (3/5) * f2;
RF2(5) = (4/5) * f2;
RF2(6) = (5/5) * f2;
RF2(7) = (4/5) * f2;
RF2(8) = (3/5) * f2;
RF2(9) = (2/5) * f2;
RF2(10)= (1/5) * f2;
RF2(11) = 0;

[linFit_RF0,S_RF0] = polyfit(RloadCase,RF0,1);
[y_RF0,delta_RF0] = polyval(linFit_RF0,RloadCase,S_RF0);

[linFit_RF1,S_RF1] = polyfit(RloadCase,RF1,1);
[y_RF1,delta_RF1] = polyval(linFit_RF1,RloadCase,S_RF1);

[linFit_RF2,S_RF2] = polyfit(RloadCase,RF2,1);
[y_RF2,delta_RF2] = polyval(linFit_RF2,RloadCase,S_RF2);

figure
plot(RloadCase,RF0,'ro','MarkerSize',8)
hold on
plot(RloadCase,RF1,'b.','MarkerSize',15)
hold on 
plot(RloadCase,RF2,'magx','MarkerSize',15)
xlabel('Load [N]')
ylabel('Reaction Forces [N]')
title('Load vs Reaction Forces of Truss (ANSYS)')
hold on 
plot(RloadCase,y_RF0,'g-')
%errorbar(RloadCase,y_RF0,delta_RF0,'r.')
hold on 
plot(RloadCase,y_RF1,'g-')
%errorbar(RloadCase,y_RF1,delta_RF1,'r.')
hold on 
plot(RloadCase,y_RF2,'g-')
%errorbar(RloadCase,y_RF2,delta_RF1,'r.')
legend('F0','F1','F2','Linear Fit')


% Int Forces
f3 = 450.14; % N ->  from ANSYS simulation

RF3(1) = 0;
RF3(2) = (1/5) * f3;
RF3(3) = (2/5) * f3;
RF3(4) = (3/5) * f3;
RF3(5) = (4/5) * f3;
RF3(6) = (5/5) * f3;
RF3(7) = (4/5) * f3;
RF3(8) = (3/5) * f3;
RF3(9) = (2/5) * f3;
RF3(10)= (1/5) * f3;
RF3(11) = 0;


[linFit_RF3,S_RF3] = polyfit(RloadCase,RF3,1);
[y_RF3,delta_RF3] = polyval(linFit_RF3,RloadCase,S_RF3);

figure
plot(RloadCase,RF3,'mag.','MarkerSize',15)
xlabel('Load [N]')
ylabel('Internal Force [N]')
title('Load vs Internal Force of Truss (ANSYS)')
hold on 
plot(RloadCase,y_RF3,'g-')
%errorbar(RloadCase,y_RF3,delta_RF3,'r.')
legend('F3D','Linear Fit','Error')

