%3112 Lab 1
% Matthew Pabin

clc 
clear all

%% CTW
c_data = readtable('Torsion400inlbf_142.txt');
c_data = table2array(c_data);

c_phi = c_data(:,2) - 48.0018;
c_gamma = c_data(:,3);
c_epsilon = c_data(:,5);
c_torq = c_data(:,4);

figure
scatter(c_torq,c_gamma)
title('Torque vs Shear Strain of CTW section')
xlabel('Torque (in-lbf)')
ylabel('Shear Strain')

c_leastsq = polyfit(c_torq,c_gamma,1);
c_y = c_leastsq(1)*c_torq + c_leastsq(2);
c_GJerr = sqrt(diag((SExten1.R)\inv(SExten1.R'))./SExten1.normr.^2/SExten1.df);
hold on
plot(c_torq,c_y,'linewidth',3)
hold off





%% OTW

o_data = readtable('Torsion20inlbf_144.txt');
o_data = table2array(o_data);

o_phi = o_data(:,2);
o_gamma = o_data(:,3);
o_epsilon = o_data(:,5);
o_torq = o_data(:,4);

figure
scatter(o_torq,o_gamma)
title('Torque vs Shear Strain of OTW section')
xlabel('Torque (in-lbf)')
ylabel('Shear Strain')

o_leastsq = polyfit(o_torq,o_gamma,1);
o_y = o_leastsq(1)*o_torq + o_leastsq(2);
hold on
plot(o_torq,o_y,'linewidth',3)



