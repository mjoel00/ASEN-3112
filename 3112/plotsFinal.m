%% File Header
% 

%% Housekeeping
clear
clc
close all
%% Read in Files
array1 = table2array(readtable('Torsion400inlbf_142.txt'));
array2 = table2array(readtable('Torsion20inlbf_144.txt'));

%% Analysis of the Closed Thin Wall Specimen
R_e = 3/8; % Exterior radius [in]
L_ext = 10; % Extensometer gauge length [in]

torque1 = array1(:,4);
shearStrainExten1 = (array1(:,3) * (pi/180));
phi1 = (array1(:,2) - array2(1,2)) * (pi/180);

[leastSquaresExten1,SExten1] = polyfit(shearStrainExten1,torque1,1);
fitExten1 = polyval(leastSquaresExten1, shearStrainExten1);
GJExten1Err = sqrt(diag((SExten1.R)\inv(SExten1.R'))./SExten1.normr.^2/SExten1.df);
GJExten1Err = GJExten1Err(1);
GJ_CTW_Exten = leastSquaresExten1(1);

% Plot torque vs. shearStrain by Extensometer
figure(1)
plot(shearStrainExten1,torque1);
hold on;
plot(shearStrainExten1,fitExten1,'LineWidth', 2);
title('Torque vs. Shear Strain Extensometer (CTW)');
xlabel('Shear Strain (rad/in)');
ylabel('Torque (lbf*in)');
legend('Extensometer', 'Polyfit', 'Location', 'southeast');

shearStrainCalc1 = (phi1 * R_e) / L_ext;

[leastSquaresCalc1, SCalc1] = polyfit(shearStrainCalc1,torque1,1);
fitCalc1 = polyval(leastSquaresCalc1,shearStrainCalc1);
GJCalc1Err = sqrt(diag((SCalc1.R)\inv(SCalc1.R'))./SCalc1.normr.^2/SCalc1.df);
GJCalc1Err = GJCalc1Err(1);
GJ_CTW_Calc = leastSquaresCalc1(1);

% Plot torque vs. shearStrain Calculated (using total rotation angle)
figure(2)
plot(shearStrainCalc1,torque1);
hold on;
plot(shearStrainCalc1,fitCalc1, 'LineWidth', 2);
title('Torque vs. Shear Strain Calculated (CTW)');
xlabel('Shear Strain (rad/in)');
ylabel('Torque (lbf*in)');
legend('Calculated (Machine)', 'Polyfit', 'Location', 'southeast');


%% Analysis of the Open Thin Wall Specimen
thickness = 1/16; % Wall thickness [in]

torque2 = array2(:,4);
shearStrainExten2 = (array2(:,3) * (pi/180));
phi2 = (array2(:,2) - array2(1,2)) * (pi/180);

[leastSquaresExten2,SExten2] = polyfit(shearStrainExten2,torque2,1);
fitExten2 = polyval(leastSquaresExten2,shearStrainExten2);
GJExten2Err = sqrt(diag((SExten2.R)\inv(SExten2.R'))./SExten2.normr.^2/SExten2.df);
GJ_OTW_Exten_Err = GJExten2Err(1);
GJ_OTW_Exten = leastSquaresExten2(1);

figure(3)
plot(shearStrainExten2,torque2);
hold on
plot(shearStrainExten2,fitExten2, 'LineWidth', 2)
title('Torque vs. Shear Strain Extensometer (OTW)');
xlabel('Shear Strain (rad/in)');
ylabel('Torque (lbf*in)');
legend('Extensometer', 'Polyfit', 'Location', 'southeast');

shearStrainCalc2 = (phi2 * thickness) / L_ext;

[leastSquaresCalc2, SCalc2] = polyfit(shearStrainCalc2,torque2,1);
fitCalc2 = polyval(leastSquaresCalc2,shearStrainCalc2);
GJCalc2Err = sqrt(diag((SCalc2.R)\inv(SCalc2.R'))./SCalc2.normr.^2/SCalc2.df);
GJ_OTW_Calc_Err = GJCalc2Err(1);
GJ_OTW_Calc = leastSquaresCalc2(1);

figure(4)
plot(shearStrainCalc2,torque2);
hold on
plot(shearStrainCalc2,fitCalc2, 'LineWidth', 2)
title('Torque vs. Shear Strain Calculated (OTW)');
xlabel('Shear Strain (rad/in)');
ylabel('Torque (lbf*in)');
legend('Calculated (Machine)', 'Polyfit', 'Location', 'southeast');