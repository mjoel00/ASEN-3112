%ASEN 3112 Lab 3 
%Problem IV.2
%Maklen Estrada
%Mathew Davis
%Tyler Soper
%Samuel Hatton
%Matthew Pabin
%Elena B
%Hrithik Hiranandani
clear ;clc
%% Variables and Constants
L = 12; %in
L_E = 4.5; %in
L_R = 5;
w = 1;
h = 1/8;
h_E = 1/4;
h_R = 0.040;
E = 10175000; %psi
rho = 0.0002505; %lb-sec^2/in^4
M_T = 1.131*rho;
S_T = 0.5655*rho;
I_T = 23.124*rho;
A = w*h;
I_z= (w*h^3)/12;
C_M4 = (rho*A*L)/806400;
C_k4 = (8*E*I_z)/L^3;

%% Master Stiff
MasterStiff = [77088 2916*L 23712 -1284*L 0 0 0 0 0 0;
2916*L 172*L^2 1284*L -73*L^2 0 0 0 0 0 0;
23712 1284*L 154176 0 23712 -1284*L 0 0 0 0;
-1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0 0 0;
0 0 23712 1284*L 154176 0 23712 -1284*L 0 0;
0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0;
0 0 0 0 23712 1284*L 154176 0 23712 -1284*L;
0 0 0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2;
0 0 0 0 0 0 23712 1284*L 77088 -2916*L;
0 0 0 0 0 0 -1284*L -73*L^2 -2916*L 172*L^2];

%% Add
Add= [0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 M_T S_T;
0 0 0 0 0 0 0 0 S_T I_T];

%% K4
K4 = [96 12*L -96 12*L 0 0 0 0 0 0;
12*L 2*L^2 -12*L L^2 0 0 0 0 0 0;
-96 -12*L 192 0 -96 12*L 0 0 0 0;
12*L L^2 0 4*L^2 -12*L L^2 0 0 0 0;
0 0 -96 -12*L 192 0 -96 12*L 0 0;
0 0 12*L L^2 0 4*L^2 -12*L L^2 0 0;
0 0 0 0 -96 -12*L 192 0 -96 12*L;
0 0 0 0 12*L L^2 0 4*L^2 -12*L L^2;
0 0 0 0 0 0 -96 -12*L 96 -12*L;
0 0 0 0 0 0 12*L L^2 -12*L 2*L^2];

%% Eigenvalue Prob
syms  w
f(w) = K4 - w^2*(MasterStiff + Add);
s = det(f);
syms x
x = sqrt(w^2);
s = s(x) ;
eqn = s == 0;
sol = solve(eqn);
%% W's
w1 = sol(1);
w2 = sol(2);
w3 = sol(3);
w4 = sol(4);
w5 = sol(5);
w6 = sol(6);
w7 = sol(7);
w8 = sol(8);
w9 = sol(9);
w10 = sol(10);
warray = [w1 w2 w3 w4 w5 w6 w7 w8 w9 w10]';
Wsol = double(warray);

%% Frequencies
Freq = Wsol/(2*pi); %Hz
%Three Smallest Freq
F1 = Freq(3);
F2 = Freq(4);
F3 = Freq(5);

%% Q3

data = readtable('test_2min_all_4');
data = table2array(data);
l = 1:1000:12;
d1 = data(:,7)/100;
d2 = data(:,8)/100;
d3 = data(:,9)/100;

[U,D] = eig(K4,(MasterStiff + Add));

new_U = [[0;0;U(:,8)]/max(abs(U(:,8))),[0;0;U(:,7)]/max(abs(U(:,7))),...
        [0;0;U(:,6)]/max(abs(U(:,6)))];
    
ploteigenvector(L,new_U(:,1),4,10,1)
hold on
% plot(l,d1,'r')
% title('Mode 1: 12.02 Hz')
% xlabel('Length [mm]')
% ylabel('Amplitude')

% ploteigenvector(L,new_U(:,2),4,10,1)
% hold on
% plot(l,d2)
% title('Mode 2: 15.05 Hz')
% xlabel('Length [mm]')
% ylabel('Amplitude')
% 
% ploteigenvector(L,new_U(:,3),4,10,1)
% hold on
% plot(l,d3)
% title('Mode 3: 203.45 Hz')
% xlabel('Length [mm]')
% ylabel('Amplitude')
