%ASEN 3112 Lab 3 
%Problem IV.2
%Maklen Estrada
%Mathew Davis
%Tyler Soper
%Samuel Hatton
%Matthew Pabin
%Elena B
%Hrithik Hiranandani
clear ;clc; close all
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
MasterStiff = C_M4*[77088 2916*L 23712 -1284*L 0 0 0 0 0 0;
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
K4 = C_k4*[96 12*L -96 12*L 0 0 0 0 0 0;
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
% syms x
% x = sqrt(w);
% s = s(x) ;
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

K4_new = K4(3:end,3:end);
M4_new = MasterStiff(3:end, 3:end);
[U,el] = eig(K4_new,M4_new);
omega =real(sqrt(diag(el)));
omega = omega/(2*pi);
nat_freq = omega(1:3);
%[frq,num] = sort(1/2/pi*omega);

U_temp = U(:,num);

U = zeros(size(K4));

U(3:end,3:end) = U_temp;
ev = U;

nsub=10;
scale=1.0;
ne = 4;
nv=ne*nsub+1;
Le=L/ne;
dx=Le/nsub;
figure(1);

for t=3:5 
    k=0;
    x=zeros(nv,1);
    v=zeros(nv,1);
    
    for e = 1:ne
        xxi=Le*(e-1);
        vi=ev(2*e-1,t);
        teti=ev(2*e,t);
        vj=ev(2*e+1,t);
        tetj=ev(2*e+2,t);
        
        if (e==1); ni=0;else; ni=1;end
        for n=ni:nsub
            xk=xxi+dx*n;
            xi=(2*n-nsub)/nsub;
            vk=scale*(0.125*(4*(vi+vj)+2*(vi-vj)*(xi^2-3)*xi+Le*(xi^2-1)*(tetj-teti+(teti+tetj)*xi)));
            k=k+1;
            x(k)=xk;
            v(k)=vk;
        end
    end
    plot(x,v);hold on; grid on;
    xlabel('Length [in]')
    ylabel('Normalized Displacement in %')
    title('Modes for 4 Element FEM')
    legend('Mode1', 'Mode 2', 'Mode 3')
end

data = importdata('test_2min_all_4');
data = double(data.data);
