clear;
clc;
close all;
%% Parameters
E = 2.41e9;
sigma_max = 72.3e6;
t = 2e-3;
r = 50e-3;
w = 50e-3;
l = 67.88e-3;
b = 5e-3;
n = r/w;
I = b*t^3/12;
S1 = 0.189394+0.899845*n-0.4333*n^2+0.097866*n^3-0.00839*n^4;
S2 = -0.09799+0.982995*n-0.96184*n^2+0.413319*n^3-0.08387*n^4+0.00653*n^5;
%% Find Theta Max
syms theta_max
eq = (E*t/2/r*(S1*theta_max+S2*theta_max^2) == sigma_max);
Sol_theta = double(vpasolve(eq,theta_max));
Sol_theta = atan(tan(Sol_theta));
disp('Theta Max = '+string(Sol_theta*180/pi));
%% Moment for 15 degrees
theta = 15*pi/180;
K_Theta = 5.300185-1.6866*n+0.885356*n^2-0.2094*n^3+0.018385*n^4;
K = K_Theta*E*I/2/l;
M = theta*K;
disp('M = ' + string(M));