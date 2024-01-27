%% Initialize
clear;
clc;
close all;
%% Parameters
n = 1;
t = 0.001*2;
E = ;
r = 0.001*50;
S1 = 0.189394 + 0:899845*n - 0.4333*n^2 + 0:097866*n^3 - 0.00839*n^4;
S2 = -0.09799 + 0.982995*n - 0.96184*n^2 + 0.413319*n^3 - 0.08387*n^4 + 0.006530*n^5
%% Stress and Angle
theta = 1;
Sigma = E*t/(2*r)*(S1*theta+S2*theta^2);


