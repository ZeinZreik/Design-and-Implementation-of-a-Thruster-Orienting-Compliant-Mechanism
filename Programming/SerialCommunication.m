clear;
clc;
close all;
instrreset;
%% Start
mech = serial('COM6','BaudRate',115200);
fopen(mech);
%% 
t = 0;
sec = 0.1;
distance = 0.1;
side = 0.1;
tic;
while t<=1*sec
    desired = traj(t*sec, sec, side, distance);
    theta = inv_kinematics(desired);
    out = string(theta(1))+','+string(theta(2));
    fprintf(mech,out);
    t = toc;
end
tic;
while 1==1
    desired = traj(2*sec+t*sec, sec, side, distance);
    theta = inv_kinematics(desired);
    out = string(theta(1))+','+string(theta(2));
    fprintf(mech,out);
    t = toc;
    if t>=4
        t = 0;
        tic;
    end
end
%% End
fclose(mech);
%% 
t = 0;
tic;
while t<=10
    desired = traj(t);
    theta2 = sin(pi*t)*2*pi;
    theta1 = sin(pi*t)*2*pi;
    out = string(theta1)+','+string(theta2);
    fprintf(mech,out);
    t = toc;
end
%% 
fprintf(mech,"on");
distance = 0.16;
x = -0.04;
y = 0;
alpha = atan2(x,distance);
gamma = atan2(y,distance);
desired = [alpha gamma];
theta = inv_kinematics(desired);
fprintf(mech, string(theta(1))+','+string(theta(2)));