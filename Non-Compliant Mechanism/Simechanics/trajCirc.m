function [desired] = trajCirc(t)
distance = 0.112;
radius = 0.01;
sec = 6;
x = radius*cos(2*pi*t/sec);
y = radius*sin(2*pi*t/sec);
alpha = atan2(x,distance);
gamma = atan2(y,distance);
desired = [alpha gamma];
end

