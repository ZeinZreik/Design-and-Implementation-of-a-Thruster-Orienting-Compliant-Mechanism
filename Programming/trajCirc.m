function [desired] = trajCirc(t, sec, radius, distance)
x = radius*cos(2*pi*t/sec);
y = radius*sin(2*pi*t/sec);
alpha = atan2(x,distance);
gamma = atan2(y,distance);
desired = [alpha gamma];
end

