function [desired] = traj(t, sec, side, distance)
xlength = side;
ylength = side;
x0 = 0;
y0 = 0;
xd = 0;
yd = 0;
if t<=1*sec
    t = t*sec;
    x0 = -xlength/2;
    y0 = -ylength/2;
    xd = +xlength/2;
    yd = -ylength/2;
elseif t<=2*sec
    t = t-1*sec;
    x0 = +xlength/2;
    y0 = -ylength/2;
    xd = +xlength/2;
    yd = +ylength/2;
elseif t<=3*sec
    t = t-2*sec;
    x0 = +xlength/2;
    y0 = +ylength/2;
    xd = -xlength/2;
    yd = +ylength/2;
else 
    t = t-3*sec;
    x0 = -xlength/2;
    y0 = +ylength/2;
    xd = -xlength/2;
    yd = -ylength/2;
end
% y0=0;yd=0;
% vx = sqrt((xd-x0)^2+(yd-y0)^2)/sec;
% vy = vx;
vx = (xd-x0)/sec;
vy = (yd-y0)/sec;
x = x0 + vx*t;
y = y0 + vy*t;
alpha = atan2(x,distance);
gamma = atan2(y,distance);
desired = [alpha gamma];
end

