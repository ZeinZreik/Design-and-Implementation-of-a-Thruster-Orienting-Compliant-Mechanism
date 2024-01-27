function angles = angle_convert(q)
q = reshape(q, 1, 4);
[roll, pitch, yaw] = quat2angle(q, "XYZ");
angles = [roll pitch yaw]*180/pi;
angles = round(angles,3);
end

