function [theta] = inv_kinematics(desired)
alpha = desired(1);
gamma = desired(2);
theta1 = atan(-tan(gamma));
theta2 = atan(-tan(alpha)*cos(gamma));
theta = [theta1 theta2];
end

