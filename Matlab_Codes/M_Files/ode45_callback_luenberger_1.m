function output = ode45_callback_luenberger_1(t,x)
M = 1000;
m1 = 100;
m2 = 100;
theta1 = x(3);
theta2 = x(5);
theta_dot_1 = x(4);
theta_dot_2 = x(6);
l1 = 20;
l2 = 10;
g = 9.8;
[A,B,Q,R] = AB();
[C_1,~,~] = get_C_matrices();
[K,~,~] = lqr(A,B,Q,R);
F= -K*x(1:6);

output = zeros(12,1);

x_ddot = (F-(m1*sind(theta1))*(g*cosd(theta1)+l1*theta_dot_1*theta_dot_1)-m2*sind(theta2)*(g*cosd(theta2)+l2*theta_dot_2*theta_dot_2))/(M+m1*sind(theta1)*sind(theta1)+m2*sind(theta2)*sind(theta2));
theta1_ddot = (x_ddot*cosd(theta1)-g*sind(theta1))/l1;
theta2_ddot = (x_ddot*cosd(theta2)-g*sind(theta2))/l2;


req_poles1 = [-10;-20;-30;-40;-50;-60];
Luenberger1 = place(A',C_1',req_poles1);
L1 = Luenberger1';

est = (A-L1*C_1)*x(7:12);

output(1) = x(2); %because initial state has x,x_dot,theta1,theta_dot_1,theta2,theta_dot_2,Estimation states (6)
output(2) = x_ddot;
output(3) = x(4);
output(4) = theta1_ddot;
output(5) = x(6);
output(6) = theta2_ddot;
output(7) = est(1);
output(8) = est(2);
output(9) = est(3);
output(10) = est(4);
output(11) = est(5);
output(12) = est(6);
end