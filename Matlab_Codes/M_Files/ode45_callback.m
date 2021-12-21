function output = ode45_callback(t,x)
[A,B,Q,R] = AB();
[K,~,~] = lqr(A,B,Q,R);
%We are giving feedback U=-KX where X is the initial state Now
F = -K*x;
output = zeros(6,1);
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

x_ddot = (F-(m1*sind(theta1))*(g*cosd(theta1)+l1*theta_dot_1*theta_dot_1)-m2*sind(theta2)*(g*cosd(theta2)+l2*theta_dot_2*theta_dot_2))/(M+m1*sind(theta1)*sind(theta1)+m2*sind(theta2)*sind(theta2));
theta1_ddot = (x_ddot*cosd(theta1)-g*sind(theta1))/l1;
theta2_ddot = (x_ddot*cosd(theta2)-g*sind(theta2))/l2;
%look at the non linear system representation in Chapter X
output(1) = x(2); %because initial state has x,x_dot,theta1,theta_dot_1,theta2,theta_dot_2
output(2) = x_ddot;
output(3) = x(4);
output(4) = theta1_ddot;
output(5) = x(6);
output(6) = theta2_ddot;
end