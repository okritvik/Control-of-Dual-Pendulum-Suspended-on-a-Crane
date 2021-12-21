% Jacobian Linearization of the Non Linear System of Dual Pendulum
% Suspended on a crane.
clc;
clear all;
close all;
syms x x_dot theta1 theta_dot_1 theta2 theta_dot_2 F M m1 m2 g l1 l2


x_ddot = (F-(m1*sin(theta1))*(g*cos(theta1)+l1*theta_dot_1*theta_dot_1)-m2*sin(theta2)*(g*cos(theta2)+l2*theta_dot_2*theta_dot_2))/(M+m1*sin(theta1)*sin(theta1)+m2*sin(theta2)*sin(theta2));
theta1_ddot = (x_ddot*cos(theta1)-g*sin(theta1))/l1
theta2_ddot = (x_ddot*cos(theta2)-g*sin(theta2))/l2;


A_intermediate = [diff(x_dot,x) diff(x_dot,x_dot) diff(x_dot,theta1) diff(x_dot,theta_dot_1) diff(x_dot,theta2) diff(x_dot,theta_dot_2);
    diff(x_ddot,x) diff(x_ddot,x_dot) diff(x_ddot,theta1) diff(x_ddot,theta_dot_1) diff(x_ddot,theta2) diff(x_ddot,theta_dot_2);
    diff(theta_dot_1,x) diff(theta_dot_1,x_dot) diff(theta_dot_1,theta1) diff(theta_dot_1,theta_dot_1) diff(theta_dot_1,theta2) diff(theta_dot_1,theta_dot_2);
    diff(theta1_ddot,x) diff(theta1_ddot,x_dot) diff(theta1_ddot,theta1) diff(theta1_ddot,theta_dot_1) diff(theta1_ddot,theta2) diff(theta1_ddot,theta_dot_2);
    diff(theta_dot_2,x) diff(theta_dot_2,x_dot) diff(theta_dot_2,theta1) diff(theta_dot_2,theta_dot_1) diff(theta_dot_2,theta2) diff(theta_dot_2,theta_dot_2);
    diff(theta2_ddot,x) diff(theta2_ddot,x_dot) diff(theta2_ddot,theta1) diff(theta2_ddot,theta_dot_1) diff(theta2_ddot,theta2) diff(theta2_ddot,theta_dot_2);
    ]
A = subs(A_intermediate,[x,x_dot,theta1,theta_dot_1,theta2,theta_dot_2],[0,0,0,0,0,0])
B = [x_dot;x_ddot;theta_dot_1;theta1_ddot;theta_dot_2;theta2_ddot];
B = subs(B,[x,x_dot,theta1,theta_dot_1,theta2,theta_dot_2,F],[0,0,0,0,0,0,1])