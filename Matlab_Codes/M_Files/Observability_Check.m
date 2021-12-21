% Observability 
% Jacobian Linearization of the Non Linear System of Dual Pendulum
% Suspended on a crane.
clc;
clear all;
close all;
syms x x_dot theta1 theta_dot_1 theta2 theta_dot_2 F M m1 m2 g l1 l2


x_ddot = (F-(m1*sin(theta1))*(g*cos(theta1)+l1*theta_dot_1*theta_dot_1)-m2*sin(theta2)*(g*cos(theta2)+l2*theta_dot_2*theta_dot_2))/(M+m1*sin(theta1)*sin(theta1)+m2*sin(theta2)*sin(theta2));
theta1_ddot = (x_ddot*cos(theta1)-g*sin(theta1))/l1;
theta2_ddot = (x_ddot*cos(theta2)-g*sin(theta2))/l2;


A_intermediate = [diff(x_dot,x) diff(x_dot,x_dot) diff(x_dot,theta1) diff(x_dot,theta_dot_1) diff(x_dot,theta2) diff(x_dot,theta_dot_2);
    diff(x_ddot,x) diff(x_ddot,x_dot) diff(x_ddot,theta1) diff(x_ddot,theta_dot_1) diff(x_ddot,theta2) diff(x_ddot,theta_dot_2);
    diff(theta_dot_1,x) diff(theta_dot_1,x_dot) diff(theta_dot_1,theta1) diff(theta_dot_1,theta_dot_1) diff(theta_dot_1,theta2) diff(theta_dot_1,theta_dot_2);
    diff(theta1_ddot,x) diff(theta1_ddot,x_dot) diff(theta1_ddot,theta1) diff(theta1_ddot,theta_dot_1) diff(theta1_ddot,theta2) diff(theta1_ddot,theta_dot_2);
    diff(theta_dot_2,x) diff(theta_dot_2,x_dot) diff(theta_dot_2,theta1) diff(theta_dot_2,theta_dot_1) diff(theta_dot_2,theta2) diff(theta_dot_2,theta_dot_2);
    diff(theta2_ddot,x) diff(theta2_ddot,x_dot) diff(theta2_ddot,theta1) diff(theta2_ddot,theta_dot_1) diff(theta2_ddot,theta2) diff(theta2_ddot,theta_dot_2);
    ];
A = subs(A_intermediate,[x,x_dot,theta1,theta_dot_1,theta2,theta_dot_2],[0,0,0,0,0,0])
B = [x_dot;x_ddot;theta_dot_1;theta1_ddot;theta_dot_2;theta2_ddot];
B = subs(B,[x,x_dot,theta1,theta_dot_1,theta2,theta_dot_2,F],[0,0,0,0,0,0,1])

% From the given problem statement, the output vectors are
% x,(theta1,theta2),(x,theta2),(x,theta1,theta2)
C_1 = [1 0 0 0 0 0]; %because Y = CX + DU. We take C matrix such that it accounts for the state variables
C_2 = [0 0 1 0 0 0;
       0 0 0 0 1 0];
C_3 = [1 0 0 0 0 0;
       0 0 0 0 1 0];
C_4 = [1 0 0 0 0 0;
       0 0 1 0 0 0;
       0 0 0 0 1 0];
% We know that the system is observable when rank[C' A'C'..... A'^nC']=n

O1 = [C_1' A'*C_1' A'*A'*C_1' A'*A'*A'*C_1' A'*A'*A'*A'*C_1' A'*A'*A'*A'*A'*C_1'];
O2 = [C_2' A'*C_2' A'*A'*C_2' A'*A'*A'*C_2' A'*A'*A'*A'*C_2' A'*A'*A'*A'*A'*C_2'];
O3 = [C_3' A'*C_3' A'*A'*C_3' A'*A'*A'*C_3' A'*A'*A'*A'*C_3' A'*A'*A'*A'*A'*C_3'];
O4 = [C_4' A'*C_4' A'*A'*C_4' A'*A'*A'*C_4' A'*A'*A'*A'*C_4' A'*A'*A'*A'*A'*C_4'];
display(['The rank of observability matrix for the choice of state x = ',num2str(rank(O1))])
display(['The rank of observability matrix for the choice of state (theta1,theta2) = ',num2str(rank(O2))])
display(['The rank of observability matrix for the choice of state (x,theta2) = ',num2str(rank(O3))])
display(['The rank of observability matrix for the choice of state (x,theta1,theta2) = ',num2str(rank(O4))])