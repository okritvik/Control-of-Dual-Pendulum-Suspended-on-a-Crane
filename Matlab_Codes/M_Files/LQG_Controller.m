% LQG For Linear Model
clc;
clear all;
close all;
[A,B,Q,R] = AB();
[C_1,C_3,C_4] = get_C_matrices();
D=0;
QXU = eye(7);
% Generate noise
w = wgn(6,1,5);
v = wgn(1,1,5);
QWV = [w;v]*[w' v'];
%Since we were asked to take only the smallest output vector, we use the
% C_1 which is having the lowest order in the observable states.
%Note that as the noise is white gaussian random noise, the output may
%change everytime we run the code.
LQG_state_space = lqg(ss(A,B,C_1,D),QXU,QWV)
initial_state = [3,0.3,20,1,10,2];
figure
initial(LQG_state_space,initial_state)
title('LQG Controller for Linear Model')
xlabel('Time')
ylabel('State Output x')
%Uncomment Below lines to get LQG controller for all states
% w = wgn(11,1,5);
% v = wgn(1,1,5);
% QWV = [w;v]*[w' v'];
% LQG_state_space = lqg(ss(A,B,eye(6),D),QXU,QWV)
% initial_state = [3,0.3,20,1,10,2];
% figure
% initial(LQG_state_space,initial_state)
% title('LQG Controller for Linear Model')
% xlabel('Time')
% ylabel('State Outputs')

%LQG (using kalman filter gain) for non-linear system model
initial_state = [3,0.3,20,1,10,2,0,0,0,0,0,0];
simulation_time = 0:1:2000;
[time,out] = ode45(@ode45_callback_lqg,simulation_time,initial_state);
figure
plot(time,out)
title('LQG Controller for Non Linear Model')
xlabel('Time')
ylabel('State Estimates')