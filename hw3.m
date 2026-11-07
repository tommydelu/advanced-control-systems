%% EQUATIONS OF MOTION (dynamic model)

clearvars -except T_kin U B q1 q2 q3; clc; close all;

syms t Q1(t) Q2(t) Q3(t) tau1(t) tau2(t) tau3(t)
syms dq1 dq2 dq3 ddq1 ddq2 ddq3 real

q = [q1; q2; q3];
Q = [Q1; Q2; Q3];

dq  = [dq1;  dq2;  dq3];
ddq = [ddq1; ddq2; ddq3];
tau = [tau1; tau2; tau3];

% So I can derive B and U w.r.t time (use B and U to derive w.r.t q)
B_Q = subs(B, [q1 q2 q3], [Q1 Q2 Q3]); 
U_Q = subs(U, [q1 q2 q3], [Q1 Q2 Q3]);

fprintf("--------------------------------------\n")
fprintf("B SYMBOLIC FORMULATION:\n\n");
disp(vpa(B_Q,4))
fprintf("--------------------------------------\n\n")

lagrangian = T_kin - U;
fprintf("LAGRANGIAN SYMBOLIC FORMULATION:\n\n");
disp(vpa(lagrangian,4))
fprintf("--------------------------------------\n\n")

% Derivative of lagrangian w.r.t dq
term1 = B_Q*dq;
term1 = subs(term1, [q1 q2 q3], [Q1 Q2 Q3]);
fprintf("TERM 1:\n\n");
disp(vpa(term1,4))
fprintf("--------------------------------------\n\n")

% To derive term1 w.r.t time first I need to find B_d which is the derivative of B w.r.t time
B_d = diff(B_Q,t);
fprintf("Time derivative of B:\n\n");
disp(vpa(B_d,4))
fprintf("--------------------------------------\n\n")

% Time derivative of term1
term2 = B_Q*ddq + B_d*dq;
fprintf("TERM 2:\n\n");
disp(vpa(term2,4))
fprintf("--------------------------------------\n\n")

% Derivative of lagrangian w.r.t q = [q1 q2 q3]
term3 = 0.5.* [diff(dq'*B*dq, q1).'; diff(dq'*B*dq, q2).'; diff(dq'*B*dq, q3).'] - [ diff(U,q1).'; diff(U,q2).'; diff(U,q3).'];
term3 = subs(term3, [q1 q2 q3], [Q1 Q2 Q3]);
fprintf("TERM 3:\n\n");
disp(vpa(term3,4))
fprintf("--------------------------------------\n\n")

%%
