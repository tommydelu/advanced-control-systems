%% COMPUTE RNE FORMULATION AND COMPUTE THE REGRESSOR Y FOR THE LINEAR-IN-THE-DYNAMIC-PARAMETERS MODEL AND ESTIMATE THE PARAMETERS
clear; close all; clc;
syms a alpha d theta q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real

A = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha)  a*cos(theta);
    sin(theta) cos(theta)*cos(alpha)  -cos(theta)*sin(alpha) a*sin(theta); 
    0          sin(alpha)             cos(alpha)             d;
    0          0                      0                      1];

AB_0 = [0  0  1     0;    
    0 -1  0     0;
    1  0  0  0.15;
    0  0  0     1];

A3_EE = [0 0 1 0;
    1 0 0 0;
    0 1 0 0;
    0 0 0 1];

dh_table = [  0  pi/2  0.5+q1   0;  
    0.4 -pi/2       0  q2;  
    0.4     0       0  q3]; 

A0_1 = subs(A,[a alpha d theta],dh_table(1,:));
A1_2 = subs(A,[a alpha d theta],dh_table(2,:));
A2_3 = subs(A,[a alpha d theta],dh_table(3,:));
AB_EE = AB_0*A0_1*A1_2*A2_3*A3_EE;
AB_1 = AB_0*A0_1;
AB_2 = AB_1*A1_2;
AB_3 = AB_2*A2_3;

%% DATA FOR THE RECURSION

% Initial conditions (base frame does not move)
w_B   = [0; 0; 0];
wB_d  = [0; 0; 0];
pB_dd = [0; 0; 9.81];

z0 = [0; 0; 1];                 % joint axis in its own frame (steps 2 and 3)
z0_B = AB_0(1:3,3);             % joint 1 axis in the base frame = [1;0;0] (step 1)

% Rotations between consecutive frames
R1_B = (AB_1(1:3,1:3)).'; % base -> frame 1
R1_2 = A1_2(1:3,1:3);  R2_1 = R1_2.';
R2_3 = A2_3(1:3,1:3);  R3_2 = R2_3.';

% Vector from origin i-1 to origin i, in frame i
r1_01 = [0; 0.5 + q1; 0];
r2_12 = [0.4; 0; 0];
r3_23 = [0.4; 0; 0];

% Vector from origin i to the CoM of link i, in frame i
r1_1c1 = [0; -0.25; 0];
r2_2c2 = [-0.2; 0; 0];
r3_3c3 = [-0.2; 0; 0];

% Motors: no data in the URDF (symbolic, or set to 0 to ignore them)
syms kr1 kr2 kr3 real
z0_m1 = z0_B;  z1_m2 = z0;  z2_m3 = z0;   % rotor axes = joint axes

%% FORWARD EQUATIONS
% Link 1 (prismatic): input vectors are in the base frame
w1_1     = R1_B * w_B;
w1_1_d   = R1_B * wB_d;
p1_dd    = R1_B*(pB_dd + ddq1*z0_B) + 2*dq1*cross(w1_1, R1_B*z0_B) ...
    + cross(w1_1_d, r1_01) + cross(w1_1, cross(w1_1, r1_01));
p1_c1_dd = p1_dd + cross(w1_1_d, r1_1c1) + cross(w1_1, cross(w1_1, r1_1c1));
w0_m1_d  = wB_d + kr1*ddq1*z0_m1 + kr1*dq1*cross(w_B, z0_m1);

% Link 2 (revolute)
w2_2     = R2_1 * (w1_1 + dq2*z0);
w2_2_d   = R2_1 * (w1_1_d + ddq2*z0 + dq2*cross(w1_1, z0));
p2_dd    = R2_1*p1_dd + cross(w2_2_d, r2_12) + cross(w2_2, cross(w2_2, r2_12));
p2_c2_dd = p2_dd + cross(w2_2_d, r2_2c2) + cross(w2_2, cross(w2_2, r2_2c2));
w1_m2_d  = w1_1_d + kr2*ddq2*z1_m2 + kr2*dq2*cross(w1_1, z1_m2);

% Link 3 (revolute)
w3_3     = R3_2 * (w2_2 + dq3*z0);
w3_3_d   = R3_2 * (w2_2_d + ddq3*z0 + dq3*cross(w2_2, z0));
p3_dd    = R3_2*p2_dd + cross(w3_3_d, r3_23) + cross(w3_3, cross(w3_3, r3_23));
p3_c3_dd = p3_dd + cross(w3_3_d, r3_3c3) + cross(w3_3, cross(w3_3, r3_3c3));
w2_m3_d  = w2_2_d + kr3*ddq3*z2_m3 + kr3*dq3*cross(w2_2, z2_m3);

%% DYNAMIC PARAMETERS (symbolic)
syms m1 m2 m3 real                    % masses of the links
syms Im1 Im2 Im3 real                 % inertia of the rotors (around their axis)

% Inertia of each link w.r.t. its CoM, expressed in frame i (slide 9 of S3: minus signs out of the diagonal)
syms I1xx I1yy I1zz I1xy I1xz I1yz real
syms I2xx I2yy I2zz I2xy I2xz I2yz real
syms I3xx I3yy I3zz I3xy I3xz I3yz real
I1 = [ I1xx -I1xy -I1xz; -I1xy  I1yy -I1yz; -I1xz -I1yz  I1zz];
I2 = [ I2xx -I2xy -I2xz; -I2xy  I2yy -I2yz; -I2xz -I2yz  I2zz];
I3 = [ I3xx -I3xy -I3xz; -I3xy  I3yy -I3yz; -I3xz -I3yz  I3zz];

%% BACKWARD EQUATIONS
% Final conditions: no force and no moment from the environment on the EE (free motion)
f4_4  = [0; 0; 0];
mu4_4 = [0; 0; 0];
R3_4  = A3_EE(1:3,1:3);               % frame 3 -> EE frame (multiplies zero vectors)

% Link 3 (no motor after it)
f3_3  = R3_4*f4_4 + m3*p3_c3_dd;
mu3_3 = -cross(f3_3, r3_23 + r3_3c3) + R3_4*mu4_4 + cross(R3_4*f4_4, r3_3c3) ...
    + I3*w3_3_d + cross(w3_3, I3*w3_3);

% Link 2 (it carries the motor of joint 3)
f2_2  = R2_3*f3_3 + m2*p2_c2_dd;
mu2_2 = -cross(f2_2, r2_12 + r2_2c2) + R2_3*mu3_3 + cross(R2_3*f3_3, r2_2c2) ...
    + I2*w2_2_d + cross(w2_2, I2*w2_2) ...
    + kr3*ddq3*Im3*z2_m3 + kr3*dq3*Im3*cross(w2_2, z2_m3);

% Link 1 (it carries the motor of joint 2)
f1_1  = R1_2*f2_2 + m1*p1_c1_dd;
mu1_1 = -cross(f1_1, r1_01 + r1_1c1) + R1_2*mu2_2 + cross(R1_2*f2_2, r1_1c1) ...
    + I1*w1_1_d + cross(w1_1, I1*w1_1) ...
    + kr2*ddq2*Im2*z1_m2 + kr2*dq2*Im2*cross(w1_1, z1_m2);

%% JOINT TORQUES (slide 38, line 9, without friction)
% Joint axis z_{i-1} written in frame i
tau1 = f1_1.'  * (R1_B*z0_B) + kr1*Im1*(w0_m1_d.'*z0_m1);   % prismatic: force along the axis
tau2 = mu2_2.' * (R2_1*z0)   + kr2*Im2*(w1_m2_d.'*z1_m2);   % revolute: moment along the axis
tau3 = mu3_3.' * (R3_2*z0)   + kr3*Im3*(w2_m3_d.'*z2_m3);

tau_NE = simplify([tau1; tau2; tau3]);

