%% COMPUTE KINETIC AND POTENTIAL ENERGY OF THE 3-DOF MANIPULATOR (RIGID BODY)

clearvars -except robot AB_0 AB_1 AB_2 AB_3 q1 q2 q3; close all; clc

debug = 1;
if debug
    config = homeConfiguration(robot);                         
    figure;
    show(robot,config);              
    figure;
    dh_frames = imread("dh_frames.HEIC");
    imshow(rot90(dh_frames)); 
end


%% LINKS INERTIA MATRICES (used in rotational part of the kinetic energy)

m1 = 1; m2 = 1; m3 = 1; % mass of links (see URDF)
w = 0.036;  L1 = 0.5;   % link 1: box w x w x L1, long side along y1
r = 0.019;  h = 0.4;    % links 2 and 3: cylinders, axis along x2 and x3

% Inertia w.r.t. the CoM, expressed with the axes of frame i (constant and diagonal)
I1_1 = diag([1/12*m1*(L1^2 + w^2), 1/12*m1*(w^2 + w^2), 1/12*m1*(w^2 + L1^2)]); % Link 1- Look at frame 1
I2_2 = diag([1/2*m2*r^2, 1/12*m2*(3*r^2 + h^2), 1/12*m2*(3*r^2 + h^2)]);        % Link 2- Look at frame 2
I3_3 = diag([1/2*m3*r^2, 1/12*m3*(3*r^2 + h^2), 1/12*m3*(3*r^2 + h^2)]);        % Link 3- Look at frame 3

% Rotation of each link frame w.r.t. the base.
% Use these to transform inertia w.r.t base (see slides)
R1 = AB_1(1:3,1:3);  R2 = AB_2(1:3,1:3);  R3 = AB_3(1:3,1:3);
fprintf("Inertia of link 1:\n")
disp(I1_1)
fprintf("Inertia of link 2:\n")
disp(I2_2)
fprintf("Inertia of link 3:\n\n")
disp(I3_3)

%% POSITION OF THE CoM OF EACH LINK W.R.T. THE BASE

% CoM of each link expressed in its own frame (constant)
pc1_1 = [0; -L1/2; 0; 1];
pc2_2 = [-h/2; 0; 0; 1];
pc3_3 = [-h/2; 0; 0; 1];

% CoM of each link expressed in the base frame (depends on q). Used in the jacobian
pc1 = AB_1*pc1_1;  pc1 = simplify(pc1(1:3)); % this one depends only on q1
pc2 = AB_2*pc2_2;  pc2 = simplify(pc2(1:3)); % this one depends on q1 and q2
pc3 = AB_3*pc3_3;  pc3 = simplify(pc3(1:3)); % depends on q1 q2 q3

fprintf("POSITION OF THE CoM OF EACH LINK\n")
disp(pc1); disp(pc2); disp(pc3)
fprintf("--------------------\n")

%% PARTIAL JACOBIANS OF THE CoM OF EACH LINK

z0 = AB_0(1:3,3);  z1 = AB_1(1:3,3);  z2 = AB_2(1:3,3);
p1 = AB_1(1:3,4);  p2 = AB_2(1:3,4);

% Link 1: only joint 1 (prismatic) moves it
JP1 = [z0, zeros(3,1), zeros(3,1)];
JO1 = zeros(3,3);

% Link 2: joints 1 and 2
JP2 = simplify([z0, cross(z1, pc2 - p1), zeros(3,1)]);
JO2 = [zeros(3,1), z1, zeros(3,1)];

% Link 3: all the joints
JP3 = simplify([z0, cross(z1, pc3 - p1), cross(z2, pc3 - p2)]);
JO3 = [zeros(3,1), z1, z2];

% Check: the linear part must be equal to the derivative of the CoM position
fprintf("ERROR BETWEEN PARTIAL JACOBIANS AND DERIVATIVE OF CoM POSITION\n")
disp(simplify(JP1 - jacobian(pc1, [q1 q2 q3])))
disp(simplify(JP2 - jacobian(pc2, [q1 q2 q3])))
disp(simplify(JP3 - jacobian(pc3, [q1 q2 q3])))
fprintf("--------------------\n")

%% KINETIC ENERGY

syms dq1 dq2 dq3 real
dq = [dq1; dq2; dq3];

% Contribution of each link to the inertia matrix B(q): translational + rotational part. It's the summatory part.
B1 = m1*(JP1.')*JP1 + (JO1.')*R1*I1_1*(R1.')*JO1;
B2 = m2*(JP2.')*JP2 + (JO2.')*R2*I2_2*(R2.')*JO2;
B3 = m3*(JP3.')*JP3 + (JO3.')*R3*I3_3*(R3.')*JO3;

B = simplify(B1 + B2 + B3);
T_kin = simplify(1/2*(dq.')*B*dq);

fprintf("INERTIA MATRIX B(q)\n")
disp(vpa(B,4))
fprintf("KINETIC ENERGY\n")
disp(vpa(T_kin,4))
fprintf("--------------------\n")

%% Let's check that the inertial matrix computed is symmetric, positive definite and configuration-dependent

fprintf("Size of B matrix: %d x %d\n", size(B,1), size(B,2))
if (size(B,1) ~= size(B,2))
    error("B is not a square matrix!")
end

if ~isequal(simplify(B - B.'), sym(zeros(3)))
    error("B is not symmetric!")
end

B_test1 = double(subs(B,[q1 q2 q3], [0 0 0]));
B_test2 = double(subs(B,[q1 q2 q3], [-0.2 0.3 0.5]));
B_test3 = double(subs(B,[q1 q2 q3], [-0.4 1.2 -2]));

eig_values_test1 = eig(B_test1);
fprintf("Eigenvalues of B in q = [0 0 0]:\n"); disp(eig_values_test1)
if any(eig_values_test1 <= 0)
    error("B is not positive definite!")
end

eig_values_test2 = eig(B_test2);
fprintf("Eigenvalues of B in q = [-0.2 0.3 0.5]:\n"); disp(eig_values_test2)
if any(eig_values_test2 <= 0)
    error("B is not positive definite!")
end

eig_values_test3 = eig(B_test3);
fprintf("Eigenvalues of B in q = [-0.4 1.2 -2]:\n"); disp(eig_values_test3)
if any(eig_values_test3 <= 0)
    error("B is not positive definite!")
end

if isequal(B_test1, B_test2)
    error("B does not depend on the configuration!")
end
fprintf("B is symmetric, positive definite and configuration-dependent\n\n")
%% Check that if v = [0 0 0] then T_kin = 0
T_kin_zero_vels = double(subs(T_kin,[dq1 dq2 dq3], [0 0 0]));

if (T_kin_zero_vels ~= 0)
    error("Kinetic energy is not zero when it should be!")
end

%% POTENTIAL ENERGY (see slide 20)

g0 = [0; 0; -9.81];   % gravity in the base frame (column vector)

U1 = -m1*(g0.')*pc1;
U2 = -m2*(g0.')*pc2;
U3 = -m3*(g0.')*pc3;
U = simplify(U1 + U2 + U3);

fprintf("POTENTIAL ENERGY\n")
disp(vpa(U,4))
fprintf("--------------------\n")

fprintf("Total energy of the manipulator\n")
disp(vpa(T_kin + U, 4))


%% Steiner theorem: use this to compute Inertia matrix

% Scegli t_start and t_end e fai polinomio quinto per vedere come varia
% T_kin e U e fai vedere movimento robot non so come si fa