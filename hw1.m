clear; close all; clc;

robot = importrobot("00001_PRR_xyx TOMMASO DELUCCHI.urdf"); % Load URDF
showdetails(robot);
config = homeConfiguration(robot);                          % Create robot config with joint values = 0
figure;
show(robot,config);                                         % Display robot in home config

% Try with different joint values (each joint has lower and upper limit)
config(1).JointPosition = -0.4; % [-0.5,0]
config(2).JointPosition = 1.2;  % [-3.14, 3.14]
config(3).JointPosition = 0;    % [-3.14, 3.14]

if (config(1).JointPosition > 0 || config(1).JointPosition < -0.5 || config(2).JointPosition > 3.14 || config(2).JointPosition < -3.14 || config(3).JointPosition > 3.14 || config(3).JointPosition < -3.14)
     error("One of the joints exceeds its limits. Try again");
else
    figure;
    show(robot,config);
    fprintf("CONFIGURATION DISPLAYED" + ...
            "\nValue of q1: %.2f\n" + ...
            "Value of q2: %.2f\n" + ...
            "Value of q3: %.2f\n" + ...
            "--------------------\n", config(1).JointPosition,config(2).JointPosition,config(3).JointPosition);
end



%% DIRECT KINEMATIC: SOLUTION BY HAND AND BY USING ROBOTICS TOOLBOX

% This figure shows how I placed the frames following DH standard
figure;
dh_frames = imread("dh_frames.HEIC");
imshow(rot90(dh_frames)); % Just to straighten the photo

syms a alpha d theta q1 q2 q3 real

% Define the general transformation between two frames
A = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha)  a*cos(theta);
     sin(theta) cos(theta)*cos(alpha)  -cos(theta)*sin(alpha) a*sin(theta); 
     0          sin(alpha)             cos(alpha)             d;
     0          0                      0                      1];

% Define constant transformation that brings a point in frame 0 to frame base
AB_0 = [0  0  1     0;    
        0 -1  0     0;
        1  0  0  0.15;
        0  0  0     1];

% Define constant transformation that brings a point in nsa to frame 3
A3_EE = [0 0 1 0;
         1 0 0 0;
         0 1 0 0;
         0 0 0 1];


% Build the DH table
dh_table = [  0  pi/2  0.5+q1   0;  % frame 0-1
            0.4 -pi/2       0  q2;  % frame 1-2
            0.4     0       0  q3]; % frame 2-3

% Build the intermediate transformations
A0_1 = subs(A,[a alpha d theta],dh_table(1,:));
A1_2 = subs(A,[a alpha d theta],dh_table(2,:));
A2_3 = subs(A,[a alpha d theta],dh_table(3,:));

% Concatenate everything to obtain AB_EE
AB_EE = AB_0*A0_1*A1_2*A2_3*A3_EE;
fprintf("The symbolic transformation from base to EE is:\n\n")
disp(AB_EE)
fprintf("--------------------\n")

% For the Jacobian I also need to build the transformation of each frame w.r.t the base
AB_1 = AB_0*A0_1;
AB_2 = AB_1*A1_2;
AB_3 = AB_2*A2_3;

% Insert the values of the joints into AB_EE
AB_EE_num = double(subs(AB_EE, [q1 q2 q3], [config(1).JointPosition config(2).JointPosition config(3).JointPosition]));
fprintf("The numerical transformation from base to EE is:\n\n")
disp(AB_EE_num)
fprintf("--------------------\n")

% Now let's check that MATLAB's built-in DK gives the same result
AB_EE_auto = getTransform(robot, config,'ee','base_link'); % base_link is the target, we want the ee expressed in the base link frame
Rz = @(angle) [cos(angle) -sin(angle) 0;
               sin(angle)  cos(angle) 0;
               0          0          1];
AB_EE_auto(1:3,1:3) = AB_EE_auto(1:3,1:3)*Rz(-pi/2); % To also consider the nsa tf
fprintf("The numerical transformation from base to EE computed automaticaly is:\n\n")
disp(AB_EE_auto)
fprintf("--------------------\n")

% Compute the difference between the manual and automatic method
err = round((AB_EE_num - AB_EE_auto),10);
fprintf("The difference between the two methods is\n")
disp(err)
fprintf("--------------------\n")



%% INVERSE KINEMATIC: SOLUTION BY HAND AND BY USING ROBOTICS TOOLBOX
desired_pose = [0.2 0.2 0.2];
x = desired_pose(1); y = desired_pose(2); z = desired_pose(3);

fprintf("EE position from DK (used to solve IK):\n")
disp(simplify(AB_EE(1:3,4)))
fprintf("--------------------\n")

% q3 from y = -0.4*sin(q3)
s3 = -y/0.4;
if abs(s3) > 1
    error("Pose not reachable: |sin(q3)| > 1");
end
c3_a =  sqrt(1 - s3^2);   % elbow solution A (cos(q3) > 0)
c3_b = -sqrt(1 - s3^2);   % elbow solution B (cos(q3) < 0)
q3_a = atan2(s3, c3_a);
q3_b = atan2(s3, c3_b);

% q2 from z = 0.15 + 0.4*cos(q2)*(1+cos(q3))
r_a = 0.4*(1 + c3_a);
r_b = 0.4*(1 + c3_b);
c2_a = (z - 0.15)/r_a;
c2_b = (z - 0.15)/r_b;
q2_a1 = atan2( sqrt(1 - c2_a^2), c2_a);
q2_a2 = atan2(-sqrt(1 - c2_a^2), c2_a);
q2_b1 = atan2( sqrt(1 - c2_b^2), c2_b);
q2_b2 = atan2(-sqrt(1 - c2_b^2), c2_b);

% q1 from x = 0.5 + q1 + 0.4*sin(q2)*(1+cos(q3))
q1_a1 = x - 0.5 - r_a*sin(q2_a1);
q1_a2 = x - 0.5 - r_a*sin(q2_a2);
q1_b1 = x - 0.5 - r_b*sin(q2_b1);
q1_b2 = x - 0.5 - r_b*sin(q2_b2);

% All 4 solutions, one per row: [q1 q2 q3]
q_sols = [q1_a1 q2_a1 q3_a;
          q1_a2 q2_a2 q3_a;
          q1_b1 q2_b1 q3_b;
          q1_b2 q2_b2 q3_b];
fprintf("ALL IK SOLUTIONS [q1 q2 q3]:\n")
disp(q_sols)

% Keep only the solutions inside the limit of q1 (q2 and q3 are always inside [-pi,pi] with atan2)
q1_min = -0.5; q1_max = 0.0;
feasible = q_sols(:,1) >= q1_min & q_sols(:,1) <= q1_max;
fprintf("FEASIBLE IK SOLUTIONS [q1 q2 q3]:\n")
disp(q_sols(feasible,:))
fprintf("--------------------\n")

% Cross check: put the IK solutions inside DK
pos_b1 = double(subs(AB_EE(1:3,4), [q1 q2 q3], q_sols(3,:))); % Put feasible solutions here
pos_b2 = double(subs(AB_EE(1:3,4), [q1 q2 q3], q_sols(4,:)));
fprintf("Error between DK of solution b1 and desired position:\n")
disp(pos_b1 - desired_pose(:))
fprintf("Error between DK of solution b2 and desired position:\n")
disp(pos_b2 - desired_pose(:))

config(1).JointPosition = q_sols(3,1);
config(2).JointPosition = q_sols(3,2);
config(3).JointPosition = q_sols(3,3);
figure;
show(robot,config);
fig = gcf;
exportgraphics(fig, "imgs/5.png");
fprintf("--------------------\n")

% IK with Robotics toolbox: see https://it.mathworks.com/help/robotics/ref/inversekinematics-system-object.html
ik = inverseKinematics("RigidBodyTree",robot);
pose_tf = [1 0 0 desired_pose(1);
           0 1 0 desired_pose(2);
           0 0 1 desired_pose(3);
           0 0 0 1];
w = [0 0 0 1 1 1];
initialGuess = homeConfiguration(robot);
[configSol, solInfo] = ik('ee', pose_tf, w, initialGuess);
fprintf("Joint values obtained automaticaly\n")
disp(configSol(1).JointPosition); disp(configSol(2).JointPosition); disp(configSol(3).JointPosition);

q_sols_toolbox = [configSol(1).JointPosition, configSol(2).JointPosition, configSol(3).JointPosition];

% Compare toolbox solution with the manual ones
diff_sols = q_sols - q_sols_toolbox;                   % one row for each manual solution
diff_sols(:,2:3) = wrapToPi(diff_sols(:,2:3));         % q2 and q3 are angles: 2pi difference = same angle
[~, idx] = min(vecnorm(diff_sols, 2, 2));              % row closest to the toolbox solution
fprintf("Toolbox solution (initial guess = home) corresponds to manual solution number %d\n", idx)
fprintf("Difference with that manual solution:\n")
disp(diff_sols(idx,:))
fprintf("--------------------\n")

% Now start from a different initial guess, near the other feasible solution
initialGuess2 = homeConfiguration(robot);
initialGuess2(1).JointPosition = -0.32;
initialGuess2(2).JointPosition = 0.4;   % negative q2 instead of positive
initialGuess2(3).JointPosition = -2.6;
[configSol2, solInfo2] = ik('ee', pose_tf, w, initialGuess2);
q_sols_toolbox2 = [configSol2(1).JointPosition, configSol2(2).JointPosition, configSol2(3).JointPosition];

diff_sols2 = q_sols - q_sols_toolbox2;
diff_sols2(:,2:3) = wrapToPi(diff_sols2(:,2:3));
[~, idx2] = min(vecnorm(diff_sols2, 2, 2));
fprintf("Toolbox solution (new initial guess) corresponds to manual solution number %d\n", idx2)
fprintf("Difference with that manual solution:\n")
disp(diff_sols2(idx2,:))
fprintf("--------------------\n")

%% ANALYTICAL AND GEOMETRIC JACOBIAN: SOLUTION BY HAND AND BY USING THE TOOLBOX

eq1 = simplify(AB_EE(1,4));  eq2 = simplify(AB_EE(2,4));  eq3 = simplify(AB_EE(3,4));

% Orientation of the EE with Euler angles ZYX: R = Rz(yaw)*Ry(pitch)*Rx(roll)
R = simplify(AB_EE(1:3,1:3));
yaw   = simplify(atan2(R(2,1), R(1,1)));
pitch = simplify(atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2)));
roll  = simplify(atan2(R(3,2), R(3,3)));

JA_man = simplify([diff(eq1,q1),   diff(eq1,q2),   diff(eq1,q3);
                   diff(eq2,q1),   diff(eq2,q2),   diff(eq2,q3);
                   diff(eq3,q1),   diff(eq3,q2),   diff(eq3,q3);
                   diff(yaw,q1),   diff(yaw,q2),   diff(yaw,q3);
                   diff(pitch,q1), diff(pitch,q2), diff(pitch,q3);
                   diff(roll,q1),  diff(roll,q2),  diff(roll,q3)]);

z0 = AB_0(1:3,3);  z1 = AB_1(1:3,3);  z2 = AB_2(1:3,3);

pEE = AB_EE(1:3,4);  p0 = AB_0(1:3,4);  p1 = AB_1(1:3,4);  p2 = AB_2(1:3,4);

JG_man = simplify([z0         cross(z1,pEE-p1)  cross(z2,pEE-p2);
          zeros(3,1)       z1                z2]);

fprintf("ANALYITICAL AND GEOMETRIC JACOBIANS FOUND MANUALLY\n")
disp(JA_man)
fprintf("\n")
disp(JG_man)

JA_man_num  = double(subs(JA_man,  [q1 q2 q3], [0 0 0]));
JG_man_num = subs(JG_man,[q1 q2 q3], [0 0 0]);

% Find JA and JG with the toolbox
JA_auto = simplify(jacobian([eq1, eq2, eq3, yaw, pitch, roll], [q1 q2 q3]));
JG_auto = geometricJacobian(robot,homeConfiguration(robot),'ee');
JG_auto = [JG_auto(4:6,:); JG_auto(1:3,:)]; % Because MATLAB's geometricJacobian function uses first w and then v

JA_auto_num = double(subs(JA_auto, [q1 q2 q3], [0 0 0]));
fprintf("--------------------\n\n")
fprintf("ERROR BETWEEN MANUAL AND AUTO GEOMETRIC JACOBIAN\n")
disp(JG_man_num - JG_auto);
fprintf("--------------------\n")
fprintf("ERROR BETWEEN MANUAL AND AUTO ANALYTICAL JACOBIAN\n")
disp(JA_man_num - JA_auto_num)
fprintf("--------------------\n")

%% LINK BETWEEN ANALYTICAL AND GEOMETRIC JACOBIAN

% Matrix that transforms the derivatives of the ZYX Euler angles into the angular velocity
T = [0  -sin(yaw)  cos(yaw)*cos(pitch);
     0   cos(yaw)  sin(yaw)*cos(pitch);
     1   0        -sin(pitch)];
TA = [eye(3) zeros(3); zeros(3) T];

JG_from_JA = simplify(TA*JA_man);

% Check in the home configuration
fprintf("ERROR BETWEEN JG AND TA*JA IN HOME CONFIGURATION\n")
disp(double(subs(JG_from_JA - JG_man, [q1 q2 q3], [0 0 0])))

% Check in another configuration
fprintf("ERROR BETWEEN JG AND TA*JA IN q = [-0.2 0.3 0.5]\n")
disp(double(subs(JG_from_JA - JG_man, [q1 q2 q3], [-0.2 0.3 0.5])))
fprintf("--------------------\n")
