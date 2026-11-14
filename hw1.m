%% HW1: FIND DH TABLE, COMPUTE DIRECT KINEMATICS, INVERSE KINEMATICS, JACOBIANS (ANALYTICAL AND GEOMETRIC)

clear; close all; clc;

save_figs = false;
fig_name = 1;

robot = importrobot("00001_PRR_xyx TOMMASO DELUCCHI.urdf"); % Load URDF
showdetails(robot);
config = homeConfiguration(robot);                          % Create robot config with joint values = 0
figure("Name","HOME CONFIGURATION");
show(robot,config);                                         % Display robot in home config

% Try with different joint values (each joint has lower and upper limit)
config(1).JointPosition = -0.4; % [-0.5,0]
config(2).JointPosition = 3.14;  % [-3.14, 3.14]
config(3).JointPosition = pi/2;    % [-3.14, 3.14]

% I do not check revolute limits because they can always be brought inside [-3.14,3.14] by adding 2*pi (collisions are neglected)
if (config(1).JointPosition > 0 || config(1).JointPosition < -0.5)
    error("Prismatic joint exceeds its limits. Try with different values");
else
    figure("Name","USER SELECTED CONFIGURATION")
    show(robot,config);
    fprintf("USER SELECTED CONFIGURATION VALUES" + ...
        "\nValue of q1: %.2f\n" + ...
        "Value of q2: %.2f\n" + ...
        "Value of q3: %.2f\n" + ...
        "--------------------\n", config(1).JointPosition,config(2).JointPosition,config(3).JointPosition);
end


%%  DIRECT KINEMATIC: SOLUTION BY HAND AND BY USING ROBOTICS TOOLBOX

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
AB_0 = [0   0  1        0;    
        0  -1  0        0;
        1   0  0     0.15;
        0   0  0        1];

% Define constant transformation that brings a point in n-s-a frame to frame 3
A3_EE = [0 0 1 0;
         1 0 0 0;
         0 1 0 0;
         0 0 0 1];

% Build the DH table
dh_table = [  0     pi/2  0.5+q1   0;  % frame 0-1: depends on q1 (link 1 translates w.r.t z0, frame 1 moves with link 1)
              0.4  -pi/2       0  q2;  % frame 1-2: depends on q2 (link 2 rotates w.r.t z1, frame 2 moves with link 2)
              0.4      0       0  q3]; % frame 2-3: depends on q3 (link 3 rotates w.r.t z2, frame 3 moves with link 3)

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

fprintf("The symbolic transformation from base to frame 1 is:\n\n")
disp(AB_1)
fprintf("--------------------\n")

fprintf("The symbolic transformation from base to frame 2 is:\n\n")
disp(AB_2)
fprintf("--------------------\n")

fprintf("The symbolic transformation from base to frame 3 is:\n\n")
disp(AB_3)
fprintf("--------------------\n")

% NUMERICAL EXAMPLE: Insert the values of the joints into AB_EE
AB_EE_num = double(subs(AB_EE, [q1 q2 q3], [config(1).JointPosition config(2).JointPosition config(3).JointPosition]));
fprintf("EXAMPLE - The numerical transformation from base to EE using q1="+config(1).JointPosition+" q2="+config(2).JointPosition+" q3="+config(3).JointPosition+" is:\n")
disp(vpa(AB_EE_num,3))
fprintf("--------------------\n")

% Now let's check that MATLAB's built-in DK gives the same result
AB_EE_auto = getTransform(robot, config,'ee','base_link'); % base_link is the target, we want the ee expressed in the base link frame
Rz = @(angle) [cos(angle) -sin(angle) 0;
    sin(angle)  cos(angle) 0;
    0          0          1];
AB_EE_auto(1:3,1:3) = AB_EE_auto(1:3,1:3)*Rz(-pi/2); % To also consider the n-s-a tf
fprintf("The numerical transformation from base to EE computed automaticaly is:\n\n")
disp(AB_EE_auto)
fprintf("--------------------\n")

% Compute the difference between the manual and automatic method
dk_err = round((AB_EE_num - AB_EE_auto),10);
fprintf("DIRECT KINEMATIC: difference between manual and automatic method\n")
disp(dk_err)
fprintf("--------------------\n")


%% INVERSE KINEMATIC: SOLUTION BY HAND AND BY USING ROBOTICS TOOLBOX

desired_pose = [0.2 0.2 0.2];
x = desired_pose(1);
y = desired_pose(2);
z = desired_pose(3);

fprintf("The symbolic position of the End Effector is:\n\n")
disp(simplify(AB_EE(1:3,4)))
fprintf("--------------------\n")

% How to find q3
s3 = -y/0.4;
if abs(s3) > 1
    error("Pose not reachable: |sin(q3)| > 1");
end

c3_a =  sqrt(1 - s3^2);   % elbow solution A (cos(q3) > 0)
c3_b = -sqrt(1 - s3^2);   % elbow solution B (cos(q3) < 0)
q3_a = atan2(s3, c3_a);
q3_b = atan2(s3, c3_b);

% How to find q2 (from equation 3) ) - BRANCH A
r_a = 0.4*(1 + c3_a);
c2_a = (z - 0.15)/r_a;
if abs(c2_a) > 1
    error("Pose not reachable: |cos(q2)| > 1");
end
q2_a1 = atan2( sqrt(1 - c2_a^2), c2_a);
q2_a2 = atan2(-sqrt(1 - c2_a^2), c2_a);
q1_a1 = x - 0.5 - r_a*sin(q2_a1);
q1_a2 = x - 0.5 - r_a*sin(q2_a2);

% BRANCH B
r_b = 0.4*(1 + c3_b);
c2_b = (z - 0.15)/r_b;
if abs(c2_b) <= 1
    q2_b1 = atan2( sqrt(1 - c2_b^2), c2_b);
    q2_b2 = atan2(-sqrt(1 - c2_b^2), c2_b);
    q1_b1 = x - 0.5 - r_b*sin(q2_b1);
    q1_b2 = x - 0.5 - r_b*sin(q2_b2);
else
    fprintf("Branch b does not exist for this pose (|cos(q2)| > 1)\n")
    q2_b1 = NaN; q2_b2 = NaN; q1_b1 = NaN; q1_b2 = NaN;
end

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

j = 0;
feasible_sols = q_sols(feasible,:);
for i = 1:size(q_sols(feasible,:),1)
    figure("Name","INVERSE KINEMATICS - FEASIBLE SOLUTION "+i)
    config(1).JointPosition = feasible_sols(i,1);
    config(2).JointPosition = feasible_sols(i,2);
    config(3).JointPosition = feasible_sols(i,3);
    show(robot,config);
    if save_figs
        fig = gcf;
        exportgraphics(fig, "imgs_hw1/"+fig_name+j+".png");
        j = j+1;
    end
end

% Cross check: put the IK solutions inside DK
for i = 1:size(q_sols(feasible,:),1)
    solutionPosition = double(subs(AB_EE(1:3,4), [q1 q2 q3], feasible_sols(i,:)));
    fprintf("Position error for feasible solution %d:\n", i)
    disp(solutionPosition - desired_pose(:))
end
fprintf("--------------------\n")

% IK with Robotics toolbox: see https://it.mathworks.com/help/robotics/ref/inversekinematics-system-object.html
ik = inverseKinematics("RigidBodyTree",robot);
pose_tf = [1 0 0 x;
    0 1 0 y;
    0 0 1 z;
    0 0 0 1];
w = [0 0 0 1 1 1];
initialGuess = homeConfiguration(robot);
[configSol, solInfo] = ik('ee', pose_tf, w, initialGuess);
fprintf("Joint values obtained automaticaly\n")
q_sol_toolbox = [configSol(1).JointPosition, configSol(2).JointPosition, configSol(3).JointPosition];
disp(q_sol_toolbox);
fprintf("--------------------\n")

% Compare toolbox solution with the manual ones
diff_sols = q_sols - q_sol_toolbox;                    % one row for each manual solution
diff_sols(:,2:3) = wrapToPi(diff_sols(:,2:3));         % q2 and q3 are angles: 2pi difference = same angle
[~, idx] = min(vecnorm(diff_sols, 2, 2));              % row closest to the toolbox solution (norm2 used)

fprintf("Toolbox solution (initial guess = home) corresponds to manual solution number %d\n", idx)
fprintf("Difference with that manual solution:\n")
disp(diff_sols(idx,:))
fprintf("--------------------\n")

% Now start from a different initial guess, near the other feasible solution
% CHANGE VALUES ACCORDING TO THE SOLUTION YOU WANT TO FIND
initialGuess2 = homeConfiguration(robot);
initialGuess2(1).JointPosition = -0.32;
initialGuess2(2).JointPosition = 0.4;   
initialGuess2(3).JointPosition = -2.6;

[configSol2, solInfo2] = ik('ee', pose_tf, w, initialGuess2);
q_sol_toolbox2 = [configSol2(1).JointPosition, configSol2(2).JointPosition, configSol2(3).JointPosition];

diff_sols2 = q_sols - q_sol_toolbox2;
diff_sols2(:,2:3) = wrapToPi(diff_sols2(:,2:3));
[~, idx2] = min(vecnorm(diff_sols2, 2, 2));
fprintf("Toolbox solution (new initial guess) corresponds to manual solution number %d\n", idx2)
fprintf("Difference with that manual solution:\n")
disp(diff_sols2(idx2,:))
fprintf("--------------------\n")

%% ANALYTICAL AND GEOMETRIC JACOBIAN: SOLUTION BY HAND AND BY USING THE TOOLBOX

eq1 = simplify(AB_EE(1,4));  eq2 = simplify(AB_EE(2,4));  eq3 = simplify(AB_EE(3,4));
R = simplify(AB_EE(1:3,1:3));

% Orientation of the EE with Euler angles ZYX: R = Rz(yaw)*Ry(pitch)*Rx(roll)
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

JG_man = simplify([z0           cross(z1,pEE-p1)     cross(z2,pEE-p2);
                   zeros(3,1)         z1                    z2]);

fprintf("ANALYITICAL AND GEOMETRIC JACOBIANS FOUND MANUALLY\n")
disp(simplify(JA_man))
fprintf("\n")
disp(simplify(JG_man))

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

fprintf("Matrix that relates angular velocity and Euler angles\n")
fprintf("--------------------\n")
disp(TA)

JG_from_JA = simplify(TA*JA_man);

% Check in the home configuration
fprintf("ERROR BETWEEN JG AND TA*JA IN HOME CONFIGURATION\n")
disp(double(subs(JG_from_JA - JG_man, [q1 q2 q3], [0 0 0])))

% Check in another configuration
fprintf("ERROR BETWEEN JG AND TA*JA IN q = [-0.2 0.3 0.5]\n")
disp(double(subs(JG_from_JA - JG_man, [q1 q2 q3], [-0.2 0.3 0.5])))
fprintf("--------------------\n")
