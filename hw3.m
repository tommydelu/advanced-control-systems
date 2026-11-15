%% EQUATIONS OF MOTION (dynamic model)

clear; close all; clc;

robot = importrobot("00001_PRR_xyx TOMMASO DELUCCHI.urdf");

syms t Q1(t) Q2(t) Q3(t) tau1(t) tau2(t) tau3(t)
syms a alpha d theta q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real


A = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha)  a*cos(theta);
    sin(theta) cos(theta)*cos(alpha)  -cos(theta)*sin(alpha) a*sin(theta); 
    0          sin(alpha)             cos(alpha)             d;
    0          0                      0                      1];

AB_0 = [0   0  1        0;    
        0  -1  0        0;
        1   0  0     0.15;
        0   0  0        1];

A3_EE = [0 0 1 0;
         1 0 0 0;
         0 1 0 0;
         0 0 0 1];

dh_table = [  0     pi/2  0.5+q1   0;  
              0.4  -pi/2       0  q2;  
              0.4      0       0  q3]; 

A0_1  = subs(A,[a alpha d theta],dh_table(1,:));
A1_2  = subs(A,[a alpha d theta],dh_table(2,:));
A2_3  = subs(A,[a alpha d theta],dh_table(3,:));
AB_EE = AB_0*A0_1*A1_2*A2_3*A3_EE;
AB_1  = AB_0*A0_1;
AB_2  = AB_1*A1_2;
AB_3  = AB_2*A2_3;

% ------------------------------ HW2 ------------------------------ %

debug = 0;
if debug
    config = homeConfiguration(robot);                         
    figure;
    show(robot,config);              
    figure;
    dh_frames = imread("dh_frames.HEIC");
    imshow(rot90(dh_frames)); 
end

m1 = 1; m2 = 1; m3 = 1; 
w = 0.036;  L1 = 0.5;   
r = 0.019;  h = 0.4;    

I1_1 = diag([1/12*m1*(L1^2 + w^2), 1/12*m1*(w^2 + w^2), 1/12*m1*(w^2 + L1^2)]); 
I2_2 = diag([1/2*m2*r^2, 1/12*m2*(3*r^2 + h^2), 1/12*m2*(3*r^2 + h^2)]);        
I3_3 = diag([1/2*m3*r^2, 1/12*m3*(3*r^2 + h^2), 1/12*m3*(3*r^2 + h^2)]);        

I1_1_steiner = I1_1 + m1*([0 L1/2 0]*[0; L1/2; 0]*eye(3,3) - [0;L1/2;0]*[0 L1/2 0]);
I2_2_steiner = I2_2 + m2*([-0.2 0 0]*[-0.2; 0; 0]*eye(3,3) - [-0.2;0;0]*[-0.2 0 0]);
I3_3_steiner = I3_3 + m3*([-0.2 0 0]*[-0.2; 0; 0]*eye(3,3) - [-0.2;0;0]*[-0.2 0 0]);

R1 = AB_1(1:3,1:3);  R2 = AB_2(1:3,1:3);  R3 = AB_3(1:3,1:3);

pc1_1 = [0; -L1/2; 0; 1];
pc2_2 = [-h/2; 0; 0; 1];
pc3_3 = [-h/2; 0; 0; 1];

pc1 = AB_1*pc1_1;  pc1 = simplify(pc1(1:3)); 
pc2 = AB_2*pc2_2;  pc2 = simplify(pc2(1:3)); 
pc3 = AB_3*pc3_3;  pc3 = simplify(pc3(1:3)); 

z0 = AB_0(1:3,3);  z1 = AB_1(1:3,3);  z2 = AB_2(1:3,3);
p1 = AB_1(1:3,4);  p2 = AB_2(1:3,4);

JP1 = [z0, zeros(3,1), zeros(3,1)];
JO1 = zeros(3,3);
JP2 = simplify([z0, cross(z1, pc2 - p1), zeros(3,1)]);
JO2 = [zeros(3,1), z1, zeros(3,1)];
JP3 = simplify([z0, cross(z1, pc3 - p1), cross(z2, pc3 - p2)]);
JO3 = [zeros(3,1), z1, z2];

dq = [dq1; dq2; dq3];

B1_steiner = m1*(JP1.')*JP1 + (JO1.')*R1*I1_1_steiner*(R1.')*JO1;
B2_steiner = m2*(JP2.')*JP2 + (JO2.')*R2*I2_2_steiner*(R2.')*JO2;
B3_steiner = m3*(JP3.')*JP3 + (JO3.')*R3*I3_3_steiner*(R3.')*JO3;
B = simplify(B1_steiner + B2_steiner + B3_steiner);
T_kin = simplify(1/2*(dq.')*B*dq);

g0 = [0; 0; -9.81];   

U1 = -m1*(g0.')*pc1;
U2 = -m2*(g0.')*pc2;
U3 = -m3*(g0.')*pc3;
U = simplify(U1 + U2 + U3);

% ------------------------------ HW3 ------------------------------ %

q   = [q1; q2; q3];
Q   = [Q1; Q2; Q3];
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

% dL/dq'
term1 = B*dq;
fprintf("TERM 1:\n\n");
disp(vpa(term1,4))
fprintf("--------------------------------------\n\n")

% d/dt * dL/dq' - 1) Find B_d which is the derivative of B w.r.t time
B_d = diff(B_Q,t);
fprintf("Time derivative of B:\n\n");
disp(vpa(B_d,4))
fprintf("--------------------------------------\n\n")

% d/dt * B*q' = B'*q' + B*q''
term2 = simplify(B_Q*ddq + B_d*dq);
fprintf("TERM 2:\n\n");
disp(vpa(term2,4))
fprintf("--------------------------------------\n\n")

% dL / dq
term3 = simplify(0.5.* [diff(dq'*B*dq, q1).'; diff(dq'*B*dq, q2).'; diff(dq'*B*dq, q3).'] - [ diff(U,q1).'; diff(U,q2).'; diff(U,q3).']);
term3 = subs(term3, [q1 q2 q3], [Q1 Q2 Q3]);
fprintf("TERM 3:\n\n");
disp(vpa(term3,4))
fprintf("--------------------------------------\n\n")

% Equation of motion
equationOfMotion = simplify(term2 - term3); 
equationOfMotion = subs(equationOfMotion, [diff(Q1,t) diff(Q2,t) diff(Q3,t)], [dq1 dq2 dq3]);
equationOfMotion = subs(equationOfMotion, [Q1 Q2 Q3], [q1 q2 q3]);

fprintf("EQUATION OF MOTION:\n\n");
disp(vpa(equationOfMotion,4))
fprintf("--------------------------------------\n\n")


% EQUATION OF MOTION USING B - C - g

% Christoffel symbols to extract C

C = sym(zeros(3,3));
for i=1:3
    for j=1:3
        for k=1:3
            cijk = 0.5 * ( diff(B(i,j),q(k)) + diff(B(i,k),q(j)) - diff(B(j,k),q(i)) );
            C(i,j) = C(i,j) + cijk*dq(k);
        end
    end
end

g = jacobian(U, q).';

equationOfMotion_option2 = B*ddq + C*dq + g == tau;

fprintf("Difference between Lagrange equations and B*ddq + C*dq + g (must be zero):\n")
disp(simplify(equationOfMotion - (B*ddq + C*dq + g)))
fprintf("--------------------------------------\n\n")


