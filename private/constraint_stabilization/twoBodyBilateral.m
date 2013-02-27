function twoBodyBilateral
% This code was tested against testFreeRigidBody.m with zero constraints
% and found to agree exactly.

% Written by Jeff Trinkle, 2012.
nconstr = 5;

% Test bilateral constraints for a two-body system
method = 1; % 0 ==> Euler,  1 ==> Bender

% Mass quantities for bodies A and B
lA = 2;
massA = 1;
% Izz = massA/12*lA*lA; 
% Iyy = Izz * 1; 
% Ixx = Izz * 1;
Izz = 1; 
Iyy = 10; 
Ixx = 3;
I_bodyA = diag([Ixx  Iyy  Izz]); % Inertia matrix of body A in {A}
Ma = zeros(6,6);
Ma(1:3, 1:3) = massA * eye(3,3);

lB = 2;
massB = 1;
% Izz = massB/12*lB*lB; 
% Iyy = Izz * 6; 
% Ixx = Izz * 2;
Izz = 1; 
Iyy = 6; 
Ixx = 2;
I_bodyB = diag([Ixx  Iyy  Izz]); % Inertia matrix of body B in {B}
Mb = zeros(6,6);
Mb(1:3, 1:3) = massB * eye(3,3);

M_body = zeros(12,12);
M_body(1:6, 1:6) = Ma;
M_body(7:12, 7:12) = Mb;

% Data for constraints
AposC = [1; 0; 0]; % Position of constraint frame {C} origin wrt {A}. 
AquatC = [1; 0; 0; 0]; % Unit quaternion form of orientation of {C} wrt {A}
aRc = ep2rot(AquatC); % Rotation matrix from {A} to {C}
aTc = [ aRc   AposC;
       0 0 0    1   ]; % hTform from {A} to {C}
   
BposC = [-1; 0; 0]; % Position of constraint frame {C} origin wrt {B}
BquatC = [1; 0; 0; 0]; % Unit quaternion form of orientation of {C} wrt {B}
bRc = ep2rot(BquatC); % Rotation matrix from {B} to {C}
bTc = [ bRc   BposC;
       0 0 0    1   ]; % hTform from {B} to {C}

% These points are the origin, z-axis, and x-axis of {C} expressed in {C}.
CconstrPts = [0 0 1;   % Don't know how aTc and bTc will be specified.
              0 0 0;   % so get the constraint points from {C} and map
              0 1 0;   % them into {A} an {B}
              1 1 1]; % Make columns homogeneous points.
aPts = aTc * CconstrPts; % Express this points in {A}.  These are constant. 
bPts = bTc * CconstrPts; % Express them in {B}.  These are constant.

% Initial conditions
qA = [-1; 2.000000; 0; 1; 0; 0; 0];  % Starting horizontal to left
nuA = [0; 0; 0; 0; 0; 0]; % Initial velocity twist expressed in {W}
% qB = [1; 2; 0; 1; 0; 0; 0];  % Starting horizontal to right
qB = [0; 1; 0; sqrt(2)/2; 0; 0; -sqrt(2)/2];  % Starting hanging straight down
nuB = [0; 0; 0; 0; 1; 0]; % Initial velocity twist expressed in {W}


q_l = [qA; qB];
nu_l = [nuA; nuB];
M_wrld = M_body;
wRa = ep2rot(q_l(4:7));
M_wrld(4:6, 4:6) = wRa * I_bodyA * wRa';
wRb = ep2rot(q_l(11:14));
M_wrld(10:12, 10:12) = wRb * I_bodyB * wRb';
ke_l = [  nu_l(1:6)' * M_wrld(1:6,1:6)   * nu_l(1:6) / 2;
         nu_l(7:12)' * M_wrld(7:12,7:12) * nu_l(7:12) / 2];
if nconstr > 0
    lam_l = zeros(nconstr, 1);
    % Compute constraints
    [C_l, Cdot_l] = constraints(q_l, nu_l, aTc, bTc, aPts, bPts);
    C_l = C_l(1:nconstr);
    Cdot_l = Cdot_l(1:nconstr);
    Jn_l = jacobian(q_l, aTc, bTc, aPts, bPts);
    Jn_l = Jn_l(1:nconstr, :);
end

dt = 0.01;
T = 10;
t_all = 0 : dt : T;
numPts = length(t_all);
ke_all = zeros(2, numPts);
ke_all(1:2,1) = ke_l;
nu_all = zeros(12, numPts);
nu_all(:,1) = nu_l;
q_all = zeros(14, numPts);
q_all(:,1) = q_l;
if nconstr > 0
    lam_all = zeros(nconstr, numPts);
    lam_all(1:nconstr,1) = lam_l(1:nconstr);
    C_all = zeros(nconstr, numPts);
    C_all(1:nconstr,1) = C_l(1:nconstr);
    Cdot_all = zeros(nconstr, numPts);
    Cdot_all(1:nconstr,1) = Cdot_l(1:nconstr);
end
i = 1; 

tic;
for t = dt : dt : T
    i = i + 1;
    [wRa, Ba] = ep2rot(q_l(4:7));
    [wRb, Bb] = ep2rot(q_l(11:14));
    I_wrldA = wRa * I_bodyA * wRa';
    I_wrldB = wRb * I_bodyB * wRb';
    M_wrld(4:6,   4:6) =   I_wrldA;
    M_wrld(10:12, 10:12) = I_wrldB;
    
    if nconstr == 0
        bigMat = M_wrld;
        bigVec = M_wrld*nu_l + gapp(massA, massB, I_bodyA, I_bodyB, q_l, nu_l, t-dt/2) * dt;
    else
        bigMat = [M_wrld                -Jn_l(1:nconstr,:)';
                  -Jn_l(1:nconstr,:)   zeros(nconstr, nconstr)];
        bigVec = [M_wrld*nu_l + gapp(massA, massB, I_bodyA, I_bodyB, q_l, nu_l, t-dt/2) * dt;
                 1 * C_l(1:nconstr) / dt - 0 * Cdot_l(1:nconstr)];
    end
    soln = bigMat \ bigVec;
%     fprintf('sphere soln = %12.8f\n', soln);

    nu_lp1 = soln(1:12);
    if nconstr > 0
        lam_lp1 = soln(13:13+nconstr-1) / dt;
    end
    Va = [ eye(3,3)    zeros(3,3);
          zeros(4,3)    Ba * wRa' ];
    Vb = [ eye(3,3)    zeros(3,3);
          zeros(4,3)    Bb * wRb' ];
    V = zeros(14, 12);
    V(1:7, 1:6) = Va;
    V(8:14, 7:12) = Vb;
    
    q_lp1 = q_l + V * nu_lp1 * dt;
    q_lp1(4:7) = q_lp1(4:7) / norm(q_lp1(4:7));
    q_lp1(11:14) = q_lp1(11:14) / norm(q_lp1(11:14));
%     fprintf('planar q_lp1 = %12.8f\n', q_lp1);

    if nconstr > 0
        % Constraint errors at end of time step
        [C_lp1, Cdot_lp1] = constraints(q_lp1, nu_lp1, aTc, bTc, aPts, bPts);
        C_lp1 = C_lp1(1:nconstr);
        Cdot_lp1 = Cdot_lp1(1:nconstr);
        % fprintf('sphere C_lp1 = %12.8f\n', C_lp1);
        % fprintf('sphere Cdot_lp1 = %12.8f\n', Cdot_lp1);
    end
    
    if method == 0
        % Do nothing.  Using Euler step above.
    elseif method == 1
        % Apply Bender's impulse-based error correction on positions first.
        MinvJtran = M_wrld \ Jn_l';
        iters = 0;
        while norm(C_lp1) > 1e-10  &&  iters < 4
            iters = iters + 1;
            % Compute impulse to apply at start of time step (therefore it
            % affects configurations, not just velocities).
            del_p = -inv(Jn_l * MinvJtran) * C_lp1 / dt;
            %         fprintf('sphere del_p = %12.8f\n', del_p);
            lam_lp1 = lam_lp1 + del_p / dt;
            del_nu = MinvJtran * del_p;
            nu_lp1 = nu_lp1 + del_nu;
            q_lp1 = q_lp1 + V * del_nu * dt;
            q_lp1(4:7) = q_lp1(4:7) / norm(q_lp1(4:7));
            q_lp1(11:14) = q_lp1(11:14) / norm(q_lp1(11:14));
            
            % Compute constraint errors again
            [C_lp1, Cdot_lp1] = constraints(q_lp1, nu_lp1, aTc, bTc, aPts, bPts);
            C_lp1 = C_lp1(1:nconstr);
            Cdot_lp1 = Cdot_lp1(1:nconstr);
        end
        
        % Next, Bender's correction is applied at the end of the time step
        % to correct the velocity error - like an instantaneous impulse at
        % the end of the time step.
        iters = 0;
        Jn_lp1 = jacobian(q_lp1, aTc, bTc, aPts, bPts);
        Jn_lp1 = Jn_lp1(1:nconstr, :);
        MinvJtran = M_wrld \ Jn_lp1';
        % Velocity correction loop
        while norm(Cdot_lp1) > 1e-6  &&  iters < 4
            iters = iters + 1;
            del_p = -inv(Jn_lp1 * MinvJtran) * Cdot_lp1;
            %         fprintf('sphere del_p = %12.8f\n', del_p);
            lam_lp1 = lam_lp1 + del_p / dt;
            del_nu = MinvJtran * del_p;
            nu_lp1 = nu_lp1 + del_nu;
            
            % Compute constraint errors again
            [C_lp1, Cdot_lp1] = constraints(q_lp1, nu_lp1, aTc, bTc, aPts, bPts);            
            C_lp1 = C_lp1(1:nconstr);
            Cdot_lp1 = Cdot_lp1(1:nconstr);
        end
    end
    
    nu_all(:,i) = nu_lp1;
    q_all(:,i) = q_lp1;
    wRa = ep2rot(q_l(4:7));
    M_wrld(4:6, 4:6) = wRa * I_bodyA * wRa';
    wRb = ep2rot(q_l(11:14));
    M_wrld(10:12, 10:12) = wRb * I_bodyB * wRb';
    ke_all(:, i) = [  nu_l(1:6)' * M_wrld(1:6,1:6)   * nu_l(1:6) / 2;
                     nu_l(7:12)' * M_wrld(7:12,7:12) * nu_l(7:12) / 2];

     if nconstr > 0
        lam_all(:,i) = lam_lp1(1:nconstr);
        C_all(:,i) = C_lp1(1:nconstr);
        Cdot_all(:,i) = Cdot_lp1(1:nconstr);
        C_l = C_lp1;
        % Need Jacobian for next time step too.
        Jn_l = jacobian(q_l, aTc, bTc, aPts, bPts);
        Jn_l = Jn_l(1:nconstr, :);
    end
    nu_l = nu_lp1;
    q_l = q_lp1;
end
toc
figure
plot(t_all, ke_all);
legend('KE_A', 'KE_B');
grid;
figure
plot(t_all, q_all(1:3, :));
legend('x_A', 'y_A', 'z_A');
grid;
figure
plot(t_all, q_all(8:10, :));
legend('x_B', 'y_B', 'z_B');
grid;
figure
plot(t_all, q_all(4:7, :));
legend('ep_{A0}', 'ep_{A1}', 'ep_{A2}', 'ep_{A3}');
grid;
figure
plot(t_all, q_all(11:14, :));
legend('ep_{B0}', 'ep_{B1}', 'ep_{B2}', 'ep_{B3}');
grid;
figure
plot(t_all, nu_all(1:6, :));
legend('v_{Ax}', 'v_{Ay}', 'v_{Az}', '\omega_{Ax}', '\omega_{Ay}', '\omega_{Az}');
grid;
figure
plot(t_all, nu_all(7:12, :));
legend('v_{Bx}', 'v_{By}', 'v_{Bz}', '\omega_{Bx}', '\omega_{By}', '\omega_{Bz}');
grid;
if nconstr > 0
    figure
    plot(t_all, lam_all);
    legend('\lambda_{1}', '\lambda_{2}', '\lambda_{3}', '\lambda_{4}', '\lambda_{5}');
    grid;
    figure
    plot(t_all, C_all(1:nconstr, :));
    legend('C_{1}', 'C_{2}', 'C_{3}', 'C_{4}', 'C_{5}');
    grid;
    figure
    plot(t_all, Cdot_all(1:nconstr, :));
    legend('Cdot_{1}', 'Cdot_{2}', 'Cdot_{3}', 'Cdot_{4}', 'Cdot_{5}');
    grid;
end
end

function [C, Cdot] = constraints(q, nu, aTc, bTc, aPts, bPts)
% Compute constraint in {C}
quatA = q(4:7);         
quatB = q(11:14);
wRa = ep2rot(quatA);     
wRb = ep2rot(quatB);
wTa = [   wRa     q(1:3);
       zeros(1,3)    1     ];
wTb = [   wRb     q(8:10);
       zeros(1,3)    1     ];
wPtsA = wTa * aPts;
wPtsB = wTb * bPts;
Wdiff = wPtsA - wPtsB;

% positions of points in {A} minus those in {B}.  Error points to {A}
% Cdiff = inv(aTc) * aPts - inv(aTc * wTa) * wTb * bPts;

% C = [Cdiff(1:2,1); Cdiff(1:2,2); Cdiff(2,3)]; % Prismatic
% C = [Cdiff(1:3,1); Cdiff(1:2,1)]; % Revolute
C = [Wdiff(1:3,1); Wdiff(1:2,1)]; % Revolute

aRc = aTc(1:3, 1:3);
bRc = bTc(1:3, 1:3);
Jn = jacobian(q, aTc, bTc, aPts, bPts);
Cdot = Jn * nu;
end

function Jn = jacobian(q, aTc, bTc, aPts, bPts)
quatA = q(4:7);                quatB = q(11:14);
wRa = ep2rot(quatA);           wRb = ep2rot(quatB);
cTa = inv(aTc);                cTb = inv(bTc);
cRa = cTa(1:3, 1:3);           cRb = cTb(1:3, 1:3);
% aJn1 = (cRa * wRa') * [eye(3,3)   -cross_matrix(aPts(1:3, 1))];
% aJn2 = (cRa * wRa') * [eye(3,3)   -cross_matrix(aPts(1:3, 2))];
% aJn3 = (cRa * wRa') * [eye(3,3)   -cross_matrix(aPts(1:3, 3))];
% bJn1 = (cRb * wRb') * [eye(3,3)   -cross_matrix(bPts(1:3, 1))];
% bJn2 = (cRb * wRb') * [eye(3,3)   -cross_matrix(bPts(1:3, 2))];
% bJn3 = (cRb * wRb') * [eye(3,3)   -cross_matrix(bPts(1:3, 3))];
W_aPts = wRa * aPts(1:3,:);            W_bPts = wRb * bPts(1:3,:);          
aJn1 = [eye(3,3)   -cross_matrix(W_aPts(1:3, 1))];
aJn2 = [eye(3,3)   -cross_matrix(W_aPts(1:3, 2))];
aJn3 = [eye(3,3)   -cross_matrix(W_aPts(1:3, 3))];
bJn1 = [eye(3,3)   -cross_matrix(W_bPts(1:3, 1))];
bJn2 = [eye(3,3)   -cross_matrix(W_bPts(1:3, 2))];
bJn3 = [eye(3,3)   -cross_matrix(W_bPts(1:3, 3))];
% % Prismatic
% Jn = [-aJn1(1:2,:)   bJn1(1:2,:);
%       -aJn2(1:2,:)   bJn2(1:2,:);
%       -aJn3(2,:)     bJn3(2,:)   ];
% Revolute
Jn = [aJn1(1:3,:)   -bJn1(1:3,:);
      aJn2(1:2,:)   -bJn2(1:2,:)];
end

function wrench = gapp(massA, massB, I_bodyA, I_bodyB, q, nu, t)
quatA = q(4:7);
wRa = ep2rot(quatA);
omegaA = nu(4:6);  % expressed in {W}
I_wrldA = wRa * I_bodyA * wRa';
omegaAcrossIomegaA = cross(omegaA, I_wrldA * omegaA);

quatB = q(11:14);
wRb = ep2rot(quatB);
omegaB = nu(10:12);  % expressed in {W}
I_wrldB = wRb * I_bodyB * wRb';
omegaBcrossIomegaB = cross(omegaB, I_wrldB * omegaB);

kd = 0;
% wrenchA = [ 0;     0.1*massA;    1;
%           -kd*nu(4:6) - omegaAcrossIomegaA + zeros(3,1)];
% wrenchB = [ 0;    -0.1*massB;    0;
%           -kd*nu(10:12) - omegaBcrossIomegaB + sin(5*t)*[1; 1; 1]];
wrenchA = [ 0;    1;    0;  -0.2 * nu(4:6) - omegaAcrossIomegaA + sin(5*t)*ones(3,1)];
wrenchB = [ 0;    -1;    0;  -0.2 * nu(10:12)  - omegaBcrossIomegaB];
wrench = [wrenchA; wrenchB];
% kd = 2;
% wrench = [ 0;    -9.81*mass;    sin(3*t);
%           -kd*nu(4:6) - omegaCrsIomega_w + 3*sin(5*t)*[1; 1; 1]];
end

function wx = cross_matrix(w)
wx = [ 0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end
