function revolutePendulum(method)
% Test bilateral constraints for a spatial pendulum with a 1-dof revolute
% joint.  In the direction of rotation, everything is identical to the
% planar pendulum used as a test in the code testPendulum.m  Results
% matched exactly using Euler and Bender, but Midpoint caused energy loss.

% Written by Jeff Trinkle, 2012.

if ~exist('method', 'var')
  method = 1; % 0 ==> Euler,  1 ==> Bender,  2 ==> Midpoint
end

if method == 0
  disp('Using Euler method')
elseif method == 1
  disp('Using Bender method')
elseif method == 2
  disp('Using Midpoint method')
end

% Set up parameters for a slender rod of length L pinned at one end.  The
% rod starts horizontal.
L = 1;  
mass = 1; 
grav = 9.81;
r1_body = [-L/2; 0; 0]; % Location of joint in {B}
r2_body = r1_body + [0; 0; 1];
Izz = mass/12*L*L; 
Iyy = Izz * 10; 
Ixx = Izz * 3;
I_body = diag([Ixx  Iyy  Izz]); % Inertia matrix in {B}
M_body = zeros(6,6);
M_body(1:3, 1:3) = mass * eye(3,3);
M_body(4:6, 4:6) = I_body;

% Initial conditions to match the planarPendulum
q_l = [L/2; 0; 0; 1; 0; 0; 0];  % Starting horizontal to right
nu_l = [0; 0; 0; 0; 0; 0]; % Zero initial velocity expressed in {W}
lam_l = [0; 0; 0; 0; 0]; % Constraint forces assumed to be zero initially

M_wrld_l = M_body;
wRb = ep2rot(q_l(4:7));
M_wrld(4:6, 4:6) = wRb * I_body * wRb';
KE_l = nu_l' * M_wrld_l * nu_l / 2;  
PE_l = mass * grav * (L/2 + q_l(2));  
TE_l = PE_l + KE_l;  

% Compute 5 constraint to make a revolute joint: 
%       (wTb * r1_body)(1:3) = 0
%       (wTb * r2_body)(1:2) = 0
[C_l, Cdot_l] = constraints(r1_body, r2_body, q_l, nu_l);
bJn_l = jacobian(r1_body, r2_body, q_l);

dt = 0.001;
T = 2;
t_all = 0 : dt : T;
numPts = length(t_all);
E_all = zeros(3,numPts);
E_all(:,1) = [KE_l; PE_l; TE_l];
nu_all = zeros(6,numPts);
nu_all(:,1) = nu_l;
q_all = zeros(7, numPts);
q_all(:,1) = q_l;
lam_all = zeros(5, numPts);
lam_all(:,1) = lam_l;
C_all = zeros(5, numPts);
C_all(:,1) = C_l;
Cdot_all = zeros(5, numPts);
Cdot_all(:,1) = Cdot_l;
i = 1; 

tic;
for t = dt : dt : T
    i = i + 1;
    [wRb_l, B_l] = ep2rot(q_l(4:7));
    I_wrld_l = wRb_l * I_body * wRb_l';
    M_wrld_l(4:6, 4:6) = I_wrld_l;
    
    bigVec = [M_wrld_l * nu_l + gapp(mass, I_body, q_l, nu_l, t-dt/2) * dt;
                 C_l / dt];
    bigMat = [  M_wrld_l     -bJn_l';
                -bJn_l    zeros(5,5)];
    soln = bigMat \ bigVec;
%     fprintf('sphere soln = %12.8f\n', soln);

    nu_lp1 = soln(1:6);
    lam_lp1 = soln(7:11) / dt;
    V_l = [ eye(3,3)   zeros(3,3);
         zeros(4,3)    B_l * wRb_l' ];
    q_lp1 = q_l + V_l * nu_lp1 * dt;
    q_lp1(4:7) = q_lp1(4:7) / norm(q_lp1(4:7));
%     fprintf('planar q_lp1 = %12.8f\n', q_lp1);

    % Constraint errors at end of time step
    [C_lp1, Cdot_lp1] = constraints(r1_body, r2_body, q_lp1, nu_lp1);
%     fprintf('sphere C_lp1 = %12.8f\n', C_lp1);
%     fprintf('sphere Cdot_lp1 = %12.8f\n', Cdot_lp1);
    
    if method == 0
        % Do nothing.  Using Euler step above.
    elseif method == 1
        % Apply Bender's impulse-based error correction on positions first.
        MinvJtran_l = M_wrld_l \ bJn_l';
        iters = 0;
        while norm(C_lp1) > 1e-6  &&  iters < 4
            iters = iters + 1;
            % Compute impulse to apply at start of time step (therefore it
            % affects configurations, not just velocities).
            del_p = -inv(bJn_l * MinvJtran_l) * C_lp1 / dt;
            %         fprintf('sphere del_p = %12.8f\n', del_p);
            lam_lp1 = lam_lp1 + del_p / dt;
            del_nu = MinvJtran_l * del_p;
            nu_lp1 = nu_lp1 + del_nu;
            q_lp1 = q_lp1 + V_l * del_nu * dt;
            q_lp1(4:7) = q_lp1(4:7) / norm(q_lp1(4:7));
            
            % Compute constraint errors again
            [C_lp1, Cdot_lp1] = constraints(r1_body, r2_body, q_lp1, nu_lp1);
        end
        
        % Next, Bender's correction is applied at the end of the time step
        % to correct the velocity error - like an instantaneous impulse at
        % the end of the time step.
        iters = 0;
        bJn_lp1 = jacobian(r1_body, r2_body, q_lp1);
        wRb_lp1 = ep2rot(q_lp1(4:7));
        I_wrld_lp1 = wRb_lp1 * I_body * wRb_lp1';
        M_wrld_lp1 = M_body;
        M_wrld_lp1(4:6, 4:6) = I_wrld_lp1;
        MinvJtran_lp1 = M_wrld_lp1 \ bJn_lp1';
        % Velocity correction loop
        while norm(Cdot_lp1) > 1e-6  &&  iters < 4
            iters = iters + 1;
            del_p = -inv(bJn_lp1 * MinvJtran_lp1) * Cdot_lp1;
            %         fprintf('sphere del_p = %12.8f\n', del_p);
            lam_lp1 = lam_lp1 + del_p / dt;
            del_nu = MinvJtran_lp1 * del_p;
            nu_lp1 = nu_lp1 + del_nu;
            
            % Compute constraint errors again
            [C_lp1, Cdot_lp1] = constraints(r1_body, r2_body, q_lp1, nu_lp1);            
        end
    elseif method == 2
        nu = (nu_l + nu_lp1) / 2; % average velocity over time step
        q = q_l + V * nu/2 * dt; % update q at midpoint
        q(4:7) = q(4:7) / norm(q(4:7));
        [wRb, B] = ep2rot(q(4:7)); % R and B at mid point
        I_wrld = wRb * I_body * wRb'; % inertia matrix at mid point
        M_wrld(4:6, 4:6) = I_wrld;
        bJn = jacobian(r1_body, r2_body, q);

        bigVec = [M_wrld*nu_l + gapp(mass, I_body, q, nu, t-dt/2)*dt;
                   C_l / dt];  % at start of time step
        bigMat = [  M_wrld     -bJn';
                  -bJn       zeros(5,5)];  % at mid point
        soln = bigMat \ bigVec;
        
        nu_lp1 = soln(1:6);
        lam_lp1 = soln(7:11) / dt;
        V = [ eye(3,3)   zeros(3,3);
             zeros(4,3)    B * wRb' ]; % at mid point
        q_lp1 = q_l + V * nu_lp1 * dt;
        q_lp1(4:7) = q_lp1(4:7) / norm(q_lp1(4:7));
        
        % Constraint errors at end of time step
        [C_lp1, Cdot_lp1] = constraints(r1_body, r2_body, q_lp1, nu_lp1);
    end
    
    wRb_lp1 = ep2rot(q_lp1(4:7));
    M_wrld_lp1(1:3, 1:3) = M_body(1:3, 1:3);
    M_wrld_lp1(4:6, 4:6) = wRb_lp1 * I_body * wRb_lp1';
    KE_lp1 = nu_lp1' * M_wrld_lp1 * nu_lp1 / 2;
    PE_lp1 = mass * grav * (L/2 + q_lp1(2));
    TE_lp1 = PE_lp1 + KE_lp1;

    E_all(:,i) = [KE_lp1; PE_lp1; TE_lp1];
    nu_all(:,i) = nu_lp1;
    q_all(:,i) = q_lp1;
    lam_all(:,i) = lam_lp1;
    C_all(:,i) = C_lp1;
    Cdot_all(:,i) = Cdot_lp1;
    
    nu_l = nu_lp1;
    q_l = q_lp1;
    C_l = C_lp1;
    % Need Jacobian for next time step too.
    bJn_l = jacobian(r1_body, r2_body, q_l);
end
toc
figure
plot(t_all, q_all(1:3, :));
legend('x', 'y', 'z');
grid;
figure
plot(t_all, q_all(4:7, :));
legend('ep0', 'ep1', 'ep2', 'ep3');
grid;
figure
plot(t_all, nu_all);
legend('v_x', 'v_y', 'v_z', '\omega_x', '\omega_y', '\omega_z');
grid;
figure
plot(t_all, lam_all);
legend('\lambda_x1', '\lambda_y1', '\lambda_z1', '\lambda_x2', '\lambda_y2');
grid;
figure
plot(t_all, C_all);
legend('C_x1', 'C_y1', 'C_z1', 'C_x2', 'C_y2');
grid;
figure
plot(t_all, Cdot_all);
legend('Cdot_x1', 'Cdot_y1', 'Cdot_z1', 'Cdot_x2', 'Cdot_y2');
grid;
figure
plot(t_all, E_all);
legend('Kinetic', 'Potential', 'Total');
grid;
end

function [C, Cdot] = constraints(r1_body, r2_body, q, nu)
quat = q(4:7);
wRb = ep2rot(quat);
wTb = [   wRb     q(1:3);
       zeros(1,3)    1     ];
% Put a point on the body at the origin of the world frame and a second
% point somewhere on the z-axis away from the origin.
r1_wrld_h = wTb * [r1_body; 1];
r2_wrld_h = wTb * [r2_body; 1];
C = [r1_wrld_h(1:3);  r2_wrld_h(1:2)];
bJn = jacobian(r1_body, r2_body, q);
Cdot = bJn * nu;
end

function bJn = jacobian(r1_body, r2_body, q)
quat = q(4:7);
wRb = ep2rot(quat);
bJn1 = [eye(3,3)   -cross_matrix(wRb * r1_body)];
bJn2 = [eye(3,3)   -cross_matrix(wRb * r2_body)];
bJn = [bJn1; bJn2(1:2, :)];
end

function wrench = gapp(mass, I_body, q, nu, t)
quat = q(4:7);
wRb = ep2rot(quat);
omega_w = nu(4:6);
I_wrld = wRb * I_body * wRb';
omegaCrsIomega_w = cross(omega_w, I_wrld * omega_w);
wrench = [ 0;    -9.81*mass;    0;   zeros(3,1)];
end

function wx = cross_matrix(w)
wx = [ 0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end
