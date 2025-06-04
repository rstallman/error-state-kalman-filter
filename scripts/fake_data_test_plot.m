clear; close all; format compact;
addpath('filters');
addpath('helper');
addpath('thirdparty/shadedErrorBar');
% -------------------------------------------------------------------------
% Fake data time limits and resolution lower resolution with noise free
% data should cause the prediction to improve.
time.tmin = 0;
time.tmax = 2;
time.dt = 1e-3;

% Covariance for sensor noise set noise.add_noise = false 
% to have perfect noise free data.
% If add_noise if false, accel_noise and gyro_noise can be anything.
%
% Q1 and Q2 are two random matrices to try to create really bad noise
% with a ton of cross corelation between states
rng(2)
noise.add_noise = true;
m = 100;
Q1 = randn(3,3);
Q2 = randn(3,3);
noise.accel_noise = (6*9.81)^2*eye(3);
noise.gyro_noise = eye(3);

% Set the frequency of the correction step (Hz)
%  - Increase the frequency to test robustness of filter
%    to slow updates
f_cor = 1;
dt_cor = 1/f_cor;

% Generate the fake data see function definition for input and output
% variable definitions.
[omega, accel, ~, ~, gt, init, wf_data] = gen_fake_data(time, noise);

% Set the time data from the accelerometer
t = accel.t;
N = length(t);


% -------------------------------------------------------------------------
% Helper functions, mostly for SO3 stuff
function ux  = skew(u)
    ux = [
        0 -u(3) u(2)
        u(3) 0 -u(1)
        -u(2) u(1) 0
    ];
end

function u = unskew(ux)
    u(1,1) = -ux(2,3);
    u(2,1) = ux(1,3);
    u(3,1) = -ux(1,2);
end

function w = Log(R)
    w = unskew(logm(R));
end

function J = J_l(theta)
    t_x = skew(theta);
    t = norm(theta);

    J = eye(3) + (1 - cos(t))/t^2 * t_x + (t - sin(t))/t^3*(t_x)^2;
end
%-------------------------------------------------------------

%---------------------------------------------------------------
function eul = rotm2eul(rotm, sequence)
    if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('rotm2eul: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    eul = zeros(3,1);

    %% Compute the Euler angles theta for the x, y and z-axis from a rotation matrix R, in
    %  dependency of the specified axis rotation sequence for the rotation factorization:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, pp. 9-10, 16-17.
    %   [2] Computing Euler angles from a rotation matrix, Gregory G. Slabaugh, <http://www.staff.city.ac.uk/~sbbh653/publications/euler.pdf>
    %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 30-33, formulas (2.19), (2.19'), (2.21) and (2.21').
    switch sequence
        case 'ZYX'
            % convention used by (*) and (**).
            % note: the final orientation is the same as in XYZ order about fixed axes ...
            if (rotm(3,1) < 1)
                if (rotm(3,1) > -1) % case 1: if r31 ~= ±1
                    % Solution with positive sign. It limits the range of the values
                    % of theta_y to (-pi/2, pi/2):
                    eul(1,1) = atan2(rotm(2,1), rotm(1,1)); % theta_z
                    eul(2,1) = asin(-rotm(3,1));            % theta_y
                    eul(3,1) = atan2(rotm(3,2), rotm(3,3)); % theta_x
                else % case 2: if r31 = -1
                    % theta_x and theta_z are linked --> Gimbal lock:
                    % There are infinity number of solutions for theta_x - theta_z = atan2(-r23, r22).
                    % To find a solution, set theta_x = 0 by convention.
                    eul(1,1) = -atan2(-rotm(2,3), rotm(2,2));
                    eul(2,1) = pi/2;
                    eul(3,1) = 0;
                end
            else % case 3: if r31 = 1
                % Gimbal lock: There is not a unique solution for
                %   theta_x + theta_z = atan2(-r23, r22), by convention, set theta_x = 0.
                eul(1,1) = atan2(-rotm(2,3), rotm(2,2));
                eul(2,1) = -pi/2;
                eul(3,1) = 0;
            end
        case 'ZYZ'
            % convention used by (*)
            if (rotm(3,3) < 1)
                if (rotm(3,3) > -1)
                    % Solution with positive sign, i.e. theta_y is in the range (0, pi):
                    eul(1,1) = atan2(rotm(2,3),  rotm(1,3)); % theta_z1
                    eul(2,1) = acos(rotm(3,3));              % theta_y (is equivalent to atan2(sqrt(r13^2 + r23^2), r33) )
                    eul(3,1) = atan2(rotm(3,2), -rotm(3,1)); % theta_z2
                else % if r33 = -1:
                    % Gimbal lock: infinity number of solutions for
                    %   theta_z2 - theta_z1 = atan2(r21, r22), --> set theta_z2 = 0.
                    eul(1,1) = -atan2(rotm(2,1), rotm(2,2)); % theta_z1
                    eul(2,1) = pi;                           % theta_y
                    eul(3,1) = 0;                            % theta_z2
                end
            else % if r33 = 1:
                % Gimbal lock: infinity number of solutions for
                %    theta_z2 + theta_z1 = atan2(r21, r22), --> set theta_z2 = 0.
                eul(1,1) = atan2(rotm(2,1), rotm(2,2)); % theta_z1
                eul(2,1) = 0;                           % theta_y
                eul(3,1) = 0;                           % theta_z2
            end
        % case 'ZYZ-'
        %     % convention used by (**)
        %     if (rotm(3,3) < 1)
        %         if (rotm(3,3) > -1)
        %             % Variant with negative sign. This is a derived solution
        %             % which produces the same effects as the solution above.
        %             % It limits the values of theta_y in the range of (-pi,0):
        %             eul(1,1) = atan2(-rotm(2,3), -rotm(1,3)); % theta_z1
        %             eul(2,1) = -acos(rotm(3,3));              % theta_y (is equivalent to atan2(-sqrt(r13^2 + r23^2), r33) )
        %             eul(3,1) = atan2(-rotm(3,2),  rotm(3,1)); % theta_z2
        %         else % if r33 = -1:
        %             % Gimbal lock: infinity number of solutions for
        %             %   theta_z2 - theta_z1 = atan2(-r12, -r11), --> set theta_z2 = 0.
        %             eul(1,1) = -atan2(-rotm(1,2), -rotm(1,1)); % theta_z1  (correct ???)
        %             eul(2,1) = -pi;                            % theta_y
        %             eul(3,1) = 0;                              % theta_z2
        %         end
        %     else % if r33 = 1:
        %         % Gimbal lock: infinity number of solutions for
        %         %    theta_z2 + theta_z1 = atan2(-r12, -r11), --> set theta_z2 = 0.
        %         eul(1,1) = atan2(-rotm(1,2), -rotm(1,1)); % theta_z1  (correct ???)
        %         eul(2,1) = 0;                             % theta_y
        %         eul(3,1) = 0;                             % theta_z2
        %     end
        otherwise
            error('rotm2eul: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end

% (*)  ... The Geometric Tools Engine (http://www.geometrictools.com),
% (**) ... The Robotics System Toolbox for Matlab (http://mathworks.com/help/robotics/index.html).

function rotm = eul2rotm(eul, sequence)
    if (size(eul,1) ~= 3)
        error('eul2rotm: %s', 'WBM.wbmErrorMsg.WRONG_VEC_DIM');
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    rotm = zeros(3,3);

    s_1 = sin(eul(1,1)); % theta_z or theta_z1
    c_1 = cos(eul(1,1));
    s_2 = sin(eul(2,1)); % theta_y
    c_2 = cos(eul(2,1));
    s_3 = sin(eul(3,1)); % theta_x or theta_z2
    c_3 = cos(eul(3,1));

    %% Convert the given Euler angles theta for the x, y and z-axis into the corresponding
    %  direction cosine rotation matrix R, in dependency of the axis rotation sequence for
    %  the multiplication order of the rotation factors:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, p. 9 & 16.
    %   [2] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>, p. 4.
    %   [3] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [4] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 31-32, formulas (2.18) and (2.20).
    switch sequence
        case 'ZYX'
            %            |c_1*c_2    c_1*s_2*s_3 - s_1*c_3    c_1*s_2*c_3 + s_1*s_3|
            % R(Theta) = |s_1*c_2    s_1*s_2*s_3 + c_1*c_3    s_1*s_2*c_3 - c_1*s_3|
            %            |   -s_2                  c_2*s_3                  c_2*c_3|
            rotm(1,1) =  c_1*c_2;
            rotm(1,2) =  c_1*s_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2*c_3 + s_1*s_3;

            rotm(2,1) =  s_1*c_2;
            rotm(2,2) =  s_1*s_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2*c_3 - c_1*s_3;

            rotm(3,1) = -s_2;
            rotm(3,2) =  c_2*s_3;
            rotm(3,3) =  c_2*c_3;
        case 'ZYZ'
            %            |c_1*c_2*c_3 - s_1*s_3   -c_1*c_2*s_3 - s_1*c_3    c_1*s_2|
            % R(Theta) = |s_1*c_2*c_3 + c_1*s_3   -s_1*c_2*s_3 + c_1*c_3    s_1*s_2|
            %            |             -s_2*c_3                  s_2*s_3        c_2|
            rotm(1,1) =  c_1*c_2*c_3 - s_1*s_3;
            rotm(1,2) = -c_1*c_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2;

            rotm(2,1) =  s_1*c_2*c_3 + c_1*s_3;
            rotm(2,2) = -s_1*c_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2;

            rotm(3,1) = -s_2*c_3;
            rotm(3,2) =  s_2*s_3;
            rotm(3,3) =  c_2;
        otherwise
            error('eul2rotm: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end
% -------------------------------------------------------------------------
% Initialize the solution / gt vectors
p_gt = [gt.x;gt.y;gt.z];
theta_gt = zeros(3, N);
theta_gt(:,1) = Log(gt.R{1});

p_ekf = zeros(3,N);
p_ekf_var = zeros(3,N);
theta_ekf = zeros(3, N);
theta_ekf_var = zeros(3,N);
p_liekf = zeros(3,N);
p_liekf_var = zeros(3,N);
theta_liekf = zeros(3, N);
theta_liekf_var = zeros(3,N);
% -------------------------------------------------------------------------
% Initialize the filter (with initial condition)
% Note: the polynomial function created by gen_fake_data almost definitely
% wont be zero at t = 0
ekf = EKF(rotm2eul(init.R0), init.p0, init.v0, eye(3) * 100);

p_ekf(:,1) = ekf.mu(1:3);
p_ekf_var(:,1) = sqrt(diag(ekf.Sigma(1:3,1:3)));
theta_ekf(:,1) = ekf.mu(7:9);
theta_ekf_var(:,1) = sqrt(diag(ekf.Sigma(7:9,7:9)));

% -------------------------------------------------------------------------
% Run the simulation on the data
t_cor = t(1);  %Time of first correction
for i = 1:N-1
    % Set dt off time data
    dt = t(i+1) - t(i);

    % Set the acceleration from the fake data
    a = [accel.x(i); accel.y(i); accel.z(i)];
    w = [omega.x(i); omega.y(i); omega.z(i)];
    
    % Run the ekf prediction step
    ekf.prediction(w, a, dt);
    
    % Run the ekf correction step
    if t(i) >= t_cor
        gps = [gt.x(i); gt.y(i); gt.z(i)];
        ekf.correction(gps)

        % Next correct at t + dt_cor
        t_cor = t(i) + dt_cor;
    end
    
    % Save the outputs (for plotting)
    variances = sqrt(diag(ekf.Sigma));
    p_ekf(:,i+1) = ekf.mu(1:3);
    theta_ekf(:,i+1) = Log(eul2rotm(ekf.mu(7:9)));
    p_ekf_var(:,i+1) = variances(1:3);
    theta_ekf_var(:,i+1) = variances(7:9);
    
    theta_gt(:,i+1) = Log(gt.R{i});
end

% -------------------------------------------------------------------------
% LIEKF
liekf = LIEKF(init.R0, init.p0, init.v0, ...
    eye(3)*.01, eye(3)*.01, eye(3)*.01, eye(3)*.01, eye(3)*.01);
lieTocartesian(liekf)
vars = sqrt(diag(liekf.sigma_cart));

p_liekf(:,1) = init.p0;
p_liekf_var(:,1) = vars(1:3);
theta_liekf(:,1) = Log(liekf.mu(1:3,1:3));
theta_liekf_var(:,1) = vars(1:3);
% Run the simulation on the data
t_cor = t(1);  %Time of first correction
for i = 1:N-1
    % Set dt off time data
    dt = t(i+1) - t(i);

    % Set the acceleration from the fake data
    a = [accel.x(i); accel.y(i); accel.z(i)];
    w = [omega.x(i); omega.y(i); omega.z(i)];
    
    % Run the ekf prediction step
    liekf.prediction(w, a, dt);
    
    % Run the ekf correction step
    if t(i) >= t_cor
        gps = [gt.x(i); gt.y(i); gt.z(i)];
        liekf.correction(gps)

        % Next correct at t + dt_cor
        t_cor = t(i) + dt_cor;
    end

    % Extract the state from the filter
    [R, p, v] = liekf.getState(); 
    lieTocartesian(liekf)
    
    % Save the outputs (for plotting)
    p_liekf(:,i+1) = p;
    theta_liekf(:,i+1) = Log(R);
    
    vars = sqrt(diag(liekf.sigma_cart));
    p_liekf_var(:,i+1) = vars(7:9);
    theta_liekf_var(:,i+1) = vars(1:3);
end
% -------------------------------------------------------------------------
% Plot position and theta data to visualize
% the operation of the filter
figure;
subplot(311)
hold('on')
plot(t, p_gt(1,:), 'k--', 'LineWidth', 2);
plot(t, p_ekf(1,:), 'g', 'LineWidth', 1);
plot(t, p_liekf(1,:), 'r', 'LineWidth', 1);
axis([0,2,-200,200])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
title("Position");
subplot(312)
hold('on')
plot(t, p_gt(2,:),  'k--', 'LineWidth', 2)
plot(t, p_ekf(2,:), 'g', 'LineWidth', 1);
plot(t, p_liekf(2,:), 'r', 'LineWidth', 1)
axis([0,2,-400,0])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
subplot(313)
hold('on')
plot(t, p_gt(3,:), 'k--', 'LineWidth', 2)
plot(t, p_ekf(3,:), 'g', 'LineWidth', 1);
plot(t, p_liekf(3,:), 'r', 'LineWidth', 1)
axis([0,2,-300,100])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
print('position_noise', '-dpng')

figure;
subplot(311)
hold('on')
plot(t, theta_gt(1,:), 'k--', 'LineWidth', 2);
plot(t, theta_ekf(1,:), 'g', 'LineWidth', 1);
plot(t, theta_liekf(1,:), 'r', 'LineWidth', 1);
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
title("Theta");
subplot(312)
hold('on')
plot(t, theta_gt(2,:), 'k--', 'LineWidth', 2)
plot(t, theta_ekf(2,:), 'g', 'LineWidth', 1);
plot(t, theta_liekf(2,:), 'r', 'LineWidth', 1)
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
subplot(313)
hold('on')
plot(t, theta_gt(3,:),  'k--', 'LineWidth', 2)
plot(t, theta_ekf(3,:), 'g', 'LineWidth', 1);
plot(t, theta_liekf(3,:), 'r', 'LineWidth', 1)
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
print('theta_noise', '-dpng')

% -------------------------------------------------------------------------
% Plot position and theta data to visualize
% the operation of the filter
figure;
subplot(311)
hold('on')
plot(t, p_gt(1,:), 'k--', 'LineWidth', 2);
%plot(t, p_ekf(1,:), 'g', 'LineWidth', 1);
%plot(t, p_liekf(1,:), 'r', 'LineWidth', 1);
%plot(t, p_ekf(1,:)+3*p_ekf_var(1,:), 'b', 'LineWidth', 1);
%plot(t, p_ekf(1,:)-3*p_ekf_var(1,:), 'b', 'LineWidth', 1);
%plot(t, p_liekf(1,:)+3*p_liekf_var(1,:), 'm', 'LineWidth', 1);
%plot(t, p_liekf(1,:)-3*p_liekf_var(1,:), 'm', 'LineWidth', 1);
shadedErrorBar(t, p_ekf(1,:), 3*p_ekf_var(1,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, p_liekf(1,:), 3*p_liekf_var(1,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-200,200])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
title("Position");
subplot(312)
hold('on')
plot(t, p_gt(2,:),  'k--', 'LineWidth', 2)
%plot(t, p_ekf(2,:), 'g', 'LineWidth', 1);
%plot(t, p_liekf(2,:), 'r', 'LineWidth', 1)
%plot(t, p_ekf(2,:)+3*p_ekf_var(2,:), 'b', 'LineWidth', 1);
%plot(t, p_ekf(2,:)-3*p_ekf_var(2,:), 'b', 'LineWidth', 1);
%plot(t, p_liekf(2,:)+3*p_liekf_var(2,:), 'm', 'LineWidth', 1);
%plot(t, p_liekf(2,:)-3*p_liekf_var(2,:), 'm', 'LineWidth', 1);
shadedErrorBar(t, p_ekf(2,:), 3*p_ekf_var(2,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, p_liekf(2,:), 3*p_liekf_var(2,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-400,0])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
subplot(313)
hold('on')
plot(t, p_gt(3,:), 'k--', 'LineWidth', 2)
% plot(t, p_ekf(3,:), 'g', 'LineWidth', 1);
% plot(t, p_liekf(3,:), 'r', 'LineWidth', 1)
% plot(t, p_ekf(3,:)+3*p_ekf_var(3,:), 'b', 'LineWidth', 1);
% plot(t, p_ekf(3,:)-3*p_ekf_var(3,:), 'b', 'LineWidth', 1);
% plot(t, p_liekf(3,:)+3*p_liekf_var(3,:), 'm', 'LineWidth', 1);
% plot(t, p_liekf(3,:)-3*p_liekf_var(3,:), 'm', 'LineWidth', 1);
shadedErrorBar(t, p_ekf(3,:), 3*p_ekf_var(3,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, p_liekf(3,:), 3*p_liekf_var(3,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-300,100])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
print('position_noise_mean_cov', '-dpng')

figure;
subplot(311)
hold('on')
plot(t, theta_gt(1,:), 'k--', 'LineWidth', 2);
plot(t, theta_ekf(1,:), 'g', 'LineWidth', 1);
plot(t, theta_liekf(1,:), 'r', 'LineWidth', 1);
%plot(t, theta_ekf(1,:)+3*theta_ekf_var(1,:), 'b', 'LineWidth', 1);
%plot(t, theta_ekf(1,:)-3*theta_ekf_var(1,:), 'b', 'LineWidth', 1);
%plot(t, theta_liekf(1,:)+3*theta_liekf_var(1,:), 'm', 'LineWidth', 1);
%plot(t, theta_liekf(1,:)-3*theta_liekf_var(1,:), 'm', 'LineWidth', 1);
shadedErrorBar(t, theta_ekf(1,:), 3*theta_ekf_var(1,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, theta_liekf(1,:), 3*theta_liekf_var(1,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
title("Theta");
subplot(312)
hold('on')
plot(t, theta_gt(2,:), 'k--', 'LineWidth', 2)
%plot(t, theta_ekf(2,:), 'g', 'LineWidth', 1);
%plot(t, theta_liekf(2,:), 'r', 'LineWidth', 1)
%plot(t, theta_ekf(2,:)+3*theta_ekf_var(2,:), 'b', 'LineWidth', 1);
%plot(t, theta_ekf(2,:)-3*theta_ekf_var(2,:), 'b', 'LineWidth', 1);
%plot(t, theta_liekf(2,:)+3*theta_liekf_var(2,:), 'm', 'LineWidth', 1);
%plot(t, theta_liekf(2,:)-3*theta_liekf_var(2,:), 'm', 'LineWidth', 1);
shadedErrorBar(t, theta_ekf(2,:), 3*theta_ekf_var(2,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, theta_liekf(2,:), 3*theta_liekf_var(2,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
subplot(313)
hold('on')
plot(t, theta_gt(3,:),  'k--', 'LineWidth', 2)
%plot(t, theta_ekf(3,:), 'g', 'LineWidth', 1);
%plot(t, theta_liekf(3,:), 'r', 'LineWidth', 1)
%plot(t, theta_ekf(3,:)+3*theta_ekf_var(3,:), 'b', 'LineWidth', 1);
%plot(t, theta_ekf(3,:)-3*theta_ekf_var(3,:), 'b', 'LineWidth', 1);
%plot(t, theta_liekf(3,:)+3*theta_liekf_var(3,:), 'm', 'LineWidth', 1);
%plot(t, theta_liekf(3,:)-3*theta_liekf_var(3,:), 'm', 'LineWidth', 1);
shadedErrorBar(t, theta_ekf(3,:), 3*theta_ekf_var(3,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, theta_liekf(3,:), 3*theta_liekf_var(3,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
print('theta_noise_mean_cov', '-dpng')


% -------------------------------------------------------------------------
% Plot position and theta data to visualize
% the operation of the filter
figure;
subplot(311)
hold('on')
plot(t, zeros(size(p_gt(1,:))), 'k--', 'LineWidth', 2);
shadedErrorBar(t, zeros(size(p_ekf(1,:))), 3*p_ekf_var(1,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, zeros(size(p_liekf(1,:))), 3*p_liekf_var(1,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-200,200])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
title("Position 3 Standard Deviations");
subplot(312)
hold('on')
plot(t, zeros(size(p_gt(2,:))), 'k--', 'LineWidth', 2);
shadedErrorBar(t, zeros(size(p_ekf(1,:))), 3*p_ekf_var(2,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, zeros(size(p_liekf(1,:))), 3*p_liekf_var(2,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-200,200])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
subplot(313)
hold('on')
plot(t, zeros(size(p_gt(3,:))), 'k--', 'LineWidth', 2);
shadedErrorBar(t, zeros(size(p_ekf(1,:))), 3*p_ekf_var(3,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, zeros(size(p_liekf(1,:))), 3*p_liekf_var(3,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-200,200])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
print('position_noise_cov', '-dpng')

figure;
subplot(311)
hold('on')
plot(t, zeros(size(theta_gt(1,:))), 'k--', 'LineWidth', 2);
shadedErrorBar(t, zeros(size(theta_ekf(1,:))), 3*theta_ekf_var(1,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, zeros(size(theta_liekf(1,:))), 3*theta_liekf_var(1,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
title("Theta 3 Standard Deviations");
subplot(312)
hold('on')
plot(t, zeros(size(theta_gt(2,:))), 'k--', 'LineWidth', 2)
shadedErrorBar(t, zeros(size(theta_ekf(2,:))), 3*theta_ekf_var(2,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, zeros(size(theta_liekf(2,:))), 3*theta_liekf_var(2,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
subplot(313)
hold('on')
plot(t, zeros(size(theta_gt(3,:))),  'k--', 'LineWidth', 2)
shadedErrorBar(t, zeros(size(theta_ekf(3,:))), 3*theta_ekf_var(3,:), 'lineProps', {'g', 'LineWidth', 1})
shadedErrorBar(t, zeros(size(theta_liekf(3,:))), 3*theta_liekf_var(3,:), 'lineProps', {'r', 'LineWidth', 1})
axis([0,2,-7,7])
legend('Ground Truth', 'EKF', 'LIEKF', 'location', 'eastoutside')
print('theta_noise_cov', '-dpng')

