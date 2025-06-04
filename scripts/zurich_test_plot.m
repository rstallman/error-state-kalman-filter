clear; close all;
addpath('filters');
addpath('helper');
addpath('thirdparty/shadedErrorBar');

% Load / process data
[T_X, omega, accel, accel_b, T_GPS, XYZ_GPS] = loadPoseGPS();
% test_N = length(omega); % Sets the number of IMU readings
test_N = 30000;
w = omega(1:test_N,:);
a = accel(1:test_N,:);
b_g = zeros(test_N,3);
b_a = accel_b(1:test_N,:);
t_x = T_X(1:test_N,:);

meas_used = T_GPS <= t_x(end);
t_gps = T_GPS(meas_used,:);
xyz_gps = XYZ_GPS(meas_used,:); 


%------------------------------------------------------------------------

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
% -------------------------------------------------------------------------
% Initialize filter

skew = @(u) [0 -u(3) u(2);
        u(3) 0 -u(1);
        -u(2) u(1) 0];

ekf = EKF();
inekf = LIEKF();

% Get first observation that happens after a prediction
obsid = 1;
while(t_gps(obsid) < t_x(1))
    obsid = obsid + 1;
end

pos_ekf = zeros(3,test_N);
pos_inekf = zeros(3,test_N);
std_ekf = zeros(3,test_N);
std_inekf = zeros(3,test_N);


for i = 2:test_N
    if i == 1
        dt = t_x; 
    else
        dt = t_x(i) - t_x(i - 1);
        
        %Assume gyro/IMU are basically synchronous
        ekf.prediction(w(i,:)',a(i,:)',dt);
        inekf.prediction(w(i,:)',a(i,:)',dt);
        
        
        %Measurement update
        if(i < test_N)
            if(t_gps(obsid) > t_x(i) && t_x(i+1) > t_gps(obsid))
                gps = [xyz_gps(obsid,1); xyz_gps(obsid,2); xyz_gps(obsid,3)];
                ekf.correction(gps);
                inekf.correction(gps);
                obsid = min(obsid + 1, length(t_gps));
            end
        end
        
        lieTocartesian(inekf)
        
        pos_ekf(:,i) = ekf.mu(1:3);
        variances = sqrt(diag(ekf.Sigma));
        std_ekf(:,i) = Log(eul2rotm(ekf.mu(1:3)));
        
        pos_inekf(:,i) = inekf.mu(1:3,5);
        vars = sqrt(diag(inekf.sigma_cart));
        std_inekf(:,i) = vars(7:9);
        
        if(mod(i,1000)==0)
           fprintf('Iteration: %d/%d\n',i,test_N); 
        end
    end
end

meas_used = T_GPS <= t_x(end);

% load gt
[~, ~, ~, ~, ~, x_gt, ~, y_gt, ~, z_gt] = loadGroundTruthAGL();
x_gt = x_gt - x_gt(1); y_gt = y_gt - y_gt(1); z_gt = z_gt - z_gt(1);
t_gt = linspace(0,T_X(end),length(x_gt));

% -------------------------------------------------------------------------
% traj plot
figure('DefaultAxesFontSize',14)
hold on;
plot3(XYZ_GPS(:,1), XYZ_GPS(:,2), XYZ_GPS(:,3),'b','LineWidth', 2);
plot3(x_gt, y_gt, z_gt,'--k','LineWidth', 4);
plot3(pos_ekf(1,:), pos_ekf(2,:), pos_ekf(3,:),'g','LineWidth', 2);
plot3(pos_inekf(1,:), pos_inekf(2,:), pos_inekf(3,:),'r','LineWidth', 2);
legend('gps', 'gt', 'EKF', 'InEKF', 'location', 'southeast')
hold off;
axis equal;

figure('DefaultAxesFontSize',14)
hold on;
plot3(XYZ_GPS(:,1), XYZ_GPS(:,2), XYZ_GPS(:,3),'b','LineWidth', 2);
plot3(x_gt, y_gt, z_gt,'--k','LineWidth', 4);
plot3(pos_inekf(1,:), pos_inekf(2,:), pos_inekf(3,:),'r','LineWidth', 2);
legend('gps', 'gt', 'InEKF', 'location', 'southeast')
hold off;
axis equal;

% -------------------------------------------------------------------------
% axis plot
figure;
subplot(3,1,1);
hold on;
plot(t_gps, XYZ_GPS(meas_used,1), 'b', 'LineWidth', 1);
plot(t_gt, x_gt, 'k--', 'LineWidth', 2);
% shadedErrorBar(T_X(1:test_N), pos_ekf(1,:), 3*std_ekf(1,:), 'lineProps', {'g', 'LineWidth', 1})
% shadedErrorBar(T_X(1:test_N), pos_inekf(1,:), 3*std_inekf(1,:), 'lineProps', {'r', 'LineWidth', 1})
legend('X_{GPS}','X_{GT}','X_{EKF}', 'X_{InEKF}', 'Location', 'eastoutside');
axis([0,T_X(test_N),-200,200])
%
subplot(3,1,2);
hold on;
plot(t_gps, XYZ_GPS(meas_used,2), 'b', 'LineWidth', 1);
plot(t_gt, y_gt, 'k--', 'LineWidth', 2);
% shadedErrorBar(T_X(1:test_N), pos_ekf(2,:), 3*std_ekf(2,:), 'lineProps', {'g', 'LineWidth', 1})
% shadedErrorBar(T_X(1:test_N), pos_inekf(2,:), 3*std_inekf(2,:), 'lineProps', {'r', 'LineWidth', 1})
legend('Y_{GPS}','Y_{GT}','Y_{EKF}', 'Y_{InEKF}', 'Location', 'eastoutside');
axis([0,T_X(test_N),-250,350])
%
subplot(3,1,3);
hold on;
plot(t_gps, XYZ_GPS(meas_used,3), 'b', 'LineWidth', 1);
plot(t_gt, z_gt, 'k--', 'LineWidth', 2);
% shadedErrorBar(T_X(1:test_N), pos_ekf(3,:), 3*std_ekf(3,:), 'lineProps', {'g', 'LineWidth', 1})
% shadedErrorBar(T_X(1:test_N), pos_inekf(3,:), 3*std_inekf(3,:), 'lineProps', {'r', 'LineWidth', 1})
legend('Z_{GPS}','Z_{GT}','Z_{EKF}', 'Z_{InEKF}', 'Location', 'eastoutside');
axis([0,T_X(test_N),-30,60])


function u = unskew(ux)
    u(1,1) = -ux(2,3);
    u(2,1) = ux(1,3);
    u(3,1) = -ux(1,2);
end

function w = Log(R)
    w = unskew(logm(R));
end