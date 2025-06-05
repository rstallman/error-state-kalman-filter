#include <eskf_lie.hpp>
#include <utils.hpp>

namespace iekf
{

ESKF_LIE::ESKF_LIE()
{
    position.setZero();
    velocity.setZero();
    attitude.setIdentity();
    accel_bias.setZero();
    gyro_bias.setZero();
    error_state.setZero();

    // Initialize covariance matrices
    P.setIdentity();
    P.block<3, 3>(0, 0) *= 1.0;      // position
    P.block<3, 3>(3, 3) *= 0.1;      // velocity
    P.block<3, 3>(6, 6) *= 0.01;     // attitude
    P.block<3, 3>(9, 9) *= 0.001;    // accel bias
    P.block<3, 3>(12, 12) *= 0.001;  // gyro bias

    Q.setIdentity();
    Q.block<3, 3>(0, 0) *= 0.01;    // process noise (position)
    Q.block<3, 3>(3, 3) *= 0.01;    // process noise (velocity)
    Q.block<3, 3>(6, 6) *= 0.001;   // process noise (attitude)
    Q.block<3, 3>(9, 9) *= 1e-5;    // process noise (accel bias)
    Q.block<3, 3>(12, 12) *= 1e-5;  // process noise (gyro bias)

    R_gps.setIdentity();
    R_gps.block<3, 3>(0, 0) *= 0.1;  // GPS position noise
    R_gps.block<3, 3>(3, 3) *= 0.5;  // GPS velocity noise
}

void ESKF_LIE::addImu(
    const Timestamp& timestamp, const Vector3d& acc, const Vector3d& gyro)
{
    // Set the time to the current timestamp
    time_ = timestamp;

    // Calculate dt based on last observed time and current timestamp
    auto duration = timestamp - time_last_predict_;
    time_last_predict_ = timestamp;

    prediction(acc, gyro, duration);
}

void ESKF_LIE::addGps(const Timestamp& timestamp, Vector3d& gps)
{
    correction(gps);
    time_last_gps_ = timestamp;
}

void ESKF_LIE::prediction(const Eigen::Vector3d& acc,
    const Eigen::Vector3d& gyro, std::chrono::duration<double> duration)
{
    const double dt = duration.count();

    if (dt <= 0) return;

    // Extract measurements (correcting for bias)
    Vector3d acc_corrected = acc - accel_bias;
    Vector3d gyro_corrected = gyro - gyro_bias;

    // Rotation matrix from quaternion
    Matrix3d R = attitude;

    // Update position and velocity using IMU measurements
    position += velocity * dt +
                0.5 * (R * acc_corrected + Vector3d(0, 0, -g)) * dt * dt;
    velocity += (R * acc_corrected + Vector3d(0, 0, -g)) * dt;

    // Update attitude using gyroscope measurements
    Matrix3d dR = expMap(gyro_corrected * dt);
    R = R * dR;
    attitude = R;

    // Update covariance
    Matrix15d F_x = Matrix15d::Identity();
    F_x.block<3, 3>(0, 3) = Matrix3d::Identity() * dt;
    F_x.block<3, 3>(3, 6) = -1.0 * R * skewSymmetric(acc_corrected) * dt;
    F_x.block<3, 3>(3, 9) = -1.0 * R * dt;
    F_x.block<3, 3>(6, 6) = dR.transpose();
    F_x.block<3, 3>(6, 12) = -Matrix3d::Identity() * dt;
    P = F_x * P * F_x.transpose() + Q;
}

// Perform the correction step of the filter using gps coordinates
void ESKF_LIE::correction(const Eigen::Vector3d& gps_lla)
{
    // Set origin if this is the first measurement
    if (!origin_set_) set_origin(gps_lla);

    // Convert from LLA (spherical) to ENU (cartesian)
    auto gps_enu = lla_to_enu(gps_lla, origin_);

    // Measurement matrix (H)
    MatrixXd H = MatrixXd::Zero(6, 15);
    H.block<3, 3>(0, 0) = Matrix3d::Identity();  // position
    H.block<3, 3>(3, 3) = Matrix3d::Identity();  // velocity

    // Measurement residual
    Vector6d z(6);
    z.head(3) = gps_enu - position;
    z.tail(3) = Vector3d::Zero() - velocity;  // Assuming zero velocity from GPS

    // Kalman gain, K is 15 x 6
    MatrixXd K = P * H.transpose() * (H * P * H.transpose() + R_gps).inverse();

    // Update error state and covariance (15 x 1)
    error_state = K * z;
    P = (Matrix15d::Identity(15, 15) - K * H) * P;

    // Inject error state into nominal state
    position += error_state.segment<3>(0);
    velocity += error_state.segment<3>(3);

    // Attitude correction (using exponential map)
    Vector3d delta_theta = error_state.segment<3>(6);
    attitude = attitude * expMap(delta_theta);

    // Update biases
    accel_bias += error_state.segment<3>(9);
    gyro_bias += error_state.segment<3>(12);

    // Reset error state
    error_state.setZero();
};

// Helper functions
Matrix3d ESKF_LIE::skewSymmetric(const Vector3d& v)
{
    Matrix3d S;
    S << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return S;
}

Matrix3d ESKF_LIE::expMap(const Vector3d& omega)
{
    double theta = omega.norm();
    if (theta < 1e-12) return I3;

    Matrix3d omega_hat = skewSymmetric(omega / theta);
    return I3 + sin(theta) * omega_hat +
           (1 - cos(theta)) * omega_hat * omega_hat;
}
void ESKF_LIE::resetFilter(const Timestamp& time, const Vector3d& origin_lla)
{
    origin_ = origin_lla;
    origin_set_ = true;
    time_ = time;
    time_last_predict_ = time;
}

}  // namespace iekf