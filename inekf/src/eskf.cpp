#include <eskf.hpp>

#include "utils.hpp"
namespace iekf
{
ESKF::ESKF() : error_state(15), P(15, 15), Q(15, 15), R_gps(6, 6)
{
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
                                     // initialize();

    position = Vector3d(0, 0, 0);
    velocity.setZero();
    attitude.setIdentity();
    accel_bias.setZero();
    gyro_bias.setZero();
    error_state.setZero();
}

// Add an IMU measurement (measurement as motion model) to the filter
void ESKF::addImu(const Timestamp &timestamp, const Eigen::Vector3d &acc_imu,
    const Eigen::Vector3d &gyro)
{
    // Set the time to the current timestamp
    time_ = timestamp;

    // Calculate dt based on last observed time and current timestamp
    auto duration = timestamp - time_last_predict_;
    time_last_predict_ = timestamp;

    const double dt = duration.count();

    if (dt <= 0) return;

    // Extract measurements (correcting for bias)
    Vector3d acc = acc_imu - accel_bias;
    Vector3d omega = gyro - gyro_bias;

    // Normalize quaternion
    attitude.normalize();

    // Rotation matrix from quaternion
    Matrix3d R = attitude.toRotationMatrix();

    // Nominal state prediction
    position += velocity * dt + 0.5 * (R * acc + Vector3d(0, 0, -g)) * dt * dt;
    velocity += (R * acc + Vector3d(0, 0, -g)) * dt;

    // Quaternion integration
    Quaternion dq;
    Vector3d theta = omega * dt;
    double theta_norm = theta.norm();
    if (theta_norm > 1e-12)
    {
        dq.w() = cos(theta_norm / 2);
        dq.vec() = sin(theta_norm / 2) * theta / theta_norm;
    }
    else
    {
        dq.w() = 1.0;
        dq.vec() = 0.5 * theta;
    }
    attitude = attitude * dq;

    // Error-state transition matrix (F)
    MatrixXd F = MatrixXd::Identity(15, 15);
    F.block<3, 3>(0, 3) = I3 * dt;
    F.block<3, 3>(3, 6) = -R * skewSymmetric(acc) * dt;
    F.block<3, 3>(3, 9) = -R * dt;
    F.block<3, 3>(6, 6) = R.transpose() * expMap(-omega * dt);
    F.block<3, 3>(6, 12) = -I3 * dt;

    // Update covariance
    P = F * P * F.transpose() + Q;

    time_last_predict_ = timestamp;
}

void ESKF::addGps(const Timestamp &timestamp, const Vector3d &gps_lla)
{
    // Set origin if this is the first measurement
    if (!origin_set_) set_origin(gps_lla);

    // Convert from lla (spherical) to enu (cartesian)
    auto gps = lla_to_enu(gps_lla, origin_);

    // Measurement matrix (H)
    MatrixXd H = MatrixXd::Zero(6, 15);
    H.block<3, 3>(0, 0) = I3;  // position
    H.block<3, 3>(3, 3) = I3;  // velocity

    // Measurement residual
    VectorXd z(6);
    z.head(3) = gps - position;
    z.tail(3) = Vector3d::Zero() - velocity;  // Assuming zero velocity from GPS

    // Kalman gain
    MatrixXd K = P * H.transpose() * (H * P * H.transpose() + R_gps).inverse();

    // Update error state and covariance
    error_state = K * z;
    P = (MatrixXd::Identity(15, 15) - K * H) * P;

    // Inject error state into nominal state
    position += error_state.segment<3>(0);
    velocity += error_state.segment<3>(3);

    // Attitude correction (using exponential map)
    Vector3d delta_theta = error_state.segment<3>(6);
    Quaternion delta_q;
    double theta_norm = delta_theta.norm();
    if (theta_norm > 1e-12)
    {
        delta_q.w() = cos(theta_norm / 2);
        delta_q.vec() = sin(theta_norm / 2) * delta_theta / theta_norm;
    }
    else
    {
        delta_q.w() = 1.0;
        delta_q.vec() = 0.5 * delta_theta;
    }
    attitude = attitude * delta_q;

    // Update biases
    accel_bias += error_state.segment<3>(9);
    gyro_bias += error_state.segment<3>(12);

    // Reset error state
    error_state.setZero();

    time_last_gps_ = timestamp;
}

void ESKF::resetFilter(const Timestamp &time, const Vector3d &origin_lla)
{
    origin_ = origin_lla;
    origin_set_ = true;
    time_ = time;
    time_last_predict_ = time;
}

// Helper functions
Matrix3d ESKF::skewSymmetric(const Vector3d &v)
{
    Matrix3d S;
    S << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return S;
}

Matrix3d ESKF::expMap(const Vector3d &omega)
{
    double theta = omega.norm();
    if (theta < 1e-12) return I3;

    Matrix3d omega_hat = skewSymmetric(omega / theta);
    return I3 + sin(theta) * omega_hat +
           (1 - cos(theta)) * omega_hat * omega_hat;
}

}  // namespace iekf