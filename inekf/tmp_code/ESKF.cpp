#include "ESKF.hpp"

#include <cmath>
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>
#include <utils.hpp>

#include "math_utils.hpp"

using namespace std::chrono;
using namespace Eigen;
using std::cos;
using std::pow;
using std::sin;

namespace practice
{
// Define a bunch of commonly used types to save typing
// and make the files more clear. These are consistent
// with Eigen naming conventions.
using Matrix18d = Eigen::Matrix<double, 18, 18>;
using Matrix12d = Eigen::Matrix<double, 12, 12>;
using Matrix9d = Eigen::Matrix<double, 9, 9>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector18d = Eigen::Matrix<double, 18, 1>;
using Matrix3d = Eigen::Matrix<double, 3, 3>;
using Vector3d = Eigen::Matrix<double, 3, 1>;
using Vector6d = Eigen::Matrix<double, 6, 1>;
using Quaternion = Eigen::Quaterniond;  // for attitude representation

// Shorthand chrono types
using Timestamp = std::chrono::time_point<std::chrono::system_clock,
    std::chrono::duration<double>>;
using Seconds = std::chrono::duration<double, std::ratio<1, 1>>;

ESKF::ESKF()
    : nominal_pos_(Vector3d::Zero()),
      nominal_vel_(Vector3d::Zero()),
      nominal_attitude_(Quaternion::Identity()),
      nominal_accel_bias_(Vector3d::Zero()),
      nominal_gyro_bias_(Vector3d::Zero()),
      P_(Matrix18d::Identity()),
      error_x_(Vector18d::Zero())
{
}

void ESKF::resetFilter(const Timestamp& time, const Vector3d& origin_lla)
{
    origin_ = origin_lla;
    origin_set_ = true;
    time_ = time;
    time_last_predict_ = time;
}

void ESKF::addImu(
    const Timestamp& timestamp, const Vector3d& acc, const Vector3d& gyro)
{
    // Set the time to the current timestamp
    time_ = timestamp;

    // Calculate dt based on last observed time and current timestamp
    auto duration = timestamp - time_last_predict_;
    time_last_predict_ = timestamp;

    const double dt = duration.count();
    const double dt2 = dt * dt;

    Matrix3d R = nominal_attitude_.toRotationMatrix();
    Vector3d accel_diff = acc - nominal_accel_bias_;
    Vector3d gyro_diff = gyro - nominal_gyro_bias_;

    nominal_pos_ +=
        nominal_vel_ * dt + 0.5 * (R * accel_diff + nominal_gravity_) * dt2;
    nominal_vel_ += (R * accel_diff + nominal_gravity_) * dt;
    Matrix3d dR = rotationVectorToRotationMatrix(gyro_diff * dt);
    nominal_attitude_ = Quaternion(R * dR);
    nominal_attitude_.normalize();

    Matrix18d F_x_ = Matrix18d::Identity();
    F_x_.block<3, 3>(0, 3) = Matrix3d::Identity() * dt;
    F_x_.block<3, 3>(3, 6) = -1.0 * R * vec3ToSkewSymmetric(accel_diff) * dt;
    F_x_.block<3, 3>(3, 9) = -1.0 * R * dt;
    F_x_.block<3, 3>(3, 15) = Matrix3d::Identity() * dt;
    F_x_.block<3, 3>(6, 6) = dR.transpose();
    F_x_.block<3, 3>(6, 12) = -1.0 * Matrix3d::Identity() * dt;
    Eigen::Matrix<double, 18, 12> F_i_ = Eigen::Matrix<double, 18, 12>::Zero();
    F_i_.block<12, 12>(3, 0) = Matrix12d::Identity();

    // imu parameters
    double accel_noise = 0.002;
    double accel_bias_stability = 0.02;
    double gyro_noise = 0.05;
    double gyro_bias_stability = 1.0;
    double update_rate = 100;

    Matrix12d Q_i_ = Matrix12d::Identity();
    Q_i_.block<3, 3>(0, 0) =
        std::pow(accel_noise, 2) * update_rate * Matrix3d::Identity();
    Q_i_.block<3, 3>(3, 3) =
        std::pow(gyro_noise, 2) * update_rate * Matrix3d::Identity();
    Q_i_.block<3, 3>(6, 6) =
        std::pow(accel_bias_stability * update_rate, 2) * Matrix3d::Identity();
    Q_i_.block<3, 3>(9, 9) =
        std::pow(gyro_bias_stability * update_rate, 2) * Matrix3d::Identity();

    P_ = F_x_ * P_ * F_x_.transpose() + F_i_ * Q_i_ * F_i_.transpose() * dt2;
}

void ESKF::addGps(const Timestamp& timestamp, const Eigen::Vector3d& gps_lla)
{
    // Set origin if this is the first measurement
    if (!origin_set_) set_origin(gps_lla);

    // Convert from lla (spherical) to enu (cartesian)
    auto gps = lla_to_enu(gps_lla, origin_);

    // TODO: fix the following
    Matrix3d position_covariance;  // position covariance
    position_covariance << 4.6778e3, 1.9437e3, 0.0858e3, 1.9437e3, 11.5621e3,
        5.8445e3, 0.0858e3, 5.8445e3, 22.4051e3;
    Matrix6d V_ = Matrix6d::Identity();
    V_.block<3, 3>(0, 0) = position_covariance;

    Eigen::Matrix<double, 19, 18> J_true_error_ =
        Eigen::Matrix<double, 19, 18>::Zero();
    J_true_error_.block<6, 6>(0, 0) = Matrix6d::Identity();
    J_true_error_.block<9, 9>(10, 9) = Matrix9d::Identity();
    J_true_error_.block<4, 3>(6, 6) = computeQuatJacobiToErrorQuat();
    J_true_error_.block<4, 3>(6, 6) = computeQuatJacobiToErrorQuat();

    Eigen::Matrix<double, 6, 19> H_x_ = Eigen::Matrix<double, 6, 19>::Zero();
    H_x_.block<6, 6>(0, 0) = Matrix6d::Identity();

    Eigen::Matrix<double, 6, 18> H_ = H_x_ * J_true_error_;
    Eigen::Matrix<double, 18, 6> H_t = H_.transpose();

    Eigen::Matrix<double, 18, 6> K_ = P_ * H_t * (H_ * P_ * H_t + V_).inverse();

    Vector6d z_ = Vector6d::Zero();

    z_.block<3, 1>(0, 0) = gps;
    // z_.block<3, 1>(3, 0) = gnssData->linearVelocity;  // do not use velocity

    error_x_ = K_ * (z_ - computeHx());
    P_ = (Matrix18d::Identity() - K_ * H_) * P_;

    injectErrorToNominal();
    resetErrorState();
    time_last_gps_ = timestamp;
}

Eigen::Matrix<double, 4, 3> ESKF::computeQuatJacobiToErrorQuat()
{
    Eigen::Matrix<double, 4, 3> quat_true_error =
        Eigen::Matrix<double, 4, 3>::Zero();
    quat_true_error.block<3, 3>(1, 0) = Matrix3d::Identity() * 0.5;
    quat_true_error =
        quatToLeftProductMatrix(nominal_attitude_) * quat_true_error;
    return quat_true_error;
}

Vector6d ESKF::computeHx()
{
    Vector6d hx;
    hx.block<3, 1>(0, 0) = nominal_pos_;
    hx.block<3, 1>(3, 0) = nominal_vel_;
    return hx;
}

void ESKF::injectErrorToNominal()
{
    nominal_pos_ += error_x_.block<3, 1>(0, 0);
    nominal_vel_ += error_x_.block<3, 1>(3, 0);
    nominal_attitude_ = nominal_attitude_ *
                        rotationVectorToQuaternion(error_x_.block<3, 1>(6, 0));
    nominal_attitude_.normalize();
    nominal_accel_bias_ += error_x_.block<3, 1>(9, 0);
    nominal_gyro_bias_ += error_x_.block<3, 1>(12, 0);
}

void ESKF::resetErrorState()
{
    Matrix18d G_ = Matrix18d::Identity();
    G_.block<3, 3>(6, 6) =
        Matrix3d::Identity() -
        0.5 * vec3ToSkewSymmetric(error_x_.block<3, 1>(6, 0));
    P_ = G_ * P_ * G_.transpose();
    error_x_.setZero();
}

}  // namespace practice