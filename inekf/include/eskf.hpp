
#pragma once
#include <Eigen/Dense>
#include <chrono>
#include <tuple>

namespace iekf
{

using Matrix3d = Eigen::Matrix<double, 3, 3>;
using Vector3d = Eigen::Matrix<double, 3, 1>;
using Quaternion = Eigen::Quaternion<double>;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using Timestamp = std::chrono::time_point<std::chrono::system_clock,
    std::chrono::duration<double>>;
using Seconds = std::chrono::duration<double, std::ratio<1, 1>>;

class ESKF
{
public:
    ESKF();
    // Reset time and origin of filter
    void resetFilter(const Timestamp& time, const Vector3d& origin_lla);

    // Add an IMU measurement to the filter
    void addImu(const Timestamp& timestamp, const Eigen::Vector3d& acc,
        const Eigen::Vector3d& gyro);

    // Add a gps measurement in LLA form
    void addGps(const Timestamp& timestamp, const Eigen::Vector3d& gps);

    // Set the origin
    // bool is set to indicate to other
    // parts of the program that the origin
    // is set and does not need to be modified
    void set_origin(Eigen::Vector3d gps_lla)
    {
        origin_ = gps_lla;
        origin_set_ = true;
    }

    // return position, velocity and attitude tuple
    std::tuple<Vector3d&, Vector3d&, Quaternion&> getState()
    {
        return std::tie(position, velocity, attitude);
    }

    // Return the covariance Sigma
    const MatrixXd& Sigma() const
    {
        return P;
    }

private:
    // Helper functions
    Matrix3d skewSymmetric(const Vector3d& v);
    Matrix3d expMap(const Vector3d& omega);

private:
    // State vectors
    Vector3d position;    // [m] in global frame
    Vector3d velocity;    // [m/s] in global frame
    Quaternion attitude;  // rotation from body to global frame
    Vector3d accel_bias;  // [m/s^2]
    Vector3d gyro_bias;   // [rad/s]

    // Error state (15x1)
    VectorXd error_state;

    // Covariance matrices
    MatrixXd P;      // Error state covariance (15x15)
    MatrixXd Q;      // Process noise covariance
    MatrixXd R_gps;  // GPS measurement noise covariance

    // Earth constants
    const double g = 9.81;  // gravity [m/s^2]
    const Matrix3d I3 = Matrix3d::Identity();

    // State mean and covaraiance variables
    // and timestamp
    Timestamp time_;               ///< filter time
    Timestamp time_last_predict_;  ///< saved time of last predict step
    Timestamp time_last_gps_;      ///< saved time for last gps addition

    Vector3d origin_;          ///< origin coordinates in lla
    bool origin_set_ = false;  ///< indicates if gps origin set
};

}  // namespace iekf