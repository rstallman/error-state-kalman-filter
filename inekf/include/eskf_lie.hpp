#pragma once

#include <Eigen/Dense>
#include <chrono>
#include <tuple>

// Matrix Lie implementation of the error state Kalman filter (ESKF)
namespace iekf
{

// Define a bunch of commonly used types to save typing
// and make the files more clear. These are consistent
// with Eigen naming conventions.
using Matrix3d = Eigen::Matrix<double, 3, 3>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector3d = Eigen::Matrix<double, 3, 1>;
using MatrixXd = Eigen::MatrixXd;

using Vector6d = Eigen::Matrix<double, 6, 1>;
using Matrix15d = Eigen::Matrix<double, 15, 15>;
using Vector15d = Eigen::Matrix<double, 15, 1>;
using Rotation = Matrix3d;

// Shorthand chrono types
using Timestamp = std::chrono::time_point<std::chrono::system_clock,
    std::chrono::duration<double>>;
using Seconds = std::chrono::duration<double, std::ratio<1, 1>>;

class ESKF_LIE
{
public:
    ESKF_LIE();

    // Reset time and origin of filter
    void resetFilter(const Timestamp& time, const Vector3d& origin_lla);

    // Add an IMU measurement to the filter
    void addImu(
        const Timestamp& timestamp, const Vector3d& acc, const Vector3d& gyro);

    // Add a GPS measurement in LLA form
    void addGps(const Timestamp& timestamp, Vector3d& gps);

    // Set the origin
    void set_origin(Vector3d gps_lla)
    {
        origin_ = gps_lla;
        origin_set_ = true;
    }

    // Return position, velocity and attitude tuple
    std::tuple<Vector3d&, Vector3d&, Rotation&> getState()
    {
        return std::tie(position, velocity, attitude);
    }

    // Return the covariance Sigma
    const Matrix15d& Sigma() const
    {
        return P;
    }

private:
    // Perform the prediction step of the filter
    void prediction(const Eigen::Vector3d& acc, const Eigen::Vector3d& gyro,
        std::chrono::duration<double> dt);
    // Perform the correction step of the filter using gps coordinates
    void correction(const Eigen::Vector3d& gps_lla);
    // Matrix lie algebra helper functions
    Eigen::Matrix3d skewSymmetric(const Vector3d& v);
    Eigen::Matrix3d expMap(const Vector3d& omega);

private:
    Vector3d position;      ///< Position vector
    Vector3d velocity;      ///< Velocity vector
    Rotation attitude;      ///< Attitude in Rotation matrix form
    Vector3d accel_bias;    ///< Accelerometer bias
    Vector3d gyro_bias;     ///< Gyroscope bias
    Vector15d error_state;  ///< Error state vector
    Matrix15d P;            ///< State covariance matrix
    Matrix15d Q;            ///< Process noise covariance matrix
    Matrix6d R_gps;         ///< GPS measurement noise covariance matrix

    Timestamp time_;               ///< Current time
    Timestamp time_last_predict_;  ///< Last prediction time
    Timestamp time_last_gps_;      ///< Last GPS measurement time

    Eigen::Vector3d origin_;   ///< Origin coordinates in LLA
    bool origin_set_ = false;  ///< Indicates if GPS origin is set

    // Earth constants
    const double g = 9.81;  // gravity [m/s^2]
    const Matrix3d I3 = Matrix3d::Identity();
};
}  // namespace iekf
