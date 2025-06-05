#pragma once

#include <chrono>
#include <eigen3/Eigen/Dense>
#include <tuple>

namespace practice
{
class ESKF
{
public:
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

public:
    // Default constructor initialize mu and Sigma to identity
    ESKF();
    ~ESKF(){};

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

    // Return the quaternion attitude
    Quaternion q() const
    {
        return nominal_attitude_;
    }

    // Return position vector p
    Vector3d p() const
    {
        return nominal_pos_;
    }

    // Return velocity vector v
    Vector3d v() const
    {
        return nominal_vel_;
    }

    // Return the covariance Sigma
    const Matrix18d& Sigma() const
    {
        return P_;
    }

private:
    void injectErrorToNominal();
    void resetErrorState();

    Eigen::Matrix<double, 4, 3> computeQuatJacobiToErrorQuat();
    Eigen::Matrix<double, 6, 1> computeHx();

    Vector3d nominal_pos_;         /// position mean
    Vector3d nominal_vel_;         /// velocity mean
    Quaternion nominal_attitude_;  /// Rotation mean
    Vector3d nominal_accel_bias_;
    Vector3d nominal_gyro_bias_;
    const Vector3d nominal_gravity_ =
        Vector3d(0, 0, -9.81);  /// gravity mean, could change with location

    Matrix18d P_;  // covariance of the error state

    Vector18d error_x_;  // error state vector

    // State mean and covaraiance variables
    // and timestamp
    Timestamp time_;               ///< filter time
    Timestamp time_last_predict_;  ///< saved time of last predict step
    Timestamp time_last_gps_;      ///< saved time for last gps addition

    Vector3d origin_;          ///< origin coordinates in lla
    bool origin_set_ = false;  ///< indicates if gps origin set
};

}  // namespace practice
