#include <DataLoader.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utils.hpp>
#include <vector>

using namespace Eigen;
// Shorthand chrono types
using Timestamp = std::chrono::time_point<std::chrono::system_clock,
    std::chrono::duration<double>>;
using Seconds = std::chrono::duration<double, std::ratio<1, 1>>;

size_t count(DataLoader::OutDataType in)
{
    size_t out{0};
    for (auto i{in.first}; i != in.second; ++i, ++out)
        ;
    return out;
}

struct ImuMeas
{
    Vector3d accel;
    Vector3d gyro;
    bool accel_set = false;
    bool gyro_set = false;
};

class AIESKF
{
private:
    // State vectors
    Vector3d position;     // [m] in global frame
    Vector3d velocity;     // [m/s] in global frame
    Quaterniond attitude;  // rotation from body to global frame
    Vector3d accel_bias;   // [m/s^2]
    Vector3d gyro_bias;    // [rad/s]

    // Error state (15x1)
    VectorXd error_state;

    // Covariance matrices
    MatrixXd P;      // Error state covariance (15x15)
    MatrixXd Q;      // Process noise covariance
    MatrixXd R_gps;  // GPS measurement noise covariance

    // Time management
    bool initialized;

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

public:
    AIESKF()
        : error_state(15), P(15, 15), Q(15, 15), R_gps(6, 6), initialized(false)
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
        initialize();
    }

    void initialize()
    {
        position = Vector3d(0, 0, 0);
        velocity.setZero();
        attitude.setIdentity();
        accel_bias.setZero();
        gyro_bias.setZero();
        error_state.setZero();
        initialized = true;
    }

    // Set the origin
    // bool is set to indicate to other
    // parts of the program that the origin
    // is set and does not need to be modified
    void set_origin(Eigen::Vector3d gps_lla)
    {
        origin_ = gps_lla;
        origin_set_ = true;
    }

    void addImu(const Timestamp &timestamp, const Eigen::Vector3d &acc_imu,
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
        position +=
            velocity * dt + 0.5 * (R * acc + Vector3d(0, 0, -g)) * dt * dt;
        velocity += (R * acc + Vector3d(0, 0, -g)) * dt;

        // Quaternion integration
        Quaterniond dq;
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

    void addGps(const Timestamp &timestamp, const Vector3d &gps_lla)
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
        z.tail(3) =
            Vector3d::Zero() - velocity;  // Assuming zero velocity from GPS

        // Kalman gain
        MatrixXd K =
            P * H.transpose() * (H * P * H.transpose() + R_gps).inverse();

        // Update error state and covariance
        error_state = K * z;
        P = (MatrixXd::Identity(15, 15) - K * H) * P;

        // Inject error state into nominal state
        position += error_state.segment<3>(0);
        velocity += error_state.segment<3>(3);

        // Attitude correction (using exponential map)
        Vector3d delta_theta = error_state.segment<3>(6);
        Quaterniond delta_q;
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

    // Helper functions
    Matrix3d skewSymmetric(const Vector3d &v)
    {
        Matrix3d S;
        S << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
        return S;
    }

    Matrix3d expMap(const Vector3d &omega)
    {
        double theta = omega.norm();
        if (theta < 1e-12) return I3;

        Matrix3d omega_hat = skewSymmetric(omega / theta);
        return I3 + sin(theta) * omega_hat +
               (1 - cos(theta)) * omega_hat * omega_hat;
    }

    // Get current state
    void getState(Vector3d &pos, Vector3d &vel, Quaterniond &att)
    {
        pos = position;
        vel = velocity;
        att = attitude;
    }
};

int main()
{
    auto filepath =
        "/Users/lijinliang/UMICH_Rob530/invariant-ekf/inekf/dataset/AGZ_subset";
    // try to open file, fail otherwise
    std::ifstream test_open(filepath);
    if (!test_open.is_open())
    {
        std::cerr << "File not found or unable to be opened\n";
        std::cerr << "Please provide a correct path to AGZ subset\n";
        exit(1);
    }
    test_open.close();

    // Create data loader and imu measurement
    DataLoader dl(filepath, false);  // false to not use ground truth gps
    ImuMeas imu;

    // Create aiekf (filter by deepseek)
    AIESKF filter;

    int cnt = 0;
    std::ofstream fused_file("fused_eskf.txt", std::ios::trunc);

    while (!dl.complete())
    {
        // Get the next set of data
        auto next_data = dl.next();
        DataLoader::Timestamp ts;

        // Loop through multiple measurements in current data set
        // if needed.
        for (auto i{next_data.first}; i != next_data.second; ++i)
        {
            // Save the timestamp to put in the imu measurement at the bottom
            // of the loop.
            ts = i->first;

            // Test based on data type
            switch (i->second.dt)
            {
                case DataLoader::DataType::omega:
                    // Logic for building imu measurements from data loader
                    if (imu.gyro_set)
                    {
                        throw std::runtime_error("multiple gyro measurements");
                    }
                    else
                    {
                        imu.gyro = i->second.datum;
                        imu.gyro_set = true;
                    }
                    break;

                case DataLoader::DataType::accel:
                    // Logic for building imu measurements from data loader
                    if (imu.accel_set)
                    {
                        throw std::runtime_error("multiple accel measurements");
                    }
                    else
                    {
                        imu.accel = i->second.datum;
                        imu.accel_set = true;
                    }
                    break;

                case DataLoader::DataType::gps:
                    // Add gps measurement to filter
                    filter.addGps(ts, i->second.datum);
                    std::cout << "Added gps measurement at time: "
                              << ts.time_since_epoch().count() << "\n";
                    std::cout << "gps: " << i->second.datum.transpose() << "\n";
                    break;
                case DataLoader::DataType::accel_bias:
                    // do nothing
                    break;
            }
        }

        // Once imu measurement was built add it to the filter
        // Really the correct thing to do would be to rewrite
        // the data loader so it returns pre-built imu measurements...
        if (imu.accel_set && imu.gyro_set)
        {
            // add imu measurement to filter
            filter.addImu(ts, imu.accel, imu.gyro);

            // Reset so that new measurements can be added
            imu.accel_set = false;
            imu.gyro_set = false;
        }

        Vector3d pos, vel;
        Quaterniond att;
        filter.getState(pos, vel, att);

        auto time = ts;      // Timestamp of the current measurement
        double px = pos[0];  // First element
        double py = pos[1];  // Second element
        double pz = pos[2];  // Third element

        fused_file << std::fixed << std::setprecision(10)
                   << time.time_since_epoch().count() << " " << px << " " << py
                   << " " << pz << " " << att.x() << " " << att.y() << " "
                   << att.z() << " " << att.w() << std::endl;
        // std::cout << "cnt: " << cnt << "\n";
        // std::cout << "Position: " << pos.transpose() << " m\n";
        // std::cout << "Velocity: " << vel.transpose() << " m/s\n";
        // std::cout << "Attitude (wxyz): " << att.w() << " "
        //           << att.vec().transpose() << "\n\n";
        // if (cnt > 20000) break;
        // cnt++;
    }
}
