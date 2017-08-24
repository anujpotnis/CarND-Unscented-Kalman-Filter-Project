#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;
    
    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;
    
    // initial state vector
    x_ = VectorXd(5);
    
    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.2;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.2        ;
    
    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;
    
    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;
    
    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;
    
    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;
    
    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;
    
    /**
     TODO:
     
     Complete the initialization. See ukf.h for other member properties.
     
     Hint: one or more values initialized above might be wildly off...
     */
    
    is_initialized_ = false;
    n_x_ = 5;
    x_ << 1,1,0,0,0;
    P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
    -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
    0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
    -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
    -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
    
    
}



UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
     TODO:
     
     Complete this function! Make sure you switch between lidar and radar
     measurements.
     */
    if(!is_initialized_) {
        previous_timestamp_ = meas_package.timestamp_;
        if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            
        }
        
        if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /*
             x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
             */
            x_ <<   5.7441,
            1.3800,
            2.2049,
            0.5015,
            0.3528;
        }
        is_initialized_ = true;
        return;
    }
    
    //compute the time elapsed between the current and previous measurements
    float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;
    
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    
    
    Prediction(dt);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    /**
     TODO:
     
     Complete this function! Estimate the object's location. Modify the state
     vector, x_. Predict sigma points, the state, and the state covariance matrix.
     */
    
    // dimension after augmentation
    n_aug_ = 7;
    
    //define spreading parameter
    double lambda = 3 - n_aug_;
    
    //create augmented mean vector x_ (size 7)
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug << x_,0,0;
    
    //create augmented state covariance P_ (size 7x7)
    Eigen::Matrix<double, 2, 2> Q_;
    Q_ << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
    
    MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
    P_aug.topLeftCorner(P_.rows(), P_.cols()) = P_;
    P_aug.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;
    // square root of P
    MatrixXd A = MatrixXd::Zero(n_aug_, n_aug_);
    
    A = P_aug.llt().matrixL();
    //create sigma point matrix Xsig_aug (size 7x15)
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    //set first column of sigma point matrix
    Xsig_aug.col(0)  = x_aug;
    
    //set remaining sigma points
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda+n_aug_) * A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * A.col(i);
    }
    std::cout << "Augmented Matrix" << std::endl;
    std::cout << Xsig_aug << std::endl;
    
    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
    
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        
        //write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }
    std::cout << "  " << std::endl;
    std::cout << "Sigma Prediction Matrix" << std::endl;
    std::cout << Xsig_pred << std::endl;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
     TODO:
     
     Complete this function! Use lidar data to update the belief about the object's
     position. Modify the state vector, x_, and covariance, P_.
     
     You'll also need to calculate the lidar NIS.
     */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
     TODO:
     
     Complete this function! Use radar data to update the belief about the object's
     position. Modify the state vector, x_, and covariance, P_.
     
     You'll also need to calculate the radar NIS.
     */
}
