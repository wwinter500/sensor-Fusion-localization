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
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.
  
  Hint: one or more values initialized above might be wildly off...
  */
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //state diemention
  n_x_ = 5;
  //argument state dimention
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(1.);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

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

  //
  time_us_ = 0;
  lambda_ = 3 - n_aug_;

  //lidar noise covirance matrix
  _lidar_R_ = MatrixXd(2, 2);
  _lidar_R_ <<    std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;

  //radar noise covirance matrix
  _radar_R_ = MatrixXd(3,3);
  _radar_R_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;

  //set weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    weights_(i) = (double)(0.5/(n_aug_+lambda_));
  }

  // augmented and predicted sigma points
  xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1);
  xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);

  is_initialized_ = false;

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
  if(!is_initialized_)//if uninitialed,
  {
    initialization(meas_package);
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }

  if(meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
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

  //generate initial augmented sigma points
  GenerateAugmentedSigmaPoint();

  PredictAugmentedSigmaPoint(delta_t);

  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 * update x_ and P_ from K+1 | K ti K+1 | K+1
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z  = 2;//update px and py;

  VectorXd z_raw = meas_package.raw_measurements_;


  //create sigma point for in measurement space
  MatrixXd z_sig_ = MatrixXd(n_z, 2*n_aug_ + 1);

  //transform sigma points to measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //2n+1 simga points

    // extract values for better readibility
    // measurement model
    z_sig_(0,i) = (double)xsig_pred_(0,i);        //r
    z_sig_(1,i) = (double)xsig_pred_(1,i);        //phi
  }

  //mean predicted measurement z_pred and convariance matrix 
  VectorXd z_mean = VectorXd(n_z);
  z_mean.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_mean = z_mean + weights_(i) * z_sig_.col(i);
  }

  MatrixXd z_convariance = MatrixXd(n_z,n_z);
  z_convariance.fill(0.);

  for(int i = 0; i < 2*n_aug_+ 1; i++)
  {
    VectorXd z_diff = z_sig_.col(i) - z_mean;
    z_convariance = z_convariance + weights_(i) * z_diff * z_diff.transpose(); 
  }

  //add noise covariance measurement 
  z_convariance = z_convariance + _lidar_R_;

  //matrix for cross correlation
  MatrixXd cc = MatrixXd(n_x_, n_z);
  cc.fill(0.);

  for(int i = 0; i < 2* n_aug_ + 1; i++)
  {
    VectorXd z_diff = z_sig_.col(i) - z_mean;

    VectorXd x_diff = xsig_pred_.col(i) - x_;

    //angle normalization
    //while (x_diff(3)> M_PI || x_diff(3)<-M_PI) 
    NormalizeAngle(x_diff(3));
    
    cc = cc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = cc * z_convariance.inverse();

  //residual
  VectorXd z_diff = z_raw - z_mean;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*z_convariance*K.transpose();

  //calculate NIS
  NIS_Lidar = (z_raw - z_mean).transpose() * z_convariance.inverse() * (z_raw - z_mean);

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
  int n_z = 3;
  VectorXd z_raw = meas_package.raw_measurements_;

  //create matrix for sigma points in measurement space
  MatrixXd z_sig_ = MatrixXd(n_z, 2*n_aug_ +  1);

  //transform sigma points to measurement space
  double eps = 1.0e-6;
  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
    const double p_x = xsig_pred_(0,i);
    const double p_y = xsig_pred_(1,i);
    const double v  = xsig_pred_(2,i);
    const double yaw = xsig_pred_(3,i);

    const double v1 = cos(yaw)*v;
    const double v2 = sin(yaw)*v;

    z_sig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    z_sig_(1,i) = atan2(p_y,p_x);                                 //phi
    z_sig_(2,i) = (p_x*v1 + p_y*v2 ) / std::max(z_sig_(0,i), eps);                 //r_dot
  }

  //mean predicted measurement 
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * z_sig_.col(i);
  }

  //measurement covariance matrix 
  MatrixXd z_convariance = MatrixXd(n_z,n_z);
  z_convariance.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //2n+1 simga points
    //residual
    VectorXd z_diff = z_sig_.col(i) - z_pred;

    //angle normalization
    //while (z_diff(1)> M_PI || z_diff(1)<-M_PI) 
    NormalizeAngle(z_diff(1));

    z_convariance = z_convariance + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  z_convariance = z_convariance + _radar_R_;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //2n+1 simga points

    //residual
    VectorXd z_diff = z_sig_.col(i) - z_pred;
    //angle normalization
    
    //while (z_diff(1)> M_PI || z_diff(1)<-M_PI) 
   	NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = xsig_pred_.col(i) - x_;
    //angle normalization
    //while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    //while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    NormalizeAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * z_convariance.inverse();

  //residual
  VectorXd z_diff = z_raw - z_pred;

  //angle normalization
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  NormalizeAngle(z_diff(1));


  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*z_convariance*K.transpose();


  //calculate NIS
  NIS_Radar = (z_raw - z_pred).transpose() * z_convariance.inverse() * (z_raw - z_pred);
}

void UKF::initialization(MeasurementPackage meas_package)
{
    //initialize state vector with the first measurement
    time_us_ = meas_package.timestamp_;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      const double rho = meas_package.raw_measurements_[0];
      const double phi = meas_package.raw_measurements_[1];
      const double px = rho*cos(phi);
      const double py = rho*sin(phi);

      x_ << px,
            py,
            fabs(meas_package.raw_measurements_[2]),
            0,
            0;
    }

    //initialize covirance matrix to identity matrix
    P_.fill(0.);
    for(int i=0;i<n_x_;i++) P_.diagonal()[i]=1.;

    is_initialized_ = true;
}

void UKF::GenerateAugmentedSigmaPoint()
{
  //generated augmented sigma point
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd p_aug = MatrixXd(n_aug_, n_aug_);

  //initial x_aug
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  //initial p_aug
  p_aug.fill(0.);
  p_aug.topLeftCorner(5,5) = P_;
  p_aug(5,5) = std_a_*std_a_;
  p_aug(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd L = p_aug.llt().matrixL();

  //create augmented sigma points
  
  xsig_aug_.col(0)  = x_aug;

  double coff = sqrt(lambda_ + n_aug_);
  for (int i = 0; i< n_aug_; i++)
  {
    xsig_aug_.col(i+1)       = x_aug + coff * L.col(i);
    xsig_aug_.col(i+1+n_aug_) = x_aug - coff * L.col(i);
  }
}

void UKF::PredictAugmentedSigmaPoint(double delta_t)
{
  for (int i = 0; i< 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    const double p_x = xsig_aug_(0,i);
    const double p_y = xsig_aug_(1,i);
    const double v = xsig_aug_(2,i);
    const double yaw = xsig_aug_(3,i);
    const double yawd = xsig_aug_(4,i);
    const double nu_a = xsig_aug_(5,i);
    const double nu_yawdd = xsig_aug_(6,i);

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
    xsig_pred_(0,i) = px_p;
    xsig_pred_(1,i) = py_p;
    xsig_pred_(2,i) = v_p;
    xsig_pred_(3,i) = yaw_p;
    xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance()
{
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = xsig_pred_.col(i) - x_;
    //angle normalization
    //while (x_diff(3)> M_PI || x_diff(3)<-M_PI) 
    NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}
