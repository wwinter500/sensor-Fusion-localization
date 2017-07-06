#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;

  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::UpdateMatrix(const VectorXd &y)
{
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  KalmanFilter::UpdateMatrix(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z0) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  const float PI = 3.1415927;

  VectorXd z(3);
  z << z0(0), z0(1), z0(2);

  if (z(1) > PI) 
    z(1) = z(1) - 2*PI;
  else if (z(1) < -PI)
    z(1) = 2*PI + z(1);

  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float phi = atan2(x_(1), x_(0)); 
  float rho_dot = 0.0;
  if (fabs(rho) > 1e-6) { 
    rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
  }
  VectorXd hx(3);
  hx << rho, phi, rho_dot;

  // Calculate y vector and normalize y(1), i.e. make -pi < y(1) < pi
  VectorXd y = z - hx;
  if (y(1) > PI) {
    y(1) = y(1) - 2*PI;
  }
  else if (y(1) < -PI) {
    y(1) = 2*PI + y(1);
  }

  KalmanFilter::UpdateMatrix(y);
}
