#include "kalman_filter.h"
#include <math.h>

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
    * predict the state
     "we are using a linear model for the prediction step""
  */
  //we got this after setting Bu = 0 and v = 0 (lesson 6). 
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**  
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred; // so it is the difference between predicted and measured 
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**  
    * update the state by using Extended Kalman Filter equations
  */
   //hx =  arctan(x) = atan(x_); 
   //test - try replacing z with x_
  /*float px = z(0);
  float py = z(1);
  float vx = z(2);
  float vy = z(3);*/

	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);

  // converting predicted data to polar co-ordinates
  float c1 = px*px+py*py;
  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1)); //sqrt(c1); (This is not working)
  float phi = atan2(x_(1), x_(0)); // atan2(py, px); this is also is not working
  float rho_dot;

  if (fabs(rho) < 0.0001) {
    rho_dot = 0;
  } else {
	  rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho; ; // (px*vx + py + vy) / rho;
  }
  
  VectorXd hx(3);
  hx << rho, phi, rho_dot;
  VectorXd y = z - hx;
  MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
