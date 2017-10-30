#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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

/*
Prediction:
x' = Fx + u (11)
P' = FPF^T + Q (12)
*/
void KalmanFilter::Predict() {
	/**
		TODO:
		* predict the state
	*/
	x_ = F_ * x_;

	P_ = F_ * P_ * F_.transpose() + Q_;
}

/*
Measurement Update:
y = z − Hx'			(13)
S = HP'(H^T) + R		(14)
K = P'(H^T)(S^−1)	(15)
x = x' + Ky			(16)
P = (I − KH)P'		(17)
*/
void KalmanFilter::Update(const VectorXd &z) {
	/**
		TODO:
		* update the state by using Kalman Filter equations
	*/

	VectorXd y = z - H_ * x_;
	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();

	x_ = x_ + (K * y);

	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - K * H_) * P_;
}

/*
To calculate y, we use the equations that map the predicted location x'
from Cartesian coordinates to polar coordinates
		| sqrt(px'^2+py'^2)					|
h(x') = | arctan(py'/px')					|
		| (px'vx'+py'vy')/sqrt(px'^2+py'^2)	|
*/
void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
		TODO:
		* update the state by using Extended Kalman Filter equations
	*/
	float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
	float phi = atan2(x_(1), x_(0));
	float rho_dot;
	if (fabs(rho) < 0.0001) {
		// make ro_dot 0.0 if too small
		rho_dot = 0;
	} else {
		rho_dot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
	}

	VectorXd Hx(3);
	Hx << rho, phi, rho_dot;
	
	VectorXd y = z - Hx;

	while (y(1)>M_PI)
	{
		y(1) -= 2 * M_PI;
	}
	while (y(1)<-M_PI)
	{
		y(1) += 2 * M_PI;
	}

	MatrixXd H_trans = H_.transpose();

	MatrixXd S = H_ * P_ * H_trans + R_;
	MatrixXd K = P_ * H_trans * S.inverse();

	x_ = x_ + (K * y);

	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - K * H_) * P_;
}