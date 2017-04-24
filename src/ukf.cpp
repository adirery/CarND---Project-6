#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.0001 // A very small number

/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() {
	// set initialization flag
	is_initialized_ = false;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initialize timestamp
	time_us_ = 0;

	// initialize delta_t
	delta_t_ = 0;

	//initialize state dimension ([pos1 pos2 vel_abs yaw_angle yaw_rate])
	n_x_ = 5;

	// initialize state vector
	x_ = VectorXd(n_x_);
	x_ = VectorXd::Ones(n_x_);

	// initialize covariance matrix
	P_ = MatrixXd(n_x_, n_x_);
	P_ = MatrixXd::Identity(n_x_, n_x_);

	// initialize sigma points matrix
	Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

	// initialize augmented state dimension [pos1 pos2 vel_abs yaw_angle yaw_rate nu_vel vu_yaw_rate]
	n_aug_ = 7;

	// initialize prediction sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	// initialize spreading parameter lambda
	lambda_ = 3 - n_aug_;

	//create vector for weights
	weights_ = VectorXd(2 * n_aug_ + 1);

	// TODO: Tune longitudinal acceleration process noise
	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3.0;

	// TODO: Tune yaw acceleration process noise
	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.8; //M_PI/3;

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
}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} measurement_pack The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
	if (!is_initialized_) {
		// first measurement
		// TODO: Tune state if needed
		cout << "UKF: " << endl;

		// set previous timestamp to current value
		time_us_ = measurement_pack.timestamp_;

		// set state for first measurement based on measurement type
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

			// convert from polar to cartesian coordinates
			double x_cart = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
			double y_cart = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
			if (x_cart == 0 || y_cart == 0) {
				return;
			}
			x_ << x_cart, y_cart, 0, 0, 0;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			if (measurement_pack.raw_measurements_[0] == 0 || measurement_pack.raw_measurements_[1] == 0) {
				return;
			}
			x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
		}

		// Print the initialization results
		cout << "UKF init state x_: " << endl;
		cout << x_ << endl;
		cout << "UKF init cov. P_: " << endl;
		cout << P_ << endl;

		// Save the initiall timestamp for dt calculation
		time_us_ = measurement_pack.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}
	else if (is_initialized_) {
		// Read timestamp and calculate delta_t for Prediction step
		delta_t_ = (measurement_pack.timestamp_ - time_us_) / 1000000.0; // convert delta t to s
		time_us_ = measurement_pack.timestamp_;
		
		// Predict
		while (delta_t_ > 0.1)
		{
			const double dt = 0.05;
			Prediction(dt);
			delta_t_ -= dt;
		}	
		Prediction(delta_t_);
		// Update
		if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			UpdateLidar(measurement_pack);
		}
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			UpdateRadar(measurement_pack);
		}
	}
}

void UKF::GenerateSigmaPoints(double delta_t) {
	// 1. Generate Sigma Points for timestep k
	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	// 2. Predict Sigma Points for timestep k+1
	for (int i = 0; i< 2 * n_aug_ + 1; i++) {
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > EPS) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
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
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
	// 1. Generate Sigma Points
	GenerateSigmaPoints(delta_t);

	// 2. Predict Mean and Covariance
	// set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i< 2 * n_aug_ + 1; i++) {
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	
	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} measurement_pack
*/
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
	// laser will have variables px, py
	int n_z = 2;

	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// 1. Predict Measurement
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		Zsig(0, i) = p_x;
		Zsig(1, i) = p_y;
	}
	
	// calculate measurement mean
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	// calculate measurement covariance and cross-covariance matrix
	// measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
																							// residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	// laser measurement covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_ * std_laspx_, 0.0,
		0.0, std_laspy_ * std_laspy_;

	S = S + R;

	// create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	// calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		// normalize angles
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

	}

	// Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// residual
	VectorXd z = VectorXd(2);
	z <<	measurement_pack.raw_measurements_[0],
				measurement_pack.raw_measurements_[1];

	VectorXd z_diff = z - z_pred;

	// update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	// calculate NIS
	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} measurement_pack
*/
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
	// laser will have variables r, phi, r_dot
	int n_z = 3;

	// create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	// 1. Predict Measurement
	// transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                       //r
		Zsig(1, i) = atan2(p_y, p_x);                               //phi
		Zsig(2, i) = ((p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y));   //r_dot
	}

	// calculate measurement mean
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	// calculate measurement covariance and cross-covariance matrix
	// measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	// laser measurement covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R <<	std_radr_ * std_radr_, 0.0, 0.0,
				0.0, std_radphi_ * std_radphi_, 0.0,
				0.0, 0.0, std_radrd_ * std_radrd_;

	S = S + R;

	// create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	// calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		// normalize angles
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}
	// Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// residual
	VectorXd z = VectorXd(2);
	z <<	measurement_pack.raw_measurements_[0],
				measurement_pack.raw_measurements_[1],
				measurement_pack.raw_measurements_[2];
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;


	// update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	// calculate NIS
	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}