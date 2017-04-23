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
	// initialize higher uncertainty for yaw and yaw_dot
	//P_(2, 2) = 10;
	//P_(3, 3) = 10;


	// initialize sigma points matrix
	Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

	// initialize augmented state dimension [pos1 pos2 vel_abs yaw_angle yaw_rate nu_vel vu_yaw_rate]
	n_aug_ = 7;

	// initialize augmented state vector
	x_aug_ = VectorXd(n_aug_);
	x_aug_ = VectorXd::Ones(n_aug_);

	// initialize augmented covariance matrix
	P_aug_ = MatrixXd(n_aug_, n_aug_);
	P_aug_ = MatrixXd::Identity(n_aug_, n_aug_);

	// initialize augmented sigma points matrix
	Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	
	// initialize prediction sigma points matrix
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

	// initialize measurement space dimension Radar(r, phi, and r_dot)
	n_z_Radar_ = 3;

	//ininitalize sigma points in measurement space matrix - Radar
	Zsig_Radar_ = MatrixXd(n_z_Radar_, 2 * n_aug_ + 1);

	//initialize mean predicted measurement - Radar
	z_pred_Radar_ = VectorXd(n_z_Radar_);

	//initialize measurement - Radar
	z_Radar_ = VectorXd(n_z_Radar_);

	// initialize measurement covariance matrix S - Radar
	S_Radar_ = MatrixXd(n_z_Radar_, n_z_Radar_);

	// initialize measurement noise covariance matrix - Radar
	R_Radar_ = MatrixXd(n_z_Radar_, n_z_Radar_);

	// initialize cross correlation matrix - Radar
	Tc_Radar_ = MatrixXd(n_x_, n_z_Radar_);

	// initialize measurement space dimension Lidar (x, y)
	n_z_Lidar_ = 2;

	//ininitalize sigma points in measurement space matrix - Radar
	Zsig_Lidar_ = MatrixXd(n_z_Lidar_, 2 * n_aug_ + 1);

	//initialize mean predicted measurement - Radar
	z_pred_Lidar_ = VectorXd(n_z_Lidar_);

	//initialize measurement - Radar
	z_Lidar_ = VectorXd(n_z_Lidar_);

	// initialize measurement covariance matrix S - Radar
	S_Lidar_ = MatrixXd(n_z_Lidar_, n_z_Lidar_);

	// initialize measurement noise covariance matrix - Radar
	R_Lidar_ = MatrixXd(n_z_Lidar_, n_z_Lidar_);
	
	// initialize cross correlation matrix - Radar
	Tc_Lidar_ = MatrixXd(n_x_, n_z_Lidar_);

	// initialize spreading parameter lambda
	lambda_ = 3 - n_aug_;

	//create vector for weights
	weights_ = VectorXd(2 * n_aug_ + 1);

	// TODO: Tune longitudinal acceleration process noise
	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.2;

	// TODO: Tune yaw acceleration process noise
	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.6; //M_PI/3;

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

	R_Radar_ << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;

	R_Lidar_ << std_laspx_ * std_laspx_, 0,
		0, std_laspy_ * std_laspy_;

}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
	if (!is_initialized_) {
		// first measurement
		// TODO: Tune state if needed
		cout << "UKF: " << endl;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			// Read measurements to rho, phi and rho_dot respectively
			double rho = measurement_pack.raw_measurements_[0];			// range
			double phi = measurement_pack.raw_measurements_[1];			// bearing
			double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho

			// Convert from polar to cartesian coordinates
			double x = rho * cos(phi);
			double y = rho * sin(phi);
			//double v = sqrt((rho_dot * cos(phi) * rho_dot * cos(phi)) + (rho_dot * sin(phi) * rho_dot * sin(phi)));
			double v = 0;
			double sigma = 0;
			double sigma_dot = 0;
			x_ << x, y, v, sigma, sigma_dot;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			double x = measurement_pack.raw_measurements_[0];
			double y = measurement_pack.raw_measurements_[1];
			double v = 0;
			double sigma = 0;
			double sigma_dot = 0;
			x_ << x, y, v, sigma, sigma_dot;
		}
		// Deal with the special case initialisation problems
		if (fabs(x_(0)) < EPS) {
			return;
			//x_(0) = 1;
			// higher uncertainty for p_x
			// P_(0, 0) = 100;
		}
		if (fabs(x_(1)) < EPS) {
			return;
			// x_(1) = 1;
			// higher uncertainty for p_y
			// P_(1, 1) = 100;
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
		delta_t_ = measurement_pack.timestamp_ - time_us_;
		delta_t_ /= 1000000.0; // convert delta t to s
		time_us_ = measurement_pack.timestamp_;
		//cout << "delta t: " << delta_t_ << endl;
		/*
		if (delta_t_ == 0) {
			delta_t_ = 1;
			//cout << "delta_t_ is zero" << endl;
		}
		*/

		// Predict
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

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {
	// 1. Generate Sigma Points for timestep k
	//create augmented mean state
	x_aug_.head(5) = x_;
	x_aug_(5) = 0;
	x_aug_(6) = 0;

	//create augmented covariance matrix
	P_aug_.fill(0.0);
	P_aug_.topLeftCorner(5, 5) = P_;
	P_aug_(5, 5) = std_a_*std_a_;
	P_aug_(6, 6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug_.llt().matrixL();

	//create augmented sigma points
	Xsig_aug_.col(0) = x_aug_;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	// 2. Predict Sigma Points for timestep k+1
	for (int i = 0; i< 2 * n_aug_ + 1; i++) {
		//extract values for better readability
		double p_x = Xsig_aug_(0, i);
		double p_y = Xsig_aug_(1, i);
		double v = Xsig_aug_(2, i);
		double yaw = Xsig_aug_(3, i);
		double yawd = Xsig_aug_(4, i);
		double nu_a = Xsig_aug_(5, i);
		double nu_yawdd = Xsig_aug_(6, i);

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


	// 3. Predict Mean and Covariance
	// set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i<2 * n_aug_ + 1; i++) {
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
	
	/*
	cout << "x_ (prediction): " << endl;
	cout << x_ << endl;
	*/
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
	// 1. Predict Measurement
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		// extract values for better readibility
		Zsig_Lidar_(0, i) = Xsig_pred_(0, i); // px
		Zsig_Lidar_(1, i) = Xsig_pred_(1, i); // py
	}
	
	// calculate measurement mean
	z_pred_Lidar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred_Lidar_ = z_pred_Lidar_ + weights_(i) * Zsig_Lidar_.col(i);
	}

	// calculate measurement covariance and cross-covariance matrix
	S_Lidar_.fill(0.0);
	Tc_Lidar_.fill(0.0);
	
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		// measurement difference
		VectorXd z_diff = Zsig_Lidar_.col(i) - z_pred_Lidar_;

		// measurement covariance
		S_Lidar_ = S_Lidar_ + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		// angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		// cross-correlation
		Tc_Lidar_ = Tc_Lidar_ + weights_(i) * x_diff * z_diff.transpose();
	}

	// Adjust measurement covariance matrix with measurement noise
	S_Lidar_ = S_Lidar_ + R_Lidar_;
	
	// 2. Update
	// Read measurements to p_x and p_y respectively
	double p_x = EPS;
	double p_y = EPS;
	if (fabs(measurement_pack.raw_measurements_[0]) > EPS) {
		p_x = measurement_pack.raw_measurements_[0];
	}
	if (fabs(measurement_pack.raw_measurements_[1]) > EPS) {
		p_y = measurement_pack.raw_measurements_[1];			
	}
	z_Lidar_ << p_x, p_y;

	//Kalman gain K;
	MatrixXd K_Lidar_ = Tc_Lidar_ * S_Lidar_.inverse();
	
	//residual
	VectorXd z_diff = z_Lidar_ - z_pred_Lidar_;
	
	//update state mean and covariance matrix	
	x_ = x_ + K_Lidar_ * z_diff;
	P_ = P_ - K_Lidar_ * S_Lidar_ * K_Lidar_.transpose();
	/*
	cout << "x_ (Lidar update): " << endl;
	cout << x_ << endl;
	*/

	// 3. Calculate NIS
	NIS_laser_ = z_diff.transpose() * S_Lidar_.inverse() * z_diff;
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
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
		Zsig_Radar_(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig_Radar_(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig_Radar_(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	// calculate measurement mean
	z_pred_Radar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred_Radar_ = z_pred_Radar_ + weights_(i) * Zsig_Radar_.col(i);
	}
	
	// calculate measurement covariance and cross-covariance matrix
	S_Radar_.fill(0.0);
	Tc_Radar_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
		// measurement difference
		VectorXd z_diff = Zsig_Radar_.col(i) - z_pred_Radar_;

		// angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
		
		// measurement covariance
		S_Radar_ = S_Radar_ + weights_(i) * z_diff * z_diff.transpose();

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		// angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		// cross-correlation
		Tc_Radar_ = Tc_Radar_ + weights_(i) * x_diff * z_diff.transpose();
	}
	
	// Adjust measurement covariance matrix with measurement noise
	S_Radar_ = S_Radar_ + R_Radar_;
	
	// 2. Update
	// Read measurements to rho, phi and rho_dot respectively
	double rho = measurement_pack.raw_measurements_[0];			// range
	double phi = measurement_pack.raw_measurements_[1];			// bearing
	double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho

	z_Radar_ << rho, phi, rho_dot;

	//Kalman gain K;
	MatrixXd K_Radar_ = Tc_Radar_ * S_Radar_.inverse();

	//residual
	VectorXd z_diff = z_Radar_ - z_pred_Radar_;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K_Radar_ * z_diff;
	P_ = P_ - K_Radar_ * S_Radar_ * K_Radar_.transpose();

	/*
	cout << "x_ (Radar update): " << endl;
	cout << x_ << endl;	
	*/

	// 3. Calcualte NIS
	NIS_radar_ = z_diff.transpose() * S_Radar_.inverse() * z_diff;
}