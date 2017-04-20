#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

	///* initially set to false, set to true in first call of ProcessMeasurement
	bool is_initialized_;

	///* if this is false, laser measurements will be ignored (except for init)
	bool use_laser_;

	///* if this is false, radar measurements will be ignored (except for init)
	bool use_radar_;

	///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	VectorXd x_;

	///* state covariance matrix
	MatrixXd P_;

	///* sigma points matrix
	MatrixXd Xsig_;

	///* augmented state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	VectorXd x_aug_;

	///* augmented state covariance matrix
	MatrixXd P_aug_;

	///* augmented sigma points matrix
	MatrixXd Xsig_aug_;

	///* predicted sigma points matrix
	MatrixXd Xsig_pred_;

	///* sigma points in measurement space matrix - Radar
	MatrixXd Zsig_Radar_;

	///* radar measurement
	VectorXd z_Radar_;
	
	///* mean predicted radar measurement
	VectorXd z_pred_Radar_;

	///* radar measurement covariance matrix S
	MatrixXd S_Radar_;

	///* radar measurement noise covariance matrix
	MatrixXd R_Radar_;

	///* radar cross correlation matrix
	MatrixXd Tc_Radar_;

	///* sigma points in measurement space matrix - Radar
	MatrixXd Zsig_Lidar_;

	///* radar measurement
	VectorXd z_Lidar_;

	///* mean predicted radar measurement
	VectorXd z_pred_Lidar_;

	///* radar measurement covariance matrix S
	MatrixXd S_Lidar_;

	///* radar measurement noise covariance matrix
	MatrixXd R_Lidar_;

	///* radar cross correlation matrix
	MatrixXd Tc_Lidar_;

	///* time when the state is true, in us
	long long time_us_;

	///* delta time between two measurements
	double delta_t_;

	///* Process noise standard deviation longitudinal acceleration in m/s^2
	double std_a_;

	///* Process noise standard deviation yaw acceleration in rad/s^2
	double std_yawdd_;

	///* Laser measurement noise standard deviation position1 in m
	double std_laspx_;

	///* Laser measurement noise standard deviation position2 in m
	double std_laspy_;

	///* Radar measurement noise standard deviation radius in m
	double std_radr_;

	///* Radar measurement noise standard deviation angle in rad
	double std_radphi_;

	///* Radar measurement noise standard deviation radius change in m/s
	double std_radrd_;

	///* Weights of sigma points
	VectorXd weights_;

	///* State dimension
	int n_x_;

	///* Augmented state dimension
	int n_aug_;

	///* Measurement state dimension Radar
	int n_z_Radar_;

	///* Measurement state dimension Lidar
	int n_z_Lidar_;

	///* Sigma point spreading parameter
	double lambda_;

	///* the current NIS for radar
	double NIS_radar_;

	///* the current NIS for laser
	double NIS_laser_;

	/**
	* Constructor
	*/
	UKF();

	/**
	* Destructor
	*/
	virtual ~UKF();

	/**
	* ProcessMeasurement
	* @param meas_package The latest measurement data of either radar or laser
	*/
	void ProcessMeasurement(MeasurementPackage measurement_pack);

	/**
	* Prediction Predicts sigma points, the state, and the state covariance
	* matrix
	* @param delta_t Time between k and k+1 in s
	*/
	void Prediction(double delta_t);

	/**
	* Updates the state and the state covariance matrix using a laser measurement
	* @param meas_package The measurement at k+1
	*/
	void UpdateLidar(MeasurementPackage measurement_pack);

	/**
	* Updates the state and the state covariance matrix using a radar measurement
	* @param meas_package The measurement at k+1
	*/
	void UpdateRadar(MeasurementPackage measurement_pack);
};

#endif /* UKF_H */