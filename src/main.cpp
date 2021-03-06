#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "FusionEKF.h"
#include "ground_truth_package.h"
#include "measurement_package.h"
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {

    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;
	
    if (sensor_type.compare("L") == 0) {
      // LASER MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      double x;
			double y;
      iss >> x;
      iss >> y;
      meas_package.raw_measurements_ << x, y;

      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } 
		else if (sensor_type.compare("R") == 0) {
      // RADAR MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
			double ro;
			double phi;
			double ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;

      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
		double x_gt;
		double y_gt;
		double vx_gt;
		double vy_gt;
		double yaw_gt;
		double yaw_rate_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
		iss >> yaw_gt;
		iss >> yaw_rate_gt;

    gt_package.gt_values_ = VectorXd(6);
		gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt, yaw_gt, yaw_rate_gt;
    gt_pack_list.push_back(gt_package);
  }

  // used to compute the RMSE later
	// estimations and ground truth for EKF
  vector<VectorXd> estimations_EKF;
	vector<VectorXd> ground_truth_EKF; 

	// estimations and ground truth for UKF
	vector<VectorXd> estimations_UKF;
	vector<VectorXd> ground_truth_UKF;

	// Create a EKF instance
	FusionEKF fusionEKF;

	// Create a UKF instance
	UKF ukf;

	// Prepare for fusion
	size_t number_of_measurements = measurement_pack_list.size();

	// column names for output file

	// column names for output file

	out_file_ << "time_stamp" << "\t";
	out_file_ << "px_state" << "\t";
	out_file_ << "py_state" << "\t";
	out_file_ << "v_state" << "\t";
	out_file_ << "yaw_angle_state" << "\t";
	out_file_ << "yaw_rate_state" << "\t";
	out_file_ << "sensor_type" << "\t";
	out_file_ << "NIS" << "\t";
	out_file_ << "px_measured" << "\t";
	out_file_ << "py_measured" << "\t";
	out_file_ << "px_ground_truth" << "\t";
	out_file_ << "py_ground_truth" << "\t";
	out_file_ << "vx_ground_truth" << "\t";
	out_file_ << "vy_ground_truth" << "\t";
	out_file_ << "v_ground_truth" << "\t";
	out_file_ << "yaw_ground_truth" << "\t";
	out_file_ << "yaw_rate_ground_truth" << "\n";

	/*
	// Executing Extended Kalman Filter for each measurement
  for (size_t k = 0; k < number_of_measurements; ++k) {
    // start filtering from the second frame (the speed is unknown in the first
    // frame)
    fusionEKF.ProcessMeasurement(measurement_pack_list[k]);
		
    // output the estimation
    //out_file_ << fusionEKF.ekf_.x_(0) << "\t";
    //out_file_ << fusionEKF.ekf_.x_(1) << "\t";
    //out_file_ << fusionEKF.ekf_.x_(2) << "\t";
    //out_file_ << fusionEKF.ekf_.x_(3) << "\t";

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // output the estimation
      //out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
      //out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      //out_file_ << ro * cos(phi) << "\t"; // p1_meas
      //out_file_ << ro * sin(phi) << "\t"; // ps_meas
    }

    // output the ground truth packages
    //out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    //out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    //out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    //out_file_ << gt_pack_list[k].gt_values_(3) << "\n";

    estimations_EKF.push_back(fusionEKF.ekf_.x_);
    ground_truth_EKF.push_back(gt_pack_list[k].gt_values_);
  }
	*/

	// Executing Unscented Kalman Filter for each measurement
	for (size_t k = 1; k < number_of_measurements; ++k) {
		// Call the UKF-based fusion
		ukf.ProcessMeasurement(measurement_pack_list[k]);

		// timestamp
		out_file_ << measurement_pack_list[k].timestamp_ << "\t";
		out_file_ << ukf.x_(0) << "\t"; // pos1 - est
		out_file_ << ukf.x_(1) << "\t"; // pos2 - est
		out_file_ << ukf.x_(2) << "\t"; // vel_abs -est
		out_file_ << ukf.x_(3) << "\t"; // yaw_angle -est
		out_file_ << ukf.x_(4) << "\t"; // yaw_rate -est

		// output lidar and radar specific data
		if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
			// sensor type
			out_file_ << "lidar" << "\t";

			// NIS value
			out_file_ << ukf.NIS_laser_ << "\t";

			// output the lidar sensor measurement px and py
			out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
			out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";

		}
		else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
			// sensor type
			out_file_ << "radar" << "\t";

			// NIS value
			out_file_ << ukf.NIS_radar_ << "\t";

			// output radar measurement in cartesian coordinates
			double ro = measurement_pack_list[k].raw_measurements_(0);
			double phi = measurement_pack_list[k].raw_measurements_(1);
			out_file_ << ro * cos(phi) << "\t"; // px measurement
			out_file_ << ro * sin(phi) << "\t"; // py measurement
		}

		// output the ground truth
		out_file_ << gt_pack_list[k].gt_values_(0) << "\t"; // pos x - ground truth
		out_file_ << gt_pack_list[k].gt_values_(1) << "\t";	// pos y - ground truth
		out_file_ << gt_pack_list[k].gt_values_(2) << "\t";	// vel_x - ground truth
		out_file_ << gt_pack_list[k].gt_values_(3) << "\t"; // vel_y - ground truth
		out_file_ << sqrt(gt_pack_list[k].gt_values_(2) * gt_pack_list[k].gt_values_(2) + gt_pack_list[k].gt_values_(3) * gt_pack_list[k].gt_values_(3)) << "\t"; // vel_y - ground truth
		out_file_ << gt_pack_list[k].gt_values_(4) << "\t"; // yaw - ground truth
		out_file_ << gt_pack_list[k].gt_values_(5) << "\n"; // yaw_rate - ground truth

		// convert ukf x vector to cartesian to compare to ground truth
		VectorXd ukf_x_cartesian_ = VectorXd(6);	

		double x_estimate_ = ukf.x_(0);
		double y_estimate_ = ukf.x_(1);
		double vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
		double vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));
		double yaw_estimate = ukf.x_(3);
		double yaw_rate_estimate = ukf.x_(4);

		ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_, yaw_estimate, yaw_rate_estimate;

		estimations_UKF.push_back(ukf_x_cartesian_);
		ground_truth_UKF.push_back(gt_pack_list[k].gt_values_);

	}



  // compute the accuracy (RMSE)
  Tools tools;
  //cout << "Accuracy - RMSE (EKF):" << endl << tools.CalculateRMSE(estimations_EKF, ground_truth_EKF) << endl;
	cout << "Accuracy - RMSE (UKF):" << endl << tools.CalculateRMSE(estimations_UKF, ground_truth_UKF) << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

	cout << "Done!" << endl;
  return 0;
}
