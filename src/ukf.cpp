#include <iostream>
#include "ukf.h"

UKF::UKF() {
	//TODO Auto-generated constructor stub
	Init();
}

UKF::~UKF() {
	//TODO Auto-generated destructor stub
}

void UKF::Init() {

}

/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/


void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

	//set state dimension
	int n_x = 5;

	//define spreading parameter
	double lambda = 3 - n_x;

	//set example state
	VectorXd x = VectorXd(n_x);
	x <<	5.7441,
				1.3800,
				2.2049,
				0.5015,
				0.3528;

	//set example covariance matrix
	MatrixXd P = MatrixXd(n_x, n_x);
	P <<	0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
				-0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
				0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
				-0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
				-0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

	//create sigma point matrix
	MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

	//calculate square root of P
	MatrixXd A = P.llt().matrixL();

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//your code goes here 

	//calculate sigma points ...
	Xsig.col(0) = x;
	for (int i = 0; i < n_x; i++)
	{
		Xsig.col(i + 1) = x + sqrt(lambda + n_x) * A.col(i);
		Xsig.col(i + 1 + n_x) = x - sqrt(lambda + n_x) * A.col(i);
	}

	//set sigma points as columns of matrix Xsig
	
	/*******************************************************************************
	* Student part end
	******************************************************************************/

	//print result
	//std::cout << "Xsig = " << std::endl << Xsig << std::endl;

	//write result
	*Xsig_out = Xsig;

	/* expected result:
	Xsig =
	5.7441  5.85768   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441
	1.38  1.34566  1.52806     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38
	2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049
	0.5015  0.44339 0.631886 0.516923 0.595227   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015
	0.3528 0.299973 0.462123 0.376339  0.48417 0.418721 0.405627 0.243477 0.329261  0.22143 0.286879
	*/

}

#include <iostream>
#include "ukf.h"

UKF::UKF() {
	//TODO Auto-generated constructor stub
	Init();
}

UKF::~UKF() {
	//TODO Auto-generated destructor stub
}

void UKF::Init() {

}


/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

	//set state dimension
	int n_x = 5;

	//set augmented dimension
	int n_aug = 7;

	//Process noise standard deviation longitudinal acceleration in m/s^2
	double std_a = 0.2;

	//Process noise standard deviation yaw acceleration in rad/s^2
	double std_yawdd = 0.2;

	//define spreading parameter
	double lambda = 3 - n_aug;

	//set example state
	VectorXd x = VectorXd(n_x);
	x << 5.7441,
			1.3800,
			2.2049,
			0.5015,
			0.3528;

	//create example covariance matrix
	MatrixXd P = MatrixXd(n_x, n_x);
	P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
		-0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
		0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
		-0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
		-0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//create augmented mean state
	x_aug << x, 0, 0;
	//create augmented covariance matrix
	P_aug <<	P, ArrayXXf::Zero(2, 2), 
						ArrayXXf::Zero(1, n_x),std_a, 0,
						ArrayXXf::Zero(1, n_x), 0, std_yawdd;
	//create square root matrix
	MatrixXd A_aug = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x;
	for (int i = 0; i < n_x; i++)
	{
		Xsig_aug.col(i + 1) = x + sqrt(lambda + n_x) * A_aug.col(i);
		Xsig_aug.col(i + 1 + n_x) = x - sqrt(lambda + n_x) * A_aug.col(i);
	}
		/*******************************************************************************
		* Student part end
		******************************************************************************/

		//print result
		std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

	//write result
	*Xsig_out = Xsig_aug;

	/* expected result:
	Xsig_aug =
	5.7441  5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441
	1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38
	2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049
	0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015
	0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528
	0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0
	0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641
	*/

}