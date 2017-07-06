#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd vec(4);
	vec << 0,0,0,0;

	if(estimations.size() != ground_truth.size() || estimations.size() == 0)//check if size correct
	{
		cout << "Invalid estimation or ground_truth data" << endl;
		return vec;
	}

	for(unsigned int i = 0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		vec += residual;
	}

	//calculate the mean
	vec = vec / estimations.size();

	vec = vec.array().sqrt();
	return vec;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  	/**
  	TODO:
    	* Calculate a Jacobian here.
  	*/
  	MatrixXd  mtx(3,4);

  	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	if(fabs(c1) < 1e-4)//exception for overflow
	{
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return mtx;
	}

	mtx << (px/c2), (py/c2), 0, 0,//line 1
			-(py/c1), (px/c1), 0, 0,//line 2
			py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;//line 3

	return mtx;
}
