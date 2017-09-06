#include "PID.h"

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {}

PID::~PID() {}

void PID::Init(double Kp, double Ki, double Kd) {
	_Kp = Kp;
	_Ki = Ki;
	_Kd = Kd;
	
	_p_error = 0.;
  	_i_error = 0.;
    _d_error = 0.;
}

void PID::UpdateError(double cte) {
	_p_error = cte;
	_i_error += cte;
	_d_error = cte - _pre_cte;
	_pre_cte = cte;
}

double PID::TotalError() {
	return -(_Kp * _p_error + _Ki * _i_error + _Kd * _d_error);
}

