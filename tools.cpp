#include "tools.h"
#include <iostream>
#include <stdexcept>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  //check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }
  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {
      for (unsigned int j=0; j < estimations[i].size(); ++j) {
        rmse[j] = rmse[j] + (estimations[i][j] - ground_truth[i][j]) * (estimations[i][j] - ground_truth[i][j]);
      }
    }
    for (unsigned int j=0; j < estimations[0].size(); ++j) {
      rmse[j] = rmse[j] / estimations.size();
      rmse[j] = sqrt(rmse[j]);
    }
  // return the result
  return rmse;
  
  }

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
   MatrixXd Hj(3,4);
  // recover state parameters, i.e. position and velocity in the x and y direction
  auto px = x_state(0);
  auto py = x_state(1);
  auto vx = x_state(2);
  auto vy = x_state(3);
  
  //compute terms which appear in the expression of the Jacobi matrix
  auto c1 = px*px+py*py;
  auto c2 = sqrt(c1);
  auto c3 = (c1*c2);
 
  // compute the Jacobian matrix
  try{
    // check division by zero
    if (fabs(c1) < 0.0001) {
      throw std::runtime_error("Math error: Attempted to divide by Zero\n");
    }
    Hj << (px/c2), (py/c2), 0, 0,
       -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  }
  catch (std::runtime_error& e){
    cout << "Exception occurred" << std::endl<< e.what();
    }
  
  return Hj;
}
