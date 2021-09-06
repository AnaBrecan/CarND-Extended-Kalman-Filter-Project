#include "kalman_filter.h"



using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * TODO: predict the state
   */
   x_ = F_*x_;
   P_ = F_*P_*(F_.transpose()) + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y;
  MatrixXd S, K, I;
  I = MatrixXd::Identity(P_.rows(), P_.cols());
  
  y = z - H_*x_;
  S = H_*P_*(H_.transpose()) + R_;
  K = P_*(H_.transpose())*(S.inverse());
    
  x_ = x_ + K*y;
  P_ = (I - K*H_)*P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  VectorXd y(3), x_polar(3);
  MatrixXd S, K, I; //, Hj;
  I = MatrixXd::Identity(P_.rows(), P_.cols());
  //Hj = tools.CalculateJacobian(H_);
  
  // compute h(x_)
  // recover state parameters, i.e. position and velocity in the x and y direction
  auto px = x_(0);
  auto py = x_(1);
  auto vx = x_(2);
  auto vy = x_(3);
  
  auto c = px*px + py*py;
 
  try{
    // check division by zero
    if (fabs(c) < 0.0001) {
      throw std::runtime_error("Math error: Attempted to divide by Zero\n");
    }
    auto x_1 = sqrt(c);
    auto x_2 = atan(py/px);
    auto x_3 = (px*vx + py*vy)/x_1;
    
    while (x_2 < -M_PI){
      x_2 += 2.0*M_PI;
    }
    
    while (x_2 > M_PI){ 
      x_2 -= 2.0*M_PI; 
    }
    
    x_polar << x_1, x_2, x_3;
  }
  catch (std::runtime_error& e){
    std::cout << "Exception occurred" << std::endl<< e.what();
    }
  
  y = z - x_polar;
  S = H_*P_*(H_.transpose()) + R_;
  K = P_*(H_.transpose())*(S.inverse());
    
  x_ = x_ + K*y;
  P_ = (I - K*H_)*P_;
  
}
