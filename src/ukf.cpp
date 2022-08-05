#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;

  //start x_ values;
  x_ << 0.0,
        0.0,
        0.0,
        0.0,
        0.0;

  P_ <<   1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0,            
          0.0, 0.0, 0.0, 1.0, 0.0, 
          0.0, 0.0, 0.0, 0.0, 1.0;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  weights_ =  VectorXd(2*n_aug_+1);

  weights_(0) = lambda_ / (lambda_ + n_aug_) ;
  for(int i = 1; i < 2* n_aug_+1 ; i++){
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
  

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_){
    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_[0] = meas_package.raw_measurements_[0]; //px
      x_[1] = meas_package.raw_measurements_[1]; //py

      //P_(0,0) = std_laspx_ * std_laspx_;
      //P_(1,1) = std_laspy_ * std_laspy_;

    }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double r = meas_package.raw_measurements_[0]; //px
      double phi = meas_package.raw_measurements_[1]; //py
      x_[0] = r * cos(phi);
      x_[1] = r* sin(phi);
      /*
      P_(0,0) = std_radr_ * std_radr_;
      P_(1,1) = std_radr_ * std_radr_;
      P_(2,2) = std_radrd_ * std_radrd_;
      P_(3,3) = std_radphi_ * std_radphi_;
      P_(4,4) = std_radrd_ * std_radrd_;
      */

    }
    // capture the time step.
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  } 
  
  long long diff = meas_package.timestamp_ - time_us_;
  time_us_ = meas_package.timestamp_; // update time_us_ !!! Important update
  double delta_t = static_cast<double>(diff) / static_cast<double>(1e6);

  Prediction(delta_t);

  if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
    UpdateLidar(meas_package);
  if(use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  //PART1 -> Generating sigma points && Augmentation 
  VectorXd x_aug_ = VectorXd(n_aug_);
  MatrixXd P_aug_ = MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(n_x_,n_x_)= std_a_*std_a_;
  P_aug_(n_x_+1,n_x_+1)=std_yawdd_*std_yawdd_;

  MatrixXd sqrtP_aug = P_aug_.llt().matrixL();
  //cout << "sqrtP_aug" << endl;
  //cout << sqrtP_aug << endl;

  Xsig_aug_.fill(0);
  Xsig_aug_.col(0) = x_aug_;

  for(int i = 0; i < n_aug_ ;i++){
    Xsig_aug_.col(i+1) = x_aug_ + (sqrt(lambda_ + n_aug_) * sqrtP_aug.col(i));
    Xsig_aug_.col(n_aug_+i+1) = x_aug_ - (sqrt(lambda_ + n_aug_) * sqrtP_aug.col(i));
  }

  //cout << "Xsig_aug_" << endl;
  //cout << Xsig_aug_ << endl;
  //PART3 -> Predict Sigma Points

  for(int i = 0; i < 2*n_aug_+1 ;i++){
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);


    double px_p, py_p;
    if(fabs(yawd) > 0.001){
      px_p = p_x + ( (v/yawd) * (sin(yaw + (yawd * delta_t)) - sin(yaw)));
      py_p = p_y + ( (v/yawd) * ( cos(yaw) - cos(yaw + (yawd * delta_t))));
    }else{
      px_p = p_x + (v * delta_t * cos(yaw) );
      py_p = p_y + (v * delta_t * sin(yaw) );
    }

    // add noise factor
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    px_p += 0.5* nu_a * delta_t * delta_t * cos(yaw);
    py_p += 0.5* nu_a * delta_t * delta_t * sin(yaw);

    double v_p = v;
    double yaw_p = yaw + (yawd * delta_t); 
    double yawd_p = yawd;

    v_p +=  nu_a * delta_t;
    yaw_p += 0.5* nu_yawdd * delta_t * delta_t;
    yawd_p += nu_yawdd* delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

  }

  //cout << "Xsig_pred_" << endl;
  //cout << Xsig_pred_ << endl;

  //Predicted Mean (5x1)

  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1 ;i++){
    x_pred = x_pred + (weights_(i)*Xsig_pred_.col(i));  
  }
 
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  //Predicted Covariance (5x5)
  for(int i = 0; i < 2*n_aug_+1 ;i++){
    VectorXd x_diff =  Xsig_pred_.col(i) - x_pred; //5x1

    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
    P_pred = P_pred + (weights_(i)*x_diff * x_diff.transpose());  
  }

  //cout << "P_pred : " << P_pred << endl;

  x_ = x_pred;
  P_ = P_pred;
  //cout << "x_ : "<< endl;
  //cout << x_ << endl;

  //cout << "P_ : " << endl;
  //cout << P_ << endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    VectorXd z = meas_package.raw_measurements_ ; 

    int n_z = 2;
    MatrixXd z_sig = MatrixXd( n_z , 2*n_aug_+1);
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S_ = MatrixXd( n_z , n_z );

    z_sig.fill(0.0);
    for(int i = 0; i < 2*n_aug_+1 ;i++){
      z_sig(0,i) =  Xsig_pred_(0,i);
      z_sig(1,i) =  Xsig_pred_(1,i);
    }
    //Predicted Mean (2x1)
    z_pred.fill(0.0);
    for(int i = 0; i < 2*n_aug_+1 ;i++){
      z_pred = z_pred + (weights_(i)*z_sig.col(i));  
    }
  
    //Predicted Covariance (2x2)
    S_.fill(0.0);
    for(int i = 0; i < 2*n_aug_+1 ;i++){
      VectorXd z_diff =  z_sig.col(i) - z_pred; //2x1
      S_ = S_ + (weights_(i)*z_diff * z_diff.transpose());  
    }
    
    MatrixXd L_ = MatrixXd( n_z , n_z );
    L_ << std_laspx_ * std_laspx_ , 0 ,
          0, std_laspy_ * std_laspy_ ; 
    
    S_ = S_ + L_ ;
    
    //cross correlation matrix
    MatrixXd Tc = MatrixXd(n_x_, n_z); // 5x2
    Tc.fill(0.0);
    for(int i = 0; i < 2*n_aug_+1 ;i++){
      VectorXd z_diff =  z_sig.col(i) - z_pred; //2x1


      VectorXd x_diff =  Xsig_pred_.col(i) - x_; //5x1

      while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
      while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

      Tc = Tc + (weights_(i)* x_diff * z_diff.transpose());  
    }

    MatrixXd K = Tc * S_.inverse(); // 5x2 * 2x2 = 5x2

    VectorXd z_res = z - z_pred; // 2x1

    //Update state vactor x_
    x_ = x_ +  K * z_res; // 5x2 * 2x1 = 5x1 

    //Update covariance mat P_
    P_ = P_ - (K * S_ * K.transpose()); 
  
    double NIS = z_res.transpose() * S_.inverse() * z_res ; // 1x1 = (1x2 2x2 2x1)
    
    cout << "NIS for LIDAR: " << NIS << endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   // Incoming radar measurement
   VectorXd z = meas_package.raw_measurements_ ; 

  // Z  5x15 ->  3x15 
  //reuse sigma points
  int n_z = 3;
  MatrixXd z_sig = MatrixXd( n_z , 2*n_aug_+1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S_ = MatrixXd( n_z , n_z );

  z_sig.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1 ;i++){
     double p_x =  Xsig_pred_(0,i);
     double p_y =  Xsig_pred_(1,i);
     double v =  Xsig_pred_(2,i);
     double yaw =  Xsig_pred_(3,i);

     double v1 = v*cos(yaw);
     double v2 = v*sin(yaw);

     z_sig(0,i) = sqrt(p_x*p_x + p_y*p_y);  
     z_sig(1,i) = atan2(p_y, p_x);
     z_sig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y); 

  }
  //Predicted Mean (3x1)
  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1 ;i++){
    z_pred = z_pred + (weights_(i)*z_sig.col(i));  
  }
 
  //Predicted Covariance (3x3)
  S_.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1 ;i++){
    VectorXd z_diff =  z_sig.col(i) - z_pred; //3x1

    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
    S_ = S_ + (weights_(i)*z_diff * z_diff.transpose());  
  }
  
  MatrixXd R_ = MatrixXd( n_z , n_z );
  R_ << std_radr_ * std_radr_ , 0 , 0 ,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_ ; 

  S_ = S_ + R_ ;
  //cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z); // 5x3
  Tc.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1 ;i++){
    VectorXd z_diff =  z_sig.col(i) - z_pred; //3x1

    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;


    VectorXd x_diff =  Xsig_pred_.col(i) - x_; //5x1

    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + (weights_(i)* x_diff * z_diff.transpose());  
  }

  MatrixXd K = Tc * S_.inverse();

  VectorXd z_res = z - z_pred; // 3x1

  while(z_res(1) > M_PI) z_res(1) -= 2.*M_PI;
  while(z_res(1) < -M_PI) z_res(1) += 2.*M_PI;

  //Update state vactor x_
  x_ = x_ +  K * z_res;

  //Update covariance mat P_
  P_ = P_ - (K * S_ * K.transpose()); 

  double NIS = z_res.transpose() * S_.inverse() * z_res ; // 1x1 = (1x3 3x3 3x1)
  cout << "NIS for RADAR: " << NIS << endl;
}