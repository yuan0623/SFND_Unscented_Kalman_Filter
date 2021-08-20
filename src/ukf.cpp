#include "ukf.h"
#include "Eigen/Dense"
#include<iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;
  is_initialized_ = false;
  // initial state vector
  x_ = VectorXd(5);

  //x_<<10,10,0,0,0.3;
  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
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
    n_x_ = x_.size();
    n_aug_ = n_x_+2;

    lambda_ = 3-n_aug_;
    Xsig_pred_ = Eigen::MatrixXd(n_x_,2*n_aug_+1);

    weights_=Eigen::VectorXd(2*n_aug_+1);
    weights_(0) = lambda_/(lambda_+n_aug_);
    for (int i=1;i<2*n_aug_+1;i++)
    {
        weights_(i) = 1/(2*(lambda_+n_aug_));
    }


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    if(meas_package.sensor_type_==MeasurementPackage::RADAR and use_radar_==true)
    {
        if (is_initialized_ == false)
        {
            is_initialized_ = true;


            P_ << std_radr_* std_radr_, 0, 0, 0, 0,
                    0, std_radr_ * std_radr_, 0, 0, 0,
                    0, 0, std_radrd_ * std_radrd_, 0, 0,
                    0, 0, 0, std_radphi_ * std_radphi_, 0,
                    0, 0, 0, 0, std_radphi_ * std_radphi_;
        }
        else
            {
                UpdateRadar(meas_package);
        }

    }
    else if(meas_package.sensor_type_==MeasurementPackage::LASER and use_laser_==true)
    {
        if (is_initialized_ == false)
        {
            is_initialized_ = true;
            x_(0) = meas_package.raw_measurements_[0];
            x_(1) = meas_package.raw_measurements_[1];
            P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
                    0, std_laspy_ * std_laspy_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1;
        }
        else
        {
            UpdateLidar(meas_package);
        }

    }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

    predictSigmaPoints(Xsig_pred_,delta_t);

    Eigen::VectorXd x_predict_mean = VectorXd(5);
    for(int i=0;i<2*n_aug_+1;i++)
    {
        x_predict_mean = x_predict_mean+weights_(i)*(Xsig_pred_.col(i));
    }

    x_ = x_predict_mean;

    Eigen::MatrixXd P = MatrixXd(5, 5);

    for (int i=0; i<2*n_aug_+1;i++)
    {
        P = P + weights_(i)*((Xsig_pred_.col(i)-x_)*((Xsig_pred_.col(i)-x_).transpose()));
    }

    P_ = P;

    //std::cout<<"Xsig_pred_"<<std::endl;
    //std::cout<<Xsig_pred_<<std::endl;





}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    int n_z = 2;
    Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z,2*n_aug_+1);
    Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
    Eigen::MatrixXd S = Eigen::MatrixXd(n_z,n_z);

    for(int i=0; i<2*n_aug_+1; i++)
    {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        double d_yaw = Xsig_pred_(4,i);


        Zsig.col(i)<<p_x,
                     p_y;
    }

    for (int i = 0;i<2*n_aug_+1;i++)
    {
        z_pred = z_pred+weights_(i)*Zsig.col(i);
    }
    for (int i = 0;i<2*n_aug_+1;i++)
    {
        S = weights_(i)*((Zsig.col(i)-z_pred)*((Zsig.col(i)-z_pred).transpose()));
    }

    Eigen::MatrixXd R = Eigen::MatrixXd(n_z,n_z);
    R(0,0) = pow(std_laspx_,2);
    R(1,1) = pow(std_laspx_,2);

    S = S+R;


    Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_,n_z);
    for(int i=0;i<2*n_aug_+1;i++)
    {
        Tc = Tc + weights_(i)*((Xsig_pred_.col(i)-x_)*((Zsig.col(i)-z_pred).transpose()));
    }

    Eigen::MatrixXd K = Tc*S.inverse();

    x_ = x_ +K*(meas_package.raw_measurements_-z_pred);
    P_ = P_-K*S*K.transpose();



}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

    int n_z = 3;
    Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z,2*n_aug_+1);
    Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
    Eigen::MatrixXd S = Eigen::MatrixXd(n_z,n_z);

    for(int i=0; i<2*n_aug_+1; i++)
    {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        double d_yaw = Xsig_pred_(4,i);


        Zsig.col(i)<<sqrt(pow(p_x,2)+pow(p_y,2)),
                atan(p_y/p_x),
                (p_x*cos(yaw)*v+p_y*sin(yaw)*v)/(sqrt(pow(p_x,2)+pow(p_y,2)));
    }

    for (int i = 0;i<2*n_aug_+1;i++)
    {
        z_pred = z_pred+weights_(i)*Zsig.col(i);
    }
    for (int i = 0;i<2*n_aug_+1;i++)
    {
        S = weights_(i)*((Zsig.col(i)-z_pred)*((Zsig.col(i)-z_pred).transpose()));
    }

    Eigen::MatrixXd R = Eigen::MatrixXd(n_z,n_z);
    R(0,0) = pow(std_radr_,2);
    R(1,1) = pow(std_radphi_,2);
    R(2,2) = pow(std_radrd_,2);
    S = S+R;


    // update state
    Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_,n_z);
    for(int i=0;i<2*n_aug_+1;i++)
    {
        Tc = Tc + weights_(i)*((Xsig_pred_.col(i)-x_)*((Zsig.col(i)-z_pred).transpose()));
    }

    Eigen::MatrixXd K = Tc*S.inverse();

    x_ = x_ +K*(meas_package.raw_measurements_-z_pred);
    P_ = P_-K*S*K.transpose();
    //std::cout<<Zsig<<std::endl;
    //std::cout<<P_<<std::endl;
}
void UKF::generateSigmaPoints(Eigen::MatrixXd &Xsig_aug){
    Eigen::VectorXd x_aug;
    Eigen::MatrixXd P_aug;
    x_aug = VectorXd(n_aug_);
    P_aug = MatrixXd(n_aug_, n_aug_);

    x_aug.block(0,0,n_x_,1) = x_;
    x_aug(n_x_,0) = 0;
    x_aug(n_x_+1,0) = 0;



    P_aug.block(0,0,n_x_,n_x_) = P_;
    P_aug(n_x_,n_x_) = pow(std_a_,2);
    P_aug(n_x_+1,n_x_+1) = pow(std_yawdd_,2);


    Eigen::MatrixXd A = P_aug.llt().matrixL();


    Xsig_aug.col(0) = x_aug;
    for (int i=0; i<n_aug_;i++)
    {
        Xsig_aug.col(i+1) = x_aug+(sqrt(lambda_+n_aug_)*A).col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug-(sqrt(lambda_+n_aug_)*A).col(i);
    }
}
void UKF::predictSigmaPoints(Eigen::MatrixXd &Xsig_pred,double delta_t){
    Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(n_aug_,2*n_aug_+1);

    generateSigmaPoints(Xsig_aug);
    for(int i =0;i<2*n_aug_+1;i++){
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double d_yaw = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_ddyaw = Xsig_aug(6,i);

        double p_x_p, p_y_p;


        if (fabs(d_yaw)<0.001)
        {
            p_x_p = p_x+v/d_yaw*(sin(yaw+d_yaw*delta_t)-sin(yaw));
            p_y_p = p_y+v/d_yaw*(-cos(yaw+d_yaw*delta_t)+cos(yaw));
        }
        else
        {
            p_x_p = p_x + v*delta_t*cos(yaw);
            p_y_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw+d_yaw*delta_t;
        double d_yaw_p = d_yaw;

        // add noise
        p_x_p = p_x_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
        p_y_p = p_y_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
        v_p = v_p + nu_a*delta_t;
        yaw_p = yaw_p + 0.5*nu_ddyaw*delta_t*delta_t;
        d_yaw_p = d_yaw_p + nu_ddyaw*delta_t;

        Xsig_pred(0,i) = p_x_p;
        Xsig_pred(1,i) = p_y_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = d_yaw_p;


    }



}