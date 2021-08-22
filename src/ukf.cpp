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
  use_radar_ = true;
  is_initialized_ = false;
  // initial state vector
  x_ = VectorXd(5);

  //x_<<10,10,0,0,0.3;
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
    n_x_ = x_.size();
    n_aug_ = n_x_+2;

    lambda_ = 3-n_aug_;
    Xsig_pred_ = Eigen::MatrixXd::Zero(n_x_,2*n_aug_+1);

    weights_=Eigen::VectorXd::Zero(2*n_aug_+1);
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
    if(is_initialized_ == false)
    {
        std::cout<<"Initialization"<<std::endl;

        previous_time_ = meas_package.timestamp_;
        if(meas_package.sensor_type_==MeasurementPackage::RADAR and use_radar_==true)
        {
            auto rho = static_cast<float>(meas_package.raw_measurements_(0));
            auto phi = static_cast<float>(meas_package.raw_measurements_(1));
            // auto rho_dot = static_cast<float>(meas_package.raw_measurements_(2));
            x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
            std::cout<<x_<<std::endl;
            P_ << std_radr_* std_radr_, 0, 0, 0, 0,
                    0, std_radr_ * std_radr_, 0, 0, 0,
                    0, 0, std_radrd_ * std_radrd_, 0, 0,
                    0, 0, 0, std_radphi_ * std_radphi_, 0,
                    0, 0, 0, 0, std_radphi_ * std_radphi_;
            //predictSigmaPoints(Xsig_pred_,0.1);
            is_initialized_ = true;
        }
        else if(meas_package.sensor_type_==MeasurementPackage::LASER and use_laser_==true)
        {
            x_(0) = meas_package.raw_measurements_[0];
            x_(1) = meas_package.raw_measurements_[1];
            P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
                    0, std_laspy_ * std_laspy_, 0, 0, 0,
                    0, 0, 1, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1;
            //predictSigmaPoints(Xsig_pred_,0.1);
            is_initialized_ = true;
        }

    }
    else{
        if(meas_package.sensor_type_==MeasurementPackage::RADAR and use_radar_==true)
        {
            std::cout<<"ha radar"<<std::endl;
            UpdateRadar(meas_package);
        }
        else if(meas_package.sensor_type_==MeasurementPackage::LASER and use_laser_==true)
        {

            UpdateLidar(meas_package);
        }
    }

    auto dt = static_cast<double>((meas_package.timestamp_- previous_time_) * 1e-6); // or divided by 1000000.0
    previous_time_ = meas_package.timestamp_;

    if(dt>0){
        std::cout<<"prediction"<<std::endl;
        Prediction(dt);
    }


}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

    predictSigmaPoints(Xsig_pred_,delta_t);

    Eigen::VectorXd x_predict_mean = VectorXd::Zero(5);
    for(int i=0;i<2*n_aug_+1;i++)
    {
        x_predict_mean = x_predict_mean+weights_(i)*(Xsig_pred_.col(i));
    }

    x_ = x_predict_mean;

    Eigen::MatrixXd P = MatrixXd::Zero(5, 5);

    for (int i=0; i<2*n_aug_+1;i++)
    {
        P = P + weights_(i)*((Xsig_pred_.col(i)-x_)*((Xsig_pred_.col(i)-x_).transpose()));
    }

    P_ = P;






}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

    n_z_ = 2;
    Eigen::MatrixXd Zsig = Eigen::MatrixXd::Zero(n_z_,2*n_aug_+1);
    Eigen::VectorXd z_pred = Eigen::VectorXd::Zero(n_z_);
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_z_,n_z_);
    std::cout<<"Xsig_pred_"<<std::endl;
    std::cout<<Xsig_pred_<<std::endl;

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
    z_pred.fill(0);
    for (int i = 0;i<2*n_aug_+1;i++)
    {
        z_pred = z_pred+weights_(i)*Zsig.col(i);
    }

    for (int i = 0;i<2*n_aug_+1;i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = weights_(i)*z_diff*z_diff.transpose();
    }

    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(n_z_,n_z_);
    R(0,0) = pow(std_laspx_,2);
    R(1,1) = pow(std_laspx_,2);

    S = S+R;
    measurementUpdate(meas_package.raw_measurements_, Zsig,  z_pred,  S);

    /*
    Eigen::MatrixXd Tc = Eigen::MatrixXd::Zero(n_z_,n_z_);
    for(int i=0;i<2*n_aug_+1;i++)
    {

        Tc = Tc + weights_(i)*((Xsig_pred_.col(i)-x_)*((Zsig.col(i)-z_pred).transpose()));
    }

    Eigen::MatrixXd K = Tc*S.inverse();



    x_ = x_ +K*(meas_package.raw_measurements_-z_pred);

    P_ = P_-K*S*K.transpose();
    */


}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */


    n_z_ = 3;
    Eigen::MatrixXd Zsig = Eigen::MatrixXd::Zero(n_z_,2*n_aug_+1);
    Eigen::VectorXd z_pred = Eigen::VectorXd::Zero(n_z_);
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_z_,n_z_);

    for(int i=0; i<2*n_aug_+1; i++)
    {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        double d_yaw = Xsig_pred_(4,i);

        double v_x = cos(yaw) * v;
        double v_y = sin(yaw) * v;

        Zsig(0,i)=sqrt(pow(p_x,2)+pow(p_y,2));
        Zsig(1,i)=atan2(p_y,p_x);
        if (Zsig(0, i) < 0.001) {
            Zsig(2, i) = (p_x * v_x + p_y * v_y) / 0.001;
        } else {
            Zsig(2, i) = (p_x * v_x + p_y * v_y) / Zsig(0, i);
        }
    }
    std::cout<<"Zsig"<<std::endl;
    std::cout<<Zsig<<std::endl;

    z_pred.fill(0);
    for (int i = 0;i<2*n_aug_+1;i++)
    {
        z_pred = z_pred+weights_(i)*Zsig.col(i);
    }

    for (int i = 0;i<2*n_aug_+1;i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = weights_(i)*z_diff*z_diff.transpose();
    }



    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(n_z_,n_z_);
    R(0,0) = pow(std_radr_,2);
    R(1,1) = pow(std_radphi_,2);
    R(2,2) = pow(std_radrd_,2);
    S = S+R;

    measurementUpdate(meas_package.raw_measurements_, Zsig,  z_pred,  S);
    /*
    // update state
    Eigen::MatrixXd Tc = Eigen::MatrixXd::Zero(n_x_,n_z_);
    for(int i=0;i<2*n_aug_+1;i++)
    {
        Tc = Tc + weights_(i)*((Xsig_pred_.col(i)-x_)*((Zsig.col(i)-z_pred).transpose()));
    }

    Eigen::MatrixXd K = Tc*S.inverse();

    x_ = x_ +K*(meas_package.raw_measurements_-z_pred);
    P_ = P_-K*S*K.transpose();
    */
    //std::cout<<P_<<std::endl;
}
void UKF::generateSigmaPoints(Eigen::MatrixXd &Xsig_aug){
    Eigen::VectorXd x_aug = VectorXd::Zero(n_aug_);;
    Eigen::MatrixXd P_aug = MatrixXd::Zero(n_aug_,n_aug_);
    x_aug = VectorXd(n_aug_);



    x_aug.block(0,0,n_x_,1) = x_;
    x_aug(n_x_,0) = 0;
    x_aug(n_x_+1,0) = 0;



    P_aug.block(0,0,n_x_,n_x_) = P_;
    P_aug(n_x_,n_x_) = pow(std_a_,2);
    P_aug(n_x_+1,n_x_+1) = pow(std_yawdd_,2);

    //std::cout<<x_aug<<std::endl;
    //std::cout<<P_aug<<std::endl;


    Eigen::MatrixXd A = P_aug.llt().matrixL();


    Xsig_aug.col(0) = x_aug;
    for (int i=0; i<n_aug_;i++)
    {
        Xsig_aug.col(i+1) = x_aug+sqrt(lambda_+n_aug_)*(A.col(i));
        Xsig_aug.col(i+1+n_aug_) = x_aug-sqrt(lambda_+n_aug_)*(A.col(i));
    }
}
void UKF::predictSigmaPoints(Eigen::MatrixXd &Xsig_pred,double delta_t){
    Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd::Zero(n_aug_,2*n_aug_+1);

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
            p_x_p = p_x + v*delta_t*cos(yaw);
            p_y_p = p_y + v*delta_t*sin(yaw);
        }
        else
        {
            //p_x_p = p_x+v/d_yaw*(sin(yaw+d_yaw*delta_t)-sin(yaw));
            //p_y_p = p_y+v/d_yaw*(-cos(yaw+d_yaw*delta_t)+cos(yaw));

            p_x_p = p_x+v/d_yaw*(sin(yaw+d_yaw*delta_t)-sin(yaw));
            p_y_p = p_y+v/d_yaw*(-cos(yaw+d_yaw*delta_t)+cos(yaw));

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

void UKF::measurementUpdate(const Eigen::VectorXd& z, Eigen::MatrixXd Zsig, Eigen::VectorXd z_pred, Eigen::MatrixXd S){
    MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        VectorXd z_diff = Zsig.col(i) - z_pred;


        while (z_diff(1) > M_PI) {
            z_diff(1) = z_diff(1) - 2. * M_PI;
        }
        while (z_diff(1) <- M_PI) {
            z_diff(1) = z_diff(1) + 2. * M_PI;
        }

        // x_diff
        while (x_diff(3) > M_PI) {
            x_diff(3) = x_diff(3) - 2. * M_PI;
        }
        while (x_diff(3) <- M_PI) {
            x_diff(3) = x_diff(3) + 2. * M_PI;
        }

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;

    while (z_diff(1) > M_PI) {
        z_diff(1) = z_diff(1) - 2. * M_PI;
    }
    while (z_diff(1) <- M_PI) {
        z_diff(1) = z_diff(1) + 2. * M_PI;
    }


    P_ = P_ - K*S*K.transpose();
    x_ = x_ + K * z_diff;
}