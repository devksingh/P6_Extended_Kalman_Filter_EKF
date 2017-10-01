#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  // my changes to the old code
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  // H Jacobian (4x4)
  /*Hj_ << 1,1,0,0,
         1,1,0,0,
         1,1,1,1;*/
  //F_ matrix (4x4)
  ekf_.F_ = MatrixXd(4,4);
  // value of F_
  ekf_.F_ << 1,0,1,0,
            0,1,0,1,
            0,0,1,0,
            0,0,0,1;
  // P Matrix (4x4)
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1,0,0,0,
            0,1,0,0,
            0,0,1000,0,
            0,0,0,1000; 


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
        // get the data from measurement pack
        float roh = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        float rohdot = measurement_pack.raw_measurements_[2];
        // Convert to cartesian x= rcos(theta), y = r sin(theta) 
         float x = roh * cos(phi);
         float y = roh * sin(phi);
         float vx = rohdot * cos(phi);
         float vy = rohdot * sin(phi);
         ekf_.x_ << x,y,vx,vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
         float x = measurement_pack.raw_measurements_[0];
         float y = measurement_pack.raw_measurements_[1];
         float vx = 0;
         float vy = 0;
         ekf_.x_ << x,y,vx,vy;
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // Elapsed time dt
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt2 = pow(dt,2);
  float dt3 = pow(dt,3);
  float dt4 = pow(dt,4);


  //changing F matrix to include time elapsed
  ekf_.F_ << 1,0,dt,0,
            0,1,0,dt,
            0,0,1,0,
            0,0,0,1;

  // Setting the noise value provided by udacity
  float noise_ax = 9; 
  float noise_ay = 9;

  // Setting the Q covariance matrix or 4x4

  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << ((dt4 / 4) * noise_ax),0,((dt3 / 2) * noise_ax),0,
             0,((dt4 / 4) * noise_ay),0,((dt3 / 2) * noise_ay),
             ((dt3 / 2) * noise_ax),0,(dt2 * noise_ax),0,
             0,((dt3 / 2) * noise_ay),0,(dt2 * noise_ay);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tools;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
